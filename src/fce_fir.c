/*
 * fce_fir.c - FIR coefficient generation (windowed-sinc method).
 *
 * Pipeline (all internal math in float64):
 *   ideal impulse response -> window -> normalize -> float32/fixed-point
 *
 * Ideal responses follow the SciPy firwin conventions:
 *   LP:  h[n] = (2 fc/fs) sinc(2 fc (n-M)/fs)
 *   HP:  h[n] = delta(n-M) - LP
 *   BP:  h[n] = LP(fc2) - LP(fc1)
 *   BS:  h[n] = delta(n-M) - BP
 *   Hilbert:      h[n] = (cos(w1 m) - cos(w2 m)) / (pi m),  m = n-M
 *   Differentiator: h[n] = -cos(pi m)/m  (full band, H(e^jw) ~= jw)
 */
#include "fce_internal.h"

/* ------------------------------------------------------------------ */
/* windows                                                             */
/* ------------------------------------------------------------------ */

/* Tukey shape parameter: the spec carries no alpha field, so a fixed
 * alpha (half-cosine taper over the outer quarters) is used. */
#define FCE_TUKEY_DEFAULT_ALPHA 0.5

double fce_window_value(fce_window_t win, uint32_t n, uint32_t N,
                        double kaiser_beta, double tukey_alpha)
{
    double x;

    switch (win)
    {
    case FCE_WIN_RECTANGULAR:
        return 1.0;

    case FCE_WIN_HANN:
        return 0.5 * (1.0 - cos(2.0 * FCE_PI * (double)n / (double)(N - 1)));

    case FCE_WIN_HAMMING:
        return 0.54 - 0.46 * cos(2.0 * FCE_PI * (double)n / (double)(N - 1));

    case FCE_WIN_BLACKMAN:
        return 0.42 - 0.5 * cos(2.0 * FCE_PI * (double)n / (double)(N - 1))
                   + 0.08 * cos(4.0 * FCE_PI * (double)n / (double)(N - 1));

    case FCE_WIN_KAISER:
    {
        double t = 2.0 * (double)n / (double)(N - 1) - 1.0;
        return fce_i0(kaiser_beta * sqrt(1.0 - t * t)) / fce_i0(kaiser_beta);
    }

    case FCE_WIN_BLACKMAN_HARRIS:
        return 0.35875
             - 0.48829 * cos(2.0 * FCE_PI * (double)n / (double)(N - 1))
             + 0.14128 * cos(4.0 * FCE_PI * (double)n / (double)(N - 1))
             - 0.01168 * cos(6.0 * FCE_PI * (double)n / (double)(N - 1));

    case FCE_WIN_BARTLETT:
        x = 2.0 * (double)n / (double)(N - 1) - 1.0;
        return 1.0 - fabs(x);

    case FCE_WIN_TUKEY:
    {
        /* scipy tukey(alpha), sym */
        double alpha = tukey_alpha;
        double nn = (double)n;
        double nmax = (double)(N - 1);
        if (alpha <= 0.0)
            return 1.0;
        if (alpha >= 1.0)
            return 0.5 * (1.0 - cos(2.0 * FCE_PI * nn / nmax));
        if (nn < alpha * nmax * 0.5)
            return 0.5 * (1.0 + cos(FCE_PI * (2.0 * nn / (alpha * nmax) - 1.0)));
        if (nn >= nmax * (1.0 - alpha * 0.5))
            return 0.5 * (1.0 + cos(FCE_PI * (2.0 * nn / (alpha * nmax)
                                             - 2.0 / alpha + 1.0)));
        return 1.0;
    }

    default:
        return 1.0;
    }
}

/*
 * Half-width of the transition band of a windowed design, in Hz: the
 * distance from the cutoff (the -6 dB midpoint) to the edge of the
 * window's main lobe (first null of the window's DTFT). The windowed
 * response is flat within its sidelobe floor only outside that region,
 * so the validation insets the passband ripple band by this amount.
 *
 * Fixed windows (first null in units of fs/M, verified against the
 * exact window DTFTs):
 *   rectangular 1, bartlett 2, hann 2, hamming 2, blackman 3,
 *   blackman-harris 4, tukey (fixed alpha = FCE_TUKEY_DEFAULT_ALPHA) 4
 * Kaiser: half of the design transition width
 *   dw = (A - 7.95) / (2.285 * (M - 1))   (same relation as
 *   fce_kaiser_taps), i.e. dw/2 in Hz.
 */
double fce_window_transition_half_hz(fce_window_t win, double atten_db,
                                     double kaiser_beta,
                                     uint32_t num_taps, double fs)
{
    double k;
    if (num_taps < 2u || fs <= 0.0)
        return 0.0;

    switch (win)
    {
    case FCE_WIN_RECTANGULAR:
        k = 1.0;
        break;

    case FCE_WIN_BARTLETT:
    case FCE_WIN_HANN:
    case FCE_WIN_HAMMING:
        k = 2.0;
        break;

    case FCE_WIN_BLACKMAN:
        k = 3.0;
        break;

    case FCE_WIN_BLACKMAN_HARRIS:
    case FCE_WIN_TUKEY:
        k = 4.0;
        break;

    case FCE_WIN_KAISER:
    {
        double a = (atten_db > 0.0) ? atten_db : 0.0;
        if (a <= 0.0 && kaiser_beta > 0.5)
            a = kaiser_beta / 0.1102 + 8.7; /* invert fce_kaiser_beta */
        if (a > 21.0)
        {
            double dw = (a - 7.95) / (2.285 * (double)(num_taps - 1u));
            return dw * 0.5 * fs / (2.0 * FCE_PI);
        }
        k = 2.0; /* unknown attenuation: Hann-like default */
        break;
    }

    default:
        k = 2.0;
        break;
    }
    return k * fs / (double)num_taps;
}

static void fce_fir_ideal(fce_fir_type_t type,
                          double fs, double fc1, double fc2,
                          uint32_t N, double* h)
{
    double M = 0.5 * (double)(N - 1);
    uint32_t n;

    for (n = 0; n < N; n++)
    {
        double m = (double)n - M;
        double v = 0.0;

        switch (type)
        {
        case FCE_FIR_LOWPASS:
            v = (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            break;

        case FCE_FIR_HIGHPASS:
            /* delta - LP. The band-limited delta is sinc(m) (= sin(pi m)/(pi m)),
             * NOT "1 at m == 0 only": for even N the samples sit at
             * half-integer m, where sinc(m) != 0. Writing the delta as
             * `if (m == 0) v += 1` silently drops it for all even tap
             * counts and turns the highpass into a negated lowpass. */
            v = fce_sinc(m) - (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            break;

        case FCE_FIR_BANDPASS:
            v = (2.0 * fc2 / fs) * fce_sinc(2.0 * fc2 * m / fs)
              - (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            break;

        case FCE_FIR_BANDSTOP:
            /* delta - BP: same sinc(m) band-limited delta as the HP above */
            v = fce_sinc(m)
              - (2.0 * fc2 / fs) * fce_sinc(2.0 * fc2 * m / fs)
              + (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            break;

        case FCE_FIR_HILBERT:
        {
            double w1 = (fc1 > 0.0 && fc2 > fc1) ? 2.0 * FCE_PI * fc1 / fs : 0.0;
            double w2 = (fc1 > 0.0 && fc2 > fc1) ? 2.0 * FCE_PI * fc2 / fs
                                                 : FCE_PI;
            if (m == 0.0)
                v = 0.0;
            else
                v = (cos(w1 * m) - cos(w2 * m)) / (FCE_PI * m);
            break;
        }

        case FCE_FIR_DIFFERENTIATOR:
        {
            double w1 = (fc1 > 0.0 && fc2 > fc1) ? 2.0 * FCE_PI * fc1 / fs : 0.0;
            double w2 = (fc1 > 0.0 && fc2 > fc1) ? 2.0 * FCE_PI * fc2 / fs
                                                 : FCE_PI;
            if (m == 0.0)
                v = 0.0;
            else
            {
                /* band-limited differentiator, H ~= jw on [w1, w2] */
                v = (w2 * cos(w2 * m) - w1 * cos(w1 * m)
                     - (sin(w2 * m) - sin(w1 * m)) / m) / (FCE_PI * m);
            }
            break;
        }

        default:
            v = 0.0;
            break;
        }

        h[n] = v;
    }
}

/* ------------------------------------------------------------------ */
/* normalization reference frequency helper                            */
/* ------------------------------------------------------------------ */

static double fce_fir_gain_at(const double* h, uint32_t N, double w)
{
    /* H(e^jw) = sum h[n] e^{-jwn} ; real part for w in {0, pi} */
    double re = 0.0, im = 0.0;
    uint32_t n;
    for (n = 0; n < N; n++)
    {
        double a = w * (double)n;
        re += h[n] * cos(a);
        im -= h[n] * sin(a);
    }
    return sqrt(re * re + im * im);
}

static double fce_fir_peak_in_band(const double* h, uint32_t N,
                                   double fs, double f_lo, double f_hi)
{
    double best = 0.0;
    uint32_t i;
    /* The golden-section refinement assumes a unimodal bracket. The
     * response ripples with period ~ fs / (N-1) Hz, so a fixed sparse
     * grid can put several ripple peaks into one bracket and the
     * refinement then lands on an arbitrary sub-peak (observed peak
     * estimate error ~5e-5 for N ~ 100). Scale the grid with the tap
     * count so one bracket always spans less than a ripple period. */
    uint32_t grid = 4u * N;
    double f_lo_c = (f_lo <= 0.0) ? 0.0 : f_lo;
    double f_hi_c = (f_hi >= 0.5 * fs) ? 0.5 * fs : f_hi;
    if (grid < 256u)
        grid = 256u;
    if (grid > 4096u)
        grid = 4096u;
    if (f_hi_c <= f_lo_c)
        f_hi_c = 0.5 * fs;

    double f_best = 0.0;
    for (i = 0; i <= grid; i++)
    {
        double f = f_lo_c + (f_hi_c - f_lo_c) * (double)i / (double)grid;
        double g = fce_fir_gain_at(h, N, 2.0 * FCE_PI * f / fs);
        if (g > best)
        {
            best = g;
            f_best = f;
        }
    }

    /* golden-section refinement around the coarse maximum */
    {
        const double gr = 0.6180339887498948482;
        double a = f_best - (f_hi_c - f_lo_c) / (double)grid;
        double b = f_best + (f_hi_c - f_lo_c) / (double)grid;
        double c, d;
        double fc, fd;
        if (a < f_lo_c)
            a = f_lo_c;
        if (b > f_hi_c)
            b = f_hi_c;
        c = b - gr * (b - a);
        d = a + gr * (b - a);
        fc = fce_fir_gain_at(h, N, 2.0 * FCE_PI * c / fs);
        fd = fce_fir_gain_at(h, N, 2.0 * FCE_PI * d / fs);
        for (i = 0; i < 60u; i++)
        {
            if (fc > fd)
            {
                b = d;
                d = c;
                fd = fc;
                c = b - gr * (b - a);
                fc = fce_fir_gain_at(h, N, 2.0 * FCE_PI * c / fs);
            }
            else
            {
                a = c;
                c = d;
                fc = fd;
                d = a + gr * (b - a);
                fd = fce_fir_gain_at(h, N, 2.0 * FCE_PI * d / fs);
            }
            if (b - a < 1e-12 * (1.0 + fabs(b)))
                break;
        }
        best = 0.5 * (fc + fd);
    }
    return best;
}

/* passband interval [f_lo, f_hi] used for passband-peak normalization */
static void fce_fir_passband(const fce_spec_t* sp, double* f_lo, double* f_hi)
{
    switch (sp->fir_type)
    {
    case FCE_FIR_BANDPASS:
        *f_lo = sp->fc1;
        *f_hi = sp->fc2;
        break;
    case FCE_FIR_HIGHPASS:
        *f_lo = sp->fc1;
        *f_hi = 0.5 * sp->fs;
        break;
    case FCE_FIR_HILBERT:
        if (sp->fc1 > 0.0 && sp->fc2 > sp->fc1)
        {
            *f_lo = sp->fc1;
            *f_hi = sp->fc2;
        }
        else
        {
            *f_lo = 0.0;
            *f_hi = 0.5 * sp->fs;
        }
        break;
    default:
        /* LP / BS / differentiator: scan the full band. For BS the old
         * code scanned [fc1, fc2] which is the STOP band - normalizing
         * the stopband peak to 1 amplified the taps into nonsense. */
        *f_lo = 0.0;
        *f_hi = 0.5 * sp->fs;
        break;
    }
}

/* is a normalization reference gain usable, relative to the tap scale? */
static int fce_norm_gain_ok(double norm_gain, double hmax)
{
    return (norm_gain > 0.0) && (norm_gain < 1e300) &&
           (norm_gain > 1e-12 * hmax);
}

/* ------------------------------------------------------------------ */
/* design entry                                                        */
/* ------------------------------------------------------------------ */

fce_status_t fce_fir_design(const fce_spec_t* sp, fce_result_t* r,
                            fce_layout_t* lay, void* base)
{
    double* work = (double*)(void*)((char*)base + lay->off_fir_work);
    double* h    = (double*)(void*)((char*)base + lay->off_fir_h);
    float*  h32  = (float*)(void*)((char*)base + lay->off_fir_h32);
    double* ideal = work;            /* [0..N)   */
    double* window = work + lay->n_taps; /* [N..2N) */
    uint32_t N = lay->n_taps;
    double beta = 0.0;
    uint32_t n;
    double norm_gain = 1.0;
    fce_norm_t norm = sp->normalization;

    /* ---- Kaiser beta ---- */
    if (sp->window == FCE_WIN_KAISER)
    {
        if (sp->kaiser_beta > 0.0)
            beta = sp->kaiser_beta;
        else if (sp->stopband_atten_db > 0.0)
            beta = fce_kaiser_beta(sp->stopband_atten_db);
        else
            beta = 0.0;
    }
    r->kaiser_beta = beta;

    /* ---- ideal response ---- */
    fce_fir_ideal(sp->fir_type, sp->fs, sp->fc1, sp->fc2, N, ideal);

    /* ---- window ---- */
    for (n = 0; n < N; n++)
        window[n] = fce_window_value(sp->window, n, N, beta,
                                     FCE_TUKEY_DEFAULT_ALPHA);

    for (n = 0; n < N; n++)
        h[n] = ideal[n] * window[n];

    /* ---- normalization ---- */
    if (norm == FCE_NORM_AUTO)
    {
        switch (sp->fir_type)
        {
        case FCE_FIR_LOWPASS:    norm = FCE_NORM_DC; break;
        case FCE_FIR_HIGHPASS:   norm = FCE_NORM_NYQUIST; break;
        case FCE_FIR_BANDPASS:   norm = FCE_NORM_PASSBAND_PEAK; break;
        case FCE_FIR_BANDSTOP:   norm = FCE_NORM_DC; break;
        case FCE_FIR_HILBERT:    norm = FCE_NORM_PASSBAND_PEAK; break;
        case FCE_FIR_DIFFERENTIATOR: norm = FCE_NORM_NYQUIST; break;
        default:                 norm = FCE_NORM_DC; break;
        }
    }

    switch (norm)
    {
    case FCE_NORM_DC:
        norm_gain = fce_fir_gain_at(h, N, 0.0);
        break;
    case FCE_NORM_NYQUIST:
        norm_gain = fce_fir_gain_at(h, N, FCE_PI);
        break;
    case FCE_NORM_PASSBAND_PEAK:
    {
        double f_lo, f_hi;
        fce_fir_passband(sp, &f_lo, &f_hi);
        norm_gain = fce_fir_peak_in_band(h, N, sp->fs, f_lo, f_hi);
        break;
    }
    case FCE_NORM_NONE:
        norm_gain = 1.0;
        break;
    default:
        norm_gain = 1.0;
        break;
    }

    /* guard degenerate normalization references: a Nyquist/DC edge
     * normalization on a filter with (numerically) zero gain at that
     * edge would otherwise silently amplify the taps by ~1e16.
     *
     * When this happens with FCE_NORM_AUTO it is a symmetry-type null
     * landing exactly on the AUTO reference (an even-tap highpass always
     * has a Nyquist null, an odd-tap differentiator always has one too).
     * AUTO then falls back to the passband peak, which is well defined
     * for every symmetry type, and raises FCE_FLAG_SYMMETRY_WARNING.
     * An explicitly requested degenerate reference is a spec error and
     * still returns FCE_ERR_NUMERICAL. */
    {
        double hmax = 0.0;
        for (n = 0; n < N; n++)
        {
            double a = fabs(h[n]);
            if (a > hmax)
                hmax = a;
        }
        if (!fce_norm_gain_ok(norm_gain, hmax) &&
            sp->normalization == FCE_NORM_AUTO &&
            norm != FCE_NORM_PASSBAND_PEAK)
        {
            double f_lo, f_hi;
            norm = FCE_NORM_PASSBAND_PEAK;
            fce_fir_passband(sp, &f_lo, &f_hi);
            norm_gain = fce_fir_peak_in_band(h, N, sp->fs, f_lo, f_hi);
            r->flags |= FCE_FLAG_SYMMETRY_WARNING;
        }
        if (!fce_norm_gain_ok(norm_gain, hmax))
            return FCE_ERR_NUMERICAL;

        /* amplification beyond ~1e8x (an explicitly requested reference
         * in a deep stopband) is honored but flagged: the gain is nearly
         * all numerical noise at that point */
        if (!(norm_gain > 1e-8 * hmax))
            r->flags |= FCE_FLAG_NUMERICAL_WARNING;
    }

    for (n = 0; n < N; n++)
        h[n] /= norm_gain;
    r->norm_factor = 1.0 / norm_gain;
    r->normalization = norm;

    /* ---- float32 output (float64 always available: internal math is
     * float64; float32 is only the final rounding) ---- */
    r->h_f64 = h;
    if (sp->precision == FCE_PRECISION_FLOAT32)
    {
        for (n = 0; n < N; n++)
            h32[n] = (float)h[n];
        r->h_f32 = h32;
    }
    else
    {
        r->h_f32 = NULL;
    }

    /* ---- symmetry info ---- */
    switch (sp->fir_type)
    {
    case FCE_FIR_LOWPASS:
    case FCE_FIR_HIGHPASS:
    case FCE_FIR_BANDPASS:
    case FCE_FIR_BANDSTOP:
        r->symmetry = (N & 1u) ? FCE_SYMMETRY_I : FCE_SYMMETRY_II;
        if (!(N & 1u) && (sp->fir_type == FCE_FIR_HIGHPASS ||
                          sp->fir_type == FCE_FIR_BANDSTOP))
            r->flags |= FCE_FLAG_SYMMETRY_WARNING; /* Type II has a Nyquist null */
        break;
    case FCE_FIR_HILBERT:
    case FCE_FIR_DIFFERENTIATOR:
        /* both are anti-symmetric: Type III (odd taps) or Type IV (even) */
        r->symmetry = (N & 1u) ? FCE_SYMMETRY_III : FCE_SYMMETRY_IV;
        break;
    default:
        r->symmetry = FCE_SYMMETRY_NONE;
        break;
    }

    r->num_taps = (uint16_t)N;

    /* ---- expose internals (no black box) ---- */
    r->fir_ideal = ideal;
    r->fir_window = window;

    return FCE_OK;
}
