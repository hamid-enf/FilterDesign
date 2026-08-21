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

/* ------------------------------------------------------------------ */
/* ideal impulse response                                              */
/* ------------------------------------------------------------------ */

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
            v = -(2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            if (m == 0.0)
                v += 1.0;
            break;

        case FCE_FIR_BANDPASS:
            v = (2.0 * fc2 / fs) * fce_sinc(2.0 * fc2 * m / fs)
              - (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs);
            break;

        case FCE_FIR_BANDSTOP:
            v = (2.0 * fc1 / fs) * fce_sinc(2.0 * fc1 * m / fs)
              - (2.0 * fc2 / fs) * fce_sinc(2.0 * fc2 * m / fs);
            if (m == 0.0)
                v += 1.0;
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
    uint32_t grid = 256u;
    double f_lo_c = (f_lo <= 0.0) ? 0.0 : f_lo;
    double f_hi_c = (f_hi >= 0.5 * fs) ? 0.5 * fs : f_hi;
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
        window[n] = fce_window_value(sp->window, n, N, beta, 0.5);

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
        switch (sp->fir_type)
        {
        case FCE_FIR_BANDPASS:
        case FCE_FIR_BANDSTOP:
            f_lo = sp->fc1; f_hi = sp->fc2; break;
        case FCE_FIR_HILBERT:
            f_lo = (sp->fc1 > 0.0 && sp->fc2 > sp->fc1) ? sp->fc1 : 0.0;
            f_hi = (sp->fc1 > 0.0 && sp->fc2 > sp->fc1) ? sp->fc2
                                                        : 0.5 * sp->fs;
            break;
        default:
            f_lo = 0.0; f_hi = 0.5 * sp->fs; break;
        }
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

    if (!(norm_gain > 0.0) || !(norm_gain < 1e300))
        return FCE_ERR_NUMERICAL;

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
        r->symmetry = FCE_SYMMETRY_III;
        break;
    case FCE_FIR_DIFFERENTIATOR:
        r->symmetry = FCE_SYMMETRY_IV;
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
