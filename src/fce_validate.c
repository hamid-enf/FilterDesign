/*
 * fce_validate.c - validation engine.
 *
 * Frequency response scans, stability analysis, spec-vs-measured checks
 * and coefficient error metrics. This module exists to prove that the
 * generated coefficients are correct - it is not a filter runtime.
 */
#include "fce_internal.h"

#if FCE_ENABLE_VALIDATION

/* ================================================================== */
/* response evaluation primitives                                      */
/* ================================================================== */

static fce_cplx_t fce_eval_fir_c(const double* h, uint32_t n, double w)
{
    /* H(e^jw) = sum h[n] e^{-jwn} */
    double re = 0.0, im = 0.0;
    double cw = cos(w), sw = sin(w);
    double zr = 1.0, zi = 0.0; /* z = e^{-jw}, z^k */
    uint32_t i;
    for (i = 0; i < n; i++)
    {
        re += h[i] * zr;
        im += h[i] * zi;
        {
            double nr = zr * cw - zi * (-sw);
            double ni = zr * (-sw) + zi * cw;
            zr = nr;
            zi = ni;
        }
    }
    return fce_cx(re, im);
}

static fce_cplx_t fce_eval_sos_c(const double* sos, uint32_t ns, double w)
{
    fce_cplx_t h = fce_cx(1.0, 0.0);
    uint32_t i;
    for (i = 0; i < ns; i++)
        h = fce_cx_mul(h, fce_eval_biquad(sos + 5u * i, w));
    return h;
}

static void fce_fill_point(fce_response_point_t* pt, double f, fce_cplx_t h)
{
    double mag = fce_cx_abs(h);
    pt->f_hz = f;
    pt->mag = mag;
    pt->mag_db = 20.0 * log10(mag > 0.0 ? mag : 1e-300);
    pt->phase_deg = atan2(h.im, h.re) * 180.0 / FCE_PI;
    pt->group_delay = 0.0;
}

fce_status_t fce_response_fir(const double* h, uint16_t n,
                              double fs, uint32_t n_points,
                              double f_start, double f_stop,
                              fce_response_cb cb, void* ctx)
{
    uint32_t i;
    if (h == NULL || cb == NULL || n == 0u || fs <= 0.0)
        return FCE_ERR_INVALID_ARGUMENT;
    if (n_points == 0u)
        n_points = 1u;
    /* f_stop == f_start probes that single frequency;
     * only an inverted range falls back to the full band */
    if (f_stop < f_start)
    {
        f_start = 0.0;
        f_stop = 0.5 * fs;
    }

    for (i = 0; i < n_points; i++)
    {
        double t = (n_points > 1u) ? (double)i / (double)(n_points - 1u) : 0.0;
        double f = f_start + (f_stop - f_start) * t;
        double w = 2.0 * FCE_PI * f / fs;
        fce_cplx_t hc = fce_eval_fir_c(h, n, w);
        fce_response_point_t pt;
        fce_fill_point(&pt, f, hc);
        /* analytic group delay: GD = -Im(H' conj(H)) / |H|^2 */
        {
            double cw = cos(w), sw = sin(w);
            double zr = 1.0, zi = 0.0;
            double dr = 0.0, di = 0.0;
            uint32_t k;
            for (k = 0; k < n; k++)
            {
                /* d/dw of z^k with z = e^{-jw}: -j k z^k */
                dr += h[k] * ((double)k * zi);
                di += h[k] * (-(double)k * zr);
                {
                    double nr = zr * cw - zi * (-sw);
                    double ni = zr * (-sw) + zi * cw;
                    zr = nr;
                    zi = ni;
                }
            }
            /* H' = dr + j di ;  GD = -Im(H' conj(H))/|H|^2
             * Im(H' conj(H)) = di*re - dr*im... check:
             * conj(H) = re - j im ; H' conj(H) = (dr+j di)(re - j im)
             * im = di*re - dr*im */
            {
                double m2 = hc.re * hc.re + hc.im * hc.im;
                /* at (or within numerical floor of) a response zero the
                 * ratio is 0/0 and the analytic value is undefined
                 * (e.g. DC of a highpass): report 0 instead of noise */
                if (m2 > 1e-24)
                    pt.group_delay = -(di * hc.re - dr * hc.im) / m2;
            }
        }
        if (!cb(ctx, &pt))
            break;
    }
    return FCE_OK;
}

fce_status_t fce_response_sos(const double* sos, uint16_t n_sections,
                              double fs, uint32_t n_points,
                              double f_start, double f_stop,
                              fce_response_cb cb, void* ctx)
{
    uint32_t i;
    if (sos == NULL || cb == NULL || n_sections == 0u || fs <= 0.0)
        return FCE_ERR_INVALID_ARGUMENT;
    if (n_points == 0u)
        n_points = 1u;
    /* f_stop == f_start probes that single frequency;
     * only an inverted range falls back to the full band */
    if (f_stop < f_start)
    {
        f_start = 0.0;
        f_stop = 0.5 * fs;
    }

    for (i = 0; i < n_points; i++)
    {
        double t = (n_points > 1u) ? (double)i / (double)(n_points - 1u) : 0.0;
        double f = f_start + (f_stop - f_start) * t;
        double w = 2.0 * FCE_PI * f / fs;
        fce_cplx_t hc = fce_eval_sos_c(sos, n_sections, w);
        fce_response_point_t pt;
        uint32_t s;
        fce_fill_point(&pt, f, hc);
        pt.group_delay = 0.0;
        for (s = 0; s < n_sections; s++)
            pt.group_delay += fce_biquad_group_delay(sos + 5u * s, w);
        /* at (or within numerical floor of) a response zero the summed
         * section delays are 0/0 noise (e.g. DC of a highpass whose
         * numerator vanishes): report 0 instead */
        {
            double m2 = hc.re * hc.re + hc.im * hc.im;
            if (m2 <= 1e-24)
                pt.group_delay = 0.0;
        }
        if (!cb(ctx, &pt))
            break;
    }
    return FCE_OK;
}

/* ================================================================== */
/* stability                                                           */
/* ================================================================== */

fce_status_t fce_stability_sos(const double* sos, uint16_t n_sections,
                               double* max_radius, double* margin)
{
    double mr = 0.0;
    uint32_t i;
    if (sos == NULL || n_sections == 0u)
        return FCE_ERR_INVALID_ARGUMENT;

    for (i = 0; i < n_sections; i++)
    {
        fce_cplx_t p1, p2;
        double r1, r2;
        fce_biquad_poles(sos[5u * i + 3], sos[5u * i + 4], &p1, &p2);
        r1 = fce_cx_abs(p1);
        r2 = fce_cx_abs(p2);
        if (r1 > mr)
            mr = r1;
        if (r2 > mr)
            mr = r2;
    }

    if (max_radius)
        *max_radius = mr;
    if (margin)
        *margin = 1.0 - mr;

    /* poles within the numerical noise band [1, 1 + tol] of extreme
     * designs are reported via margin (< 0) but do not hard-fail:
     * the float64 extraction itself cannot resolve them better, and
     * the equivalent scipy design is considered stable */
    if (mr >= 1.0 + FCE_STABILITY_RADIUS_TOL)
        return FCE_ERR_UNSTABLE;
    return FCE_OK;
}

/* ================================================================== */
/* coefficient error metrics                                           */
/* ================================================================== */

void fce_coeff_error(const double* ref, const double* test, uint32_t n,
                     double* max_abs, double* rms, double* max_rel)
{
    uint32_t i;
    double mx = 0.0, sse = 0.0, mr = 0.0;
    if (ref == NULL || test == NULL || n == 0u)
    {
        if (max_abs) *max_abs = 0.0;
        if (rms) *rms = 0.0;
        if (max_rel) *max_rel = 0.0;
        return;
    }
    for (i = 0; i < n; i++)
    {
        double e = fabs(ref[i] - test[i]);
        if (e > mx)
            mx = e;
        sse += e * e;
        if (ref[i] != 0.0)
        {
            double rel = e / fabs(ref[i]);
            if (rel > mr)
                mr = rel;
        }
    }
    if (max_abs) *max_abs = mx;
    if (rms) *rms = sqrt(sse / (double)n);
    if (max_rel) *max_rel = mr;
}

/* ================================================================== */
/* spec checks used by fce_generate                                    */
/* ================================================================== */

/* scan |H| of an SOS filter on a grid; collect min/max in a band */
void fce_scan_sos_band(const double* sos, uint32_t ns, double fs,
                              double f_lo, double f_hi,
                              double* min_db, double* max_db)
{
    uint32_t i;
    uint32_t grid = FCE_VALIDATE_GRID_POINTS;
    double mn = 1e300, mx = -1e300;
    if (f_hi <= f_lo)
    {
        f_lo = 0.0;
        f_hi = 0.5 * fs;
    }
    for (i = 0; i <= grid; i++)
    {
        double f = f_lo + (f_hi - f_lo) * (double)i / (double)grid;
        double w = 2.0 * FCE_PI * f / fs;
        double g = fce_cx_abs(fce_eval_sos_c(sos, ns, w));
        double db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        if (db < mn) mn = db;
        if (db > mx) mx = db;
    }
    *min_db = mn;
    *max_db = mx;
}

void fce_scan_fir_band(const double* h, uint32_t n, double fs,
                              double f_lo, double f_hi,
                              double* min_db, double* max_db)
{
    uint32_t i;
    uint32_t grid = FCE_VALIDATE_GRID_POINTS;
    double mn = 1e300, mx = -1e300;
    if (f_hi <= f_lo)
    {
        f_lo = 0.0;
        f_hi = 0.5 * fs;
    }
    for (i = 0; i <= grid; i++)
    {
        double f = f_lo + (f_hi - f_lo) * (double)i / (double)grid;
        double w = 2.0 * FCE_PI * f / fs;
        double g = fce_cx_abs(fce_eval_fir_c(h, n, w));
        double db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        if (db < mn) mn = db;
        if (db > mx) mx = db;
    }
    *min_db = mn;
    *max_db = mx;
}

/* find the first frequency where the gain crosses `target_db`, scanning
 * from f_lo towards f_hi.
 *   mode 0: down-crossing  (gain above target, then drops below)
 *   mode 1: up-crossing    (gain below target, then rises above)
 * Returns f_hi/f_lo if no crossing is found in the scan. */
static double fce_find_crossing_sos(const double* sos, uint32_t ns, double fs,
                                    double f_lo, double f_hi,
                                    double target_db, int mode)
{
    uint32_t i;
    uint32_t grid = FCE_VALIDATE_GRID_POINTS;
    double fa = f_lo, fb = f_hi;
    int found = 0;
    double prev_db = 0.0;

    for (i = 0; i <= grid; i++)
    {
        double f = f_lo + (f_hi - f_lo) * (double)i / (double)grid;
        double w = 2.0 * FCE_PI * f / fs;
        double g = fce_cx_abs(fce_eval_sos_c(sos, ns, w));
        double db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        if (i > 0)
        {
            if (mode == 0 && prev_db > target_db && db <= target_db)
            {
                fa = f_lo + (f_hi - f_lo) * (double)(i - 1u) / (double)grid;
                fb = f;
                found = 1;
                break;
            }
            if (mode == 1 && prev_db < target_db && db >= target_db)
            {
                fa = f_lo + (f_hi - f_lo) * (double)(i - 1u) / (double)grid;
                fb = f;
                found = 1;
                break;
            }
        }
        prev_db = db;
    }
    if (!found)
        return (mode == 0) ? f_hi : f_lo;

    for (i = 0; i < 40u; i++)
    {
        double fm = 0.5 * (fa + fb);
        double w = 2.0 * FCE_PI * fm / fs;
        double g = fce_cx_abs(fce_eval_sos_c(sos, ns, w));
        double db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        if ((mode == 0 && db <= target_db) || (mode == 1 && db >= target_db))
            fb = fm;
        else
            fa = fm;
    }
    return 0.5 * (fa + fb);
}

/* measure passband ripple, stopband attenuation and cutoff for an IIR */
void fce_validate_iir_measures(const fce_spec_t* sp, fce_result_t* r,
                               const double* sos, uint32_t ns)
{
    double fs = sp->fs;
    double nyq = 0.5 * fs;
    double pb_lo, pb_hi, sb_lo, sb_hi;
    double mn, mx;

    if (sp->iir_type == FCE_IIR_LOWPASS)
    {
        pb_lo = 0.0; pb_hi = sp->fc1;
        sb_lo = sp->fc1 + 0.5 * (nyq - sp->fc1); sb_hi = nyq;
    }
    else if (sp->iir_type == FCE_IIR_HIGHPASS)
    {
        pb_lo = sp->fc1; pb_hi = nyq;
        sb_lo = 0.0; sb_hi = 0.5 * sp->fc1;
    }
    else if (sp->iir_type == FCE_IIR_BANDPASS)
    {
        pb_lo = sp->fc1; pb_hi = sp->fc2;
        sb_lo = 0.0; sb_hi = 0.5 * sp->fc1;
    }
    else
    {
        pb_lo = 0.0; pb_hi = sp->fc1;
        sb_lo = sp->fc1; sb_hi = sp->fc2;
    }

    fce_scan_sos_band(sos, ns, fs, pb_lo, pb_hi, &mn, &mx);
    r->passband_ripple_measured_db = mx - mn;

    if (sp->iir_type == FCE_IIR_BANDPASS)
    {
        /* stopband below AND above the passband; take the worse case.
         * (The passband scan above leaves mx holding the PASSBAND max,
         *  so the lower stopband must be scanned separately.) */
        double mn1, mx1, mn2, mx2, mxw;
        fce_scan_sos_band(sos, ns, fs, sb_lo, sb_hi, &mn1, &mx1);
        fce_scan_sos_band(sos, ns, fs, sp->fc2 + 0.5 * (nyq - sp->fc2),
                          nyq, &mn2, &mx2);
        mxw = (mx1 > mx2) ? mx1 : mx2;
        r->stopband_atten_measured_db = -mxw;
    }
    else if (sp->iir_type == FCE_IIR_BANDSTOP)
    {
        /* ripple in the upper passband too; take the worse case */
        double mn2, mx2, range;
        fce_scan_sos_band(sos, ns, fs, sp->fc2, nyq, &mn2, &mx2);
        range = (mx2 - mn2) > (mx - mn) ? (mx2 - mn2) : (mx - mn);
        r->passband_ripple_measured_db = range;
        fce_scan_sos_band(sos, ns, fs, sb_lo, sb_hi, &mn, &mx);
        r->stopband_atten_measured_db = -mx;
    }
    else
    {
        fce_scan_sos_band(sos, ns, fs, sb_lo, sb_hi, &mn, &mx);
        r->stopband_atten_measured_db = -mx;
    }

    /* cutoff crossing */
    {
        double target = -3.0;
        if (sp->iir_family == FCE_IIR_CHEBYSHEV1 ||
            sp->iir_family == FCE_IIR_ELLIPTIC)
            target = -(sp->passband_ripple_db > 0.0 ? sp->passband_ripple_db
                                                    : 3.0);
        if (sp->iir_type == FCE_IIR_LOWPASS)
            r->cutoff_measured_hz = fce_find_crossing_sos(sos, ns, fs, 0.0,
                                                          nyq, target, 0);
        else if (sp->iir_type == FCE_IIR_HIGHPASS)
            r->cutoff_measured_hz = fce_find_crossing_sos(sos, ns, fs, 0.0,
                                                          nyq, target, 1);
        else if (sp->iir_type == FCE_IIR_BANDPASS)
            r->cutoff_measured_hz = fce_find_crossing_sos(sos, ns, fs, 0.0,
                                                          sp->fc2, target, 1);
        else
            r->cutoff_measured_hz = fce_find_crossing_sos(sos, ns, fs, 0.0,
                                                          sp->fc2, target, 0);
    }

    /* DC and Nyquist gains */
    {
        double g;
        fce_cplx_t h = fce_eval_sos_c(sos, ns, 0.0);
        g = fce_cx_abs(h);
        r->dc_gain_db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        h = fce_eval_sos_c(sos, ns, FCE_PI);
        g = fce_cx_abs(h);
        r->nyquist_gain_db = 20.0 * log10(g > 0.0 ? g : 1e-300);
    }
}

void fce_validate_fir_measures(const fce_spec_t* sp, fce_result_t* r,
                               const double* h, uint32_t n)
{
    double fs = sp->fs;
    double nyq = 0.5 * fs;
    double pb_lo, pb_hi, sb_lo, sb_hi;
    double mn, mx;
    /* transition half-width: the windowed response is only flat within
     * its sidelobe floor this far away from the cutoff (the -6 dB
     * midpoint).  The passband ripple band is inset by it so the
     * metric reflects the flat passband, not the transition (which
     * would always read ~6 dB). */
    double dtrans = fce_window_transition_half_hz(sp->window,
                                                  sp->stopband_atten_db,
                                                  r->kaiser_beta, n, fs);

    if (sp->fir_type == FCE_FIR_LOWPASS)
    {
        pb_lo = 0.0; pb_hi = sp->fc1 - dtrans;
        if (pb_hi < 0.0)
            pb_hi = sp->fc1; /* degenerate spec: keep the full band */
        sb_lo = sp->fc1 + 0.5 * (nyq - sp->fc1); sb_hi = nyq;
    }
    else if (sp->fir_type == FCE_FIR_HIGHPASS)
    {
        pb_lo = sp->fc1 + dtrans;
        if (pb_lo > nyq)
            pb_lo = sp->fc1;
        pb_hi = nyq;
        sb_lo = 0.0; sb_hi = 0.5 * sp->fc1;
    }
    else if (sp->fir_type == FCE_FIR_BANDPASS)
    {
        pb_lo = sp->fc1 + dtrans; pb_hi = sp->fc2 - dtrans;
        if (pb_hi < pb_lo)
        {
            pb_lo = sp->fc1; pb_hi = sp->fc2; /* degenerate spec */
        }
        sb_lo = 0.0; sb_hi = 0.5 * sp->fc1;
    }
    else if (sp->fir_type == FCE_FIR_BANDSTOP)
    {
        pb_lo = 0.0; pb_hi = sp->fc1 - dtrans;
        if (pb_hi < 0.0)
            pb_hi = sp->fc1;
        /* fc1/fc2 are MIDPOINTS of the transition bands (each sits at
         * -6 dB), so scanning the full [fc1, fc2] would report ~6 dB
         * of "attenuation".  Measure the inner half of the stopband,
         * where the window's main-lobe floor actually is. (Same
         * rationale as the stopband insets used for LP/HP above.) */
        sb_lo = sp->fc1 + 0.25 * (sp->fc2 - sp->fc1);
        sb_hi = sp->fc2 - 0.25 * (sp->fc2 - sp->fc1);
    }
    else
    {
        pb_lo = 0.0; pb_hi = nyq;
        sb_lo = 0.0; sb_hi = nyq;
    }

    fce_scan_fir_band(h, n, fs, pb_lo, pb_hi, &mn, &mx);
    r->passband_ripple_measured_db = mx - mn;

    if (sp->fir_type == FCE_FIR_BANDSTOP)
    {
        /* the UPPER passband [fc2, nyq] is a real passband too - take
         * the worse ripple of the two (parity with the IIR bandstop). */
        double mn2, mx2, up_lo, up_hi;
        up_lo = sp->fc2 + dtrans;
        if (up_lo > nyq)
            up_lo = sp->fc2;
        up_hi = nyq;
        fce_scan_fir_band(h, n, fs, up_lo, up_hi, &mn2, &mx2);
        if ((mx2 - mn2) > (r->passband_ripple_measured_db))
            r->passband_ripple_measured_db = mx2 - mn2;
    }

    if (sp->fir_type == FCE_FIR_BANDPASS)
    {
        /* stopband below AND above the passband; take the worse case.
         * (mx from the band scan above is the PASSBAND max; the lower
         *  stopband needs its own scan.) */
        double mn1, mx1, mn2, mx2, mxw;
        fce_scan_fir_band(h, n, fs, sb_lo, sb_hi, &mn1, &mx1);
        fce_scan_fir_band(h, n, fs, sp->fc2 + 0.5 * (nyq - sp->fc2),
                          nyq, &mn2, &mx2);
        mxw = (mx1 > mx2) ? mx1 : mx2;
        r->stopband_atten_measured_db = -mxw;
    }
    else
    {
        fce_scan_fir_band(h, n, fs, sb_lo, sb_hi, &mn, &mx);
        r->stopband_atten_measured_db = -mx;
    }

    {
        double g;
        fce_cplx_t hc = fce_eval_fir_c(h, n, 0.0);
        g = fce_cx_abs(hc);
        r->dc_gain_db = 20.0 * log10(g > 0.0 ? g : 1e-300);
        hc = fce_eval_fir_c(h, n, FCE_PI);
        g = fce_cx_abs(hc);
        r->nyquist_gain_db = 20.0 * log10(g > 0.0 ? g : 1e-300);
    }
}

/* max |response_db(float) - response_db(quantized)| over a grid.
 * Gains below FLOOR are clamped so that deep-null noise does not dominate
 * the metric (a null shifted by a tiny coefficient change produces an
 * arbitrarily large dB difference). */
double fce_validate_quant_response(const double* ref, const double* test,
                                   uint32_t ns, int is_sos, double fs)
{
    const double floor = 1e-3; /* -60 dB: error is measured where the
                                 * response matters (passband/transition) */
    uint32_t i;
    uint32_t grid = FCE_VALIDATE_GRID_POINTS;
    double worst = 0.0;
    for (i = 0; i <= grid; i++)
    {
        double f = 0.5 * fs * (double)i / (double)grid;
        double w = 2.0 * FCE_PI * f / fs;
        fce_cplx_t a = is_sos ? fce_eval_sos_c(ref, ns, w)
                              : fce_eval_fir_c(ref, ns, w);
        fce_cplx_t b = is_sos ? fce_eval_sos_c(test, ns, w)
                              : fce_eval_fir_c(test, ns, w);
        double ga = fce_cx_abs(a);
        double gb = fce_cx_abs(b);
        /* skip points inside the deep-stopband noise floor: a null moved
         * by a tiny coefficient change gives a meaningless dB difference */
        if (ga < floor || gb < floor)
            continue;
        {
            double da = 20.0 * log10(ga);
            double db = 20.0 * log10(gb);
            double d = fabs(da - db);
            if (d > worst)
                worst = d;
        }
    }
    return worst;
}

#endif /* FCE_ENABLE_VALIDATION */
