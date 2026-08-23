/*
 * fce_iir.c - IIR coefficient generation.
 *
 * Pipeline (mirrors the SciPy reference implementation):
 *
 *   analog lowpass prototype (butter/cheby1/cheby2/ellip/bessel)
 *        -> frequency transformation (lp2lp/lp2hp/lp2bp/lp2bs)
 *        -> bilinear transform (with prewarping)
 *        -> zpk2sos ('nearest' pole-zero pairing)
 *        -> section ordering -> gain distribution
 *        -> float32/fixed-point output + validation
 *
 * Sign convention (see filtercoeff.h):
 *     y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]
 *     sos section layout: { b0, b1, b2, a1, a2 }   (a0 = 1 implicit)
 *
 * Internal math is float64 throughout; float32 is only the final rounding.
 */
#include "fce_internal.h"

/* ================================================================== */
/* numerical helpers                                                    */
/* ================================================================== */

/*
 * Re(prod(-num)/prod(-den)) without overflow.
 *
 * The naive formulation squares the accumulated magnitudes
 * (dot / (re^2+im^2)) and overflows beyond ~1e154 - e.g. the bilinear
 * gain compensation of a BP/BS design near order 22-25 silently
 * collapsed the filter gain to 0, or to NaN a few orders higher.
 * Here the ratio is accumulated factor-by-factor as a running complex
 * division, so no giant intermediate is ever squared.
 * Signs cancel per pair: (-a)/(-b) = a/b.
 */
static double fce_ratio_neg_products(const fce_cplx_t* num, uint32_t nnum,
                                     const fce_cplx_t* den, uint32_t nden)
{
    fce_cplx_t q = fce_cx(1.0, 0.0);
    uint32_t i, n = (nnum > nden) ? nnum : nden;
    for (i = 0; i < n; i++)
    {
        if (i < nnum)
            q = fce_cx_mul(q, fce_cx_scale(num[i], -1.0));
        if (i < nden)
            q = fce_cx_div(q, fce_cx_scale(den[i], -1.0));
    }
    return q.re;
}

/*
 * Re(prod(fs2 - z)/prod(fs2 - p)) without overflow.
 * Same idea as above; every factor is additionally scaled by 1/fs2
 * (corrected by fs2^(nz-np) at the end) so each complex factor stays
 * around unit magnitude.
 */
static double fce_ratio_bilinear(const fce_cplx_t* z, uint32_t nz,
                                 const fce_cplx_t* p, uint32_t np, double fs2)
{
    fce_cplx_t q = fce_cx(1.0, 0.0);
    double inv = 1.0 / fs2;
    uint32_t i, n = (np > nz) ? np : nz;
    for (i = 0; i < n; i++)
    {
        if (i < np)
            q = fce_cx_div(q, fce_cx_scale(fce_cx_sub(fce_cx(fs2, 0.0),
                                                      p[i]), inv));
        if (i < nz)
            q = fce_cx_mul(q, fce_cx_scale(fce_cx_sub(fce_cx(fs2, 0.0),
                                                      z[i]), inv));
    }
    return q.re * pow(fs2, (double)nz - (double)np);
}

/* ================================================================== */
/* analog lowpass prototypes (cutoff = 1 rad/s)                        */
/* ================================================================== */

static void proto_butter(uint32_t n, fce_cplx_t* p, double* k)
{
    uint32_t i;
    for (i = 0; i < n; i++)
    {
        int32_t m = (int32_t)(2u * i) - (int32_t)n + 1; /* -n+1 .. n-1 step 2 */
        double th = FCE_PI * (double)m / (2.0 * (double)n);
        p[i] = fce_cx(-cos(th), -sin(th)); /* -exp(j*th) */
    }
    *k = 1.0;
}

void proto_cheb1(uint32_t n, double rp, fce_cplx_t* p, double* k)
{
    double eps = sqrt(pow(10.0, 0.1 * rp) - 1.0);
    double mu = asinh(1.0 / eps) / (double)n;
    double prod = 1.0;
    uint32_t i;
    for (i = 0; i < n; i++)
    {
        int32_t m = (int32_t)(2u * i) - (int32_t)n + 1;
        double th = FCE_PI * (double)m / (2.0 * (double)n);
        double sh = sinh(mu);
        double ch = cosh(mu);
        p[i] = fce_cx(-sh * cos(th), -ch * sin(th));
        prod *= (sh * cos(th)) * (sh * cos(th)) + (ch * sin(th)) * (ch * sin(th));
    }
    *k = sqrt(prod);
    if ((n & 1u) == 0u)
        *k /= sqrt(1.0 + eps * eps);
}

void proto_cheb2(uint32_t n, double rs, fce_cplx_t* z, uint32_t* nz,
                        fce_cplx_t* p, double* k)
{
    double de = 1.0 / sqrt(pow(10.0, 0.1 * rs) - 1.0);
    double mu = asinh(1.0 / de) / (double)n;
    uint32_t i;

    *nz = 0;
    {
        uint32_t nz_total = (n & 1u) ? (n - 1u) : n; /* odd: zero at inf */
        for (i = 0; i < nz_total; i++)
        {
            int32_t m;
            if (n & 1u)
            {
                /* odd order: m = -(n-1), -(n-3), ..., -2, 2, ..., n-1 */
                uint32_t half = (n - 1u) / 2u;
                if (i < half)
                    m = -(int32_t)(2u * (half - i));
                else
                    m = (int32_t)(2u * (i - half + 1u));
            }
            else
            {
                m = (int32_t)(2u * i) - (int32_t)n + 1;
            }
            {
                double th = FCE_PI * (double)m / (2.0 * (double)n);
                z[(*nz)++] = fce_cx(0.0, 1.0 / sin(th));
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        int32_t m = (int32_t)(2u * i) - (int32_t)n + 1;
        double th = FCE_PI * (double)m / (2.0 * (double)n);
        double sh = sinh(mu);
        double ch = cosh(mu);
        fce_cplx_t s = fce_cx(sh * cos(th), ch * sin(th));
        p[i] = fce_cx_div(fce_cx(-1.0, 0.0), s);
    }

    /* k = Re(prod(-p) / prod(-z)) */
    *k = fce_ratio_neg_products(p, n, z, *nz);
}

fce_status_t proto_ellip(uint32_t n, double rp, double rs,
                                fce_cplx_t* z, uint32_t* nz,
                                fce_cplx_t* p, uint32_t* np, double* k)
{
    double eps_sq, eps, ck1_sq, val0, val1, m, capk, r, v0;
    uint32_t j;

    if (n == 1)
    {
        p[0] = fce_cx(-sqrt(1.0 / (pow(10.0, 0.1 * rp) - 1.0)), 0.0);
        *k = -p[0].re;
        *np = 1;
        *nz = 0;
        return FCE_OK;
    }

    eps_sq = pow(10.0, 0.1 * rp) - 1.0;
    eps = sqrt(eps_sq);
    ck1_sq = eps_sq / (pow(10.0, 0.1 * rs) - 1.0);
    if (!(ck1_sq > 0.0) || ck1_sq >= 1.0)
        return FCE_ERR_INVALID_SPEC;

    val0 = fce_ellipk(ck1_sq);
    val1 = fce_ellipkm1(ck1_sq);

    m = fce_ellipdeg(n, ck1_sq);
    if (!(m > 0.0) || !(m < 1.0))
        return FCE_ERR_NUMERICAL;
    capk = fce_ellipk(m);

    /* ---- zeros: z = 1j / (sqrt(m) * sn(u, m)),  u = jv*K/N ----
     * jv = 0,2,4,.. (odd n) or 1,3,5,.. (even n), n values total,
     * matching scipy's arange(1 - N%2, N, 2) */
    *nz = 0;
    for (j = 0; j < (n + 1u) / 2u; j++)
    {
        int32_t jv = (n & 1u) ? (int32_t)(2u * j) : (int32_t)(2u * j) + 1;
        double u = (double)jv * capk / (double)n;
        double s, c, d;
        fce_ellipj(u, m, &s, &c, &d);
        (void)c; (void)d;
        if (fabs(s) > 1e-14)
            z[(*nz)++] = fce_cx(0.0, 1.0 / (sqrt(m) * s));
    }
    {
        uint32_t cnt = *nz;
        for (j = 0; j < cnt; j++)
            z[(*nz)++] = fce_cx_conj(z[j]);
    }

    /* ---- r = arc_jac_sc1(1/eps, ck1_sq) = Im( asn(1j/eps, ck1_sq) ) ---- */
    {
        fce_cplx_t w = fce_arc_jac_sn(fce_cx(0.0, 1.0 / eps), ck1_sq);
        if (!fce_cx_isfinite(w) || fabs(w.re) > 1e-9)
            return FCE_ERR_NUMERICAL;
        r = w.im;
    }
    v0 = capk * r / ((double)n * val0);

    /* ---- poles ---- */
    {
        double sv, cv, dv;
        fce_ellipj(v0, 1.0 - m, &sv, &cv, &dv);
        *np = 0;
        for (j = 0; j < (n + 1u) / 2u; j++)
        {
            int32_t jv = (n & 1u) ? (int32_t)(2u * j) : (int32_t)(2u * j) + 1;
            double u = (double)jv * capk / (double)n;
            double s, c, d;
            double den;
            fce_cplx_t num;
            fce_ellipj(u, m, &s, &c, &d);
            /* p = -(c*d*sv*cv + j*s*dv) / (1 - (d*sv)^2) */
            num = fce_cx(-(c * d * sv * cv), -(s * dv));
            den = 1.0 - (d * sv) * (d * sv);
            if (fabs(den) < 1e-300)
                return FCE_ERR_NUMERICAL;
            p[(*np)++] = fce_cx(num.re / den, num.im / den);
        }
        if (n & 1u)
        {
            /* keep the real pole + complex ones, then add conjugates */
            double norm2 = 0.0;
            uint32_t keep = 0;
            fce_cplx_t tmp[FCE_MAX_IIR_ORDER];
            for (j = 0; j < *np; j++)
                norm2 += fce_cx_abs_sq(p[j]);
            for (j = 0; j < *np; j++)
                if (fabs(p[j].im) > 1e-14 * sqrt(norm2))
                    tmp[keep++] = p[j];
            for (j = 0; j < keep; j++)
                p[(*np)++] = fce_cx_conj(tmp[j]);
        }
        else
        {
            uint32_t base = *np;
            for (j = 0; j < base; j++)
                p[(*np)++] = fce_cx_conj(p[j]);
        }
    }

    /* ---- gain: k = Re(prod(-p)/prod(-z)); even order: /sqrt(1+eps^2) ---- */
    *k = fce_ratio_neg_products(p, *np, z, *nz);
    if ((n & 1u) == 0u)
        *k /= sqrt(1.0 + eps_sq);
    (void)val1;
    return FCE_OK;
}

fce_status_t proto_bessel(uint32_t n, fce_cplx_t* p, double* k)
{
    /* poles of the reverse Bessel polynomial (delay normalization) */
    fce_status_t st = fce_bessel_poles(n, p);
    double a_last = 1.0;
    double lo, hi;
    double target = 1.0 / sqrt(2.0);
    uint32_t i, it;
    double norm_factor;

    if (st != FCE_OK)
        return st;

    /* a_last = theta_n(0) = (2n)! / (2^n n!)  ->  DC gain = 1 */
    for (i = 0; i < n; i++)
        a_last *= (double)(2u * n - i);
    for (i = 0; i < n; i++)
        a_last /= 2.0;

    /* 'mag' normalization: find w so |H(jw)| = 1/sqrt(2),
     * H(s) = a_last / prod(s - p)   (delay-normalized) */
    lo = 0.1;
    hi = 1.5;
    for (it = 0; it < 200; it++)
    {
        fce_cplx_t prod = fce_cx(1.0, 0.0);
        for (i = 0; i < n; i++)
            prod = fce_cx_mul(prod, fce_cx_sub(fce_cx(0.0, hi), p[i]));
        if (a_last / fce_cx_abs(prod) < target)
            break;
        hi *= 2.0;
    }
    for (it = 0; it < 200; it++)
    {
        double mid = 0.5 * (lo + hi);
        fce_cplx_t prod_lo = fce_cx(1.0, 0.0);
        fce_cplx_t prod_mid = fce_cx(1.0, 0.0);
        double g_lo, g_mid;
        for (i = 0; i < n; i++)
            prod_lo = fce_cx_mul(prod_lo, fce_cx_sub(fce_cx(0.0, lo), p[i]));
        g_lo = a_last / fce_cx_abs(prod_lo);
        for (i = 0; i < n; i++)
            prod_mid = fce_cx_mul(prod_mid, fce_cx_sub(fce_cx(0.0, mid), p[i]));
        g_mid = a_last / fce_cx_abs(prod_mid);
        if ((g_lo - target) * (g_mid - target) <= 0.0)
            hi = mid;
        else
            lo = mid;
        if (hi - lo < 1e-14 * hi)
            break;
    }
    norm_factor = 0.5 * (lo + hi);

    for (i = 0; i < n; i++)
        p[i] = fce_cx_scale(p[i], 1.0 / norm_factor);
    *k = a_last * pow(norm_factor, -(double)n);
    return FCE_OK;
}

/* ================================================================== */
/* analog frequency transformations (zpk)                              */
/* ================================================================== */

static void tr_lp2lp(fce_cplx_t* z, uint32_t nz, fce_cplx_t* p, uint32_t np,
                     double* k, double wo)
{
    uint32_t i;
    for (i = 0; i < nz; i++)
        z[i] = fce_cx_scale(z[i], wo);
    for (i = 0; i < np; i++)
        p[i] = fce_cx_scale(p[i], wo);
    *k *= pow(wo, (double)(np - nz));
}

static void tr_lp2hp(fce_cplx_t* z, uint32_t* nz, fce_cplx_t* p, uint32_t np,
                     double* k, double wo)
{
    uint32_t i, deg = np - *nz;

    /* gain compensation uses the ORIGINAL prototype values
     * (scipy: k_hp = k * Re(prod(-z) / prod(-p))) */
    *k *= fce_ratio_neg_products(z, *nz, p, np);

    for (i = 0; i < *nz; i++)
        z[i] = fce_cx_div(fce_cx(wo, 0.0), z[i]);
    for (i = 0; i < np; i++)
        p[i] = fce_cx_div(fce_cx(wo, 0.0), p[i]);
    for (i = 0; i < deg; i++)
        z[(*nz)++] = fce_cx(0.0, 0.0);
}

static void tr_lp2bp(const fce_cplx_t* z, uint32_t nz,
                     const fce_cplx_t* p, uint32_t np,
                     double k, double wo, double bw,
                     fce_cplx_t* zo, uint32_t* nzo,
                     fce_cplx_t* po, uint32_t* npo, double* ko)
{
    uint32_t i, deg = np - nz;
    for (i = 0; i < nz; i++)
    {
        fce_cplx_t a = fce_cx_scale(z[i], 0.5 * bw);
        fce_cplx_t s = fce_cx_sqrt(fce_cx_sub(fce_cx_mul(a, a),
                                              fce_cx(wo * wo, 0.0)));
        zo[i] = fce_cx_add(a, s);
        zo[nz + i] = fce_cx_sub(a, s);
    }
    for (i = 0; i < np; i++)
    {
        fce_cplx_t a = fce_cx_scale(p[i], 0.5 * bw);
        fce_cplx_t s = fce_cx_sqrt(fce_cx_sub(fce_cx_mul(a, a),
                                              fce_cx(wo * wo, 0.0)));
        po[i] = fce_cx_add(a, s);
        po[np + i] = fce_cx_sub(a, s);
    }
    *nzo = 2u * nz + deg;
    *npo = 2u * np;
    for (i = 2u * nz; i < *nzo; i++)
        zo[i] = fce_cx(0.0, 0.0);
    *ko = k * pow(bw, (double)deg);
}

static void tr_lp2bs(const fce_cplx_t* z, uint32_t nz,
                     const fce_cplx_t* p, uint32_t np,
                     double k, double wo, double bw,
                     fce_cplx_t* zo, uint32_t* nzo,
                     fce_cplx_t* po, uint32_t* npo, double* ko)
{
    uint32_t i, deg = np - nz;
    for (i = 0; i < nz; i++)
    {
        fce_cplx_t a = fce_cx_div(fce_cx(0.5 * bw, 0.0), z[i]);
        fce_cplx_t s = fce_cx_sqrt(fce_cx_sub(fce_cx_mul(a, a),
                                              fce_cx(wo * wo, 0.0)));
        zo[i] = fce_cx_add(a, s);
        zo[nz + i] = fce_cx_sub(a, s);
    }
    for (i = 0; i < np; i++)
    {
        fce_cplx_t a = fce_cx_div(fce_cx(0.5 * bw, 0.0), p[i]);
        fce_cplx_t s = fce_cx_sqrt(fce_cx_sub(fce_cx_mul(a, a),
                                              fce_cx(wo * wo, 0.0)));
        po[i] = fce_cx_add(a, s);
        po[np + i] = fce_cx_sub(a, s);
    }
    *nzo = 2u * nz;
    *npo = 2u * np;

    /* gain compensation uses the ORIGINAL prototype values
     * (scipy: k_bs = k * Re(prod(-z) / prod(-p))) */
    *ko = k * fce_ratio_neg_products(z, nz, p, np);

    for (i = 0; i < deg; i++)
    {
        zo[(*nzo)++] = fce_cx(0.0, wo);
        zo[(*nzo)++] = fce_cx(0.0, -wo);
    }
}

/* ================================================================== */
/* bilinear transform (Tustin), s = 2 fs (z-1)/(z+1)                   */
/* ================================================================== */

static void tr_bilinear(fce_cplx_t* z, uint32_t* nz,
                        fce_cplx_t* p, uint32_t np, double* k, double fs)
{
    uint32_t i, deg = np - *nz;
    double fs2 = 2.0 * fs;

    /* gain compensation uses the ANALOG values, so compute it first
     * (scipy: k_z = k * Re(prod(fs2 - z) / prod(fs2 - p))) */
    *k *= fce_ratio_bilinear(z, *nz, p, np, fs2);

    /* map poles and zeros: z -> (2fs + z)/(2fs - z) */
    for (i = 0; i < *nz; i++)
    {
        fce_cplx_t num = fce_cx_add(fce_cx(fs2, 0.0), z[i]);
        fce_cplx_t den = fce_cx_sub(fce_cx(fs2, 0.0), z[i]);
        z[i] = fce_cx_div(num, den);
    }
    for (i = 0; i < np; i++)
    {
        fce_cplx_t num = fce_cx_add(fce_cx(fs2, 0.0), p[i]);
        fce_cplx_t den = fce_cx_sub(fce_cx(fs2, 0.0), p[i]);
        p[i] = fce_cx_div(num, den);
    }
    for (i = 0; i < deg; i++)
        z[(*nz)++] = fce_cx(-1.0, 0.0); /* zeros at infinity -> z = -1 */
}

/* ================================================================== */
/* zpk -> SOS ('nearest' pairing, mirroring scipy zpk2sos)             */
/* ================================================================== */

#define FCE_IS_REAL(c) (fabs((c).im) <= 1e-12 * (1.0 + fabs((c).re)))

static double fce_dist2(fce_cplx_t a, fce_cplx_t b)
{
    double dr = a.re - b.re;
    double di = a.im - b.im;
    return dr * dr + di * di;
}

/* mark `v` and its conjugate mate as used in the used[] array; returns
 * the number of elements consumed (1 for real v, 2 for complex v).
 * Note: for real v, conj(v) == v, so only ONE element is consumed. */
static uint32_t fce_consume(fce_cplx_t* arr, uint8_t* used, uint32_t n,
                            fce_cplx_t v)
{
    uint32_t i;
    uint32_t consumed = 0;
    int need_conj = !FCE_IS_REAL(v);
    for (i = 0; i < n && consumed < (need_conj ? 2u : 1u); i++)
    {
        if (used[i])
            continue;
        if (fce_dist2(arr[i], v) < 1e-20)
        {
            used[i] = 1;
            consumed++;
        }
        else if (need_conj && fce_dist2(arr[i], fce_cx_conj(v)) < 1e-20)
        {
            used[i] = 1;
            consumed++;
        }
    }
    return consumed;
}

static uint32_t fce_nearest_zero(const fce_cplx_t* z, const uint8_t* used,
                                 uint32_t nz, fce_cplx_t p1, int kind)
{
    /* kind: 0 any, 1 real only, 2 complex only */
    uint32_t best = nz;
    double best_d = 0.0;
    uint32_t i;
    for (i = 0; i < nz; i++)
    {
        double d;
        if (used[i])
            continue;
        if (kind == 1 && !FCE_IS_REAL(z[i]))
            continue;
        if (kind == 2 && FCE_IS_REAL(z[i]))
            continue;
        d = fce_dist2(z[i], p1);
        if (best == nz || d < best_d)
        {
            best = i;
            best_d = d;
        }
    }
    return best;
}

/* build one section from up to 2 poles and 2 zeros:
 * b = prod (1 - z_i z^-1), a = prod (1 - p_i z^-1) -> {b0,b1,b2,a1,a2} */
static void fce_section_make(const fce_cplx_t* zs, uint32_t nzs,
                             const fce_cplx_t* ps, uint32_t nps,
                             double* sos)
{
    double b[3] = {0.0, 0.0, 1.0};
    double a[3] = {0.0, 0.0, 1.0};

    if (nzs == 1)
    {
        b[0] = 0.0;
        b[1] = 1.0;
        b[2] = -zs[0].re;
    }
    else if (nzs == 2)
    {
        b[0] = 1.0;
        b[1] = -(zs[0].re + zs[1].re);
        b[2] = zs[0].re * zs[1].re - zs[0].im * zs[1].im;
    }
    if (nps == 1)
    {
        a[0] = 0.0;
        a[1] = 1.0;
        a[2] = -ps[0].re;
    }
    else if (nps == 2)
    {
        a[0] = 1.0;
        a[1] = -(ps[0].re + ps[1].re);
        a[2] = ps[0].re * ps[1].re - ps[0].im * ps[1].im;
    }

    sos[0] = b[0];
    sos[1] = b[1];
    sos[2] = b[2];
    sos[3] = a[1];
    sos[4] = a[2];
}

/* ================================================================== */
/* automatic order selection                                           */
/* ================================================================== */

static double fce_acosh(double x)
{
    return log(x + sqrt(x * x - 1.0));
}

/* ================================================================== */
/* bandstop order minimization (scipy band_stop_obj + fminbound)       */
/* ================================================================== */

/*
 * Modern scipy (>= ~1.10) does NOT compute the bandstop order from the
 * user's edges directly: the classic min(n1, n2) transition ratio
 * systematically over-estimates the order. Instead _find_nat_freq
 * optimizes each PASSBAND edge (fminbound of band_stop_obj) so that
 * the fractional order is minimal, and derives nat from the optimized
 * edges. Mirror that here so fce orders match scipy's.
 */

typedef struct fce_bs_obj_ctx
{
    double passb0, passb1;   /* prewarped passband edges (fixed sides) */
    double stopb0, stopb1;   /* prewarped stopband edges               */
    double gpass, gstop;     /* linear (10^(db/10))                    */
    int      family;         /* fce_iir_family_t                       */
    int      ind;            /* which passband edge is optimized (0/1) */
} fce_bs_obj_ctx_t;

static double fce_bs_obj(double wp, const fce_bs_obj_ctx_t* c)
{
    double p0 = (c->ind == 0) ? wp : c->passb0;
    double p1 = (c->ind == 1) ? wp : c->passb1;
    double n1 = fabs(c->stopb0 * (p0 - p1) /
                     (c->stopb0 * c->stopb0 - p0 * p1));
    double n2 = fabs(c->stopb1 * (p0 - p1) /
                     (c->stopb1 * c->stopb1 - p0 * p1));
    double nat = (n1 < n2) ? n1 : n2;

    switch (c->family)
    {
    case FCE_IIR_BUTTERWORTH:
        return log10((c->gstop - 1.0) / (c->gpass - 1.0)) /
               (2.0 * log10(nat));
    case FCE_IIR_ELLIPTIC:
    {
        double arg1_sq = (c->gpass - 1.0) / (c->gstop - 1.0);
        double arg0 = 1.0 / nat;
        double d00 = fce_ellipk(arg0 * arg0);
        double d01 = fce_ellipkm1(arg0 * arg0);
        double d10 = fce_ellipk(arg1_sq);
        double d11 = fce_ellipkm1(arg1_sq);
        return d00 * d11 / (d01 * d10);
    }
    default: /* cheby 1 & 2 */
        return fce_acosh(sqrt((c->gstop - 1.0) / (c->gpass - 1.0))) /
               fce_acosh(nat);
    }
}

/*
 * Bounded Brent minimization: a 1:1 port of scipy's
 * _minimize_scalar_bounded (golden section + parabolic interpolation,
 * xatol = 1e-5, maxfun = 500). Returns the minimizer of fn in [x1, x2].
 */
static double fce_fminbound(fce_bs_obj_ctx_t* ctx, double x1, double x2)
{
    const double sqrt_eps = 1.4901161193847656e-08; /* sqrt(2.2e-16)   */
    const double golden_mean = 0.38196601125010510; /* 0.5*(3-sqrt(5)) */
    const double xatol = 1e-5;
    double a = x1, b = x2;
    double fulc, nfc, xf;
    double ffulc, fnfc, fx, fu = 1e300;
    double rat = 0.0, e = 0.0, x;
    double xm, tol1, tol2;
    int num = 1;
    int golden;

    fulc = nfc = xf = a + golden_mean * (b - a);
    fx = fce_bs_obj(xf, ctx);
    ffulc = fnfc = fx;
    xm = 0.5 * (a + b);
    tol1 = sqrt_eps * fabs(xf) + xatol / 3.0;
    tol2 = 2.0 * tol1;

    while (fabs(xf - xm) > (tol2 - 0.5 * (b - a)))
    {
        double p = 0.0, q = 0.0, r;
        golden = 1;
        if (fabs(e) > tol1)
        {
            golden = 0;
            r = (xf - nfc) * (fx - ffulc);
            q = (xf - fulc) * (fx - fnfc);
            p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = fabs(q);
            r = e;
            e = rat;
            if ((fabs(p) < fabs(0.5 * q * r)) && (p > q * (a - xf)) &&
                (p < q * (b - xf)))
            {
                rat = p / q;
                x = xf + rat;
                if ((x - a) < tol2 || (b - x) < tol2)
                {
                    /* scipy: si = sign(xm - xf) + (xm == xf) -> +1 on tie */
                    double si = (xm >= xf) ? 1.0 : -1.0;
                    rat = tol1 * si;
                }
            }
            else
            {
                golden = 1;
            }
        }
        if (golden)
        {
            e = (xf >= xm) ? (a - xf) : (b - xf);
            rat = golden_mean * e;
        }
        {
            /* scipy: si = sign(rat) + (rat == 0) -> +1 on zero */
            double si = (rat > 0.0) ? 1.0 : (rat < 0.0 ? -1.0 : 1.0);
            double marat = (fabs(rat) > tol1) ? fabs(rat) : tol1;
            x = xf + si * marat;
        }
        fu = fce_bs_obj(x, ctx);
        if (++num >= 500)
            break;
        if (fu <= fx)
        {
            if (x >= xf)
                a = xf;
            else
                b = xf;
            fulc = nfc; ffulc = fnfc;
            nfc = xf;   fnfc = fx;
            xf = x;     fx = fu;
        }
        else
        {
            if (x < xf)
                a = x;
            else
                b = x;
            if ((fu <= fnfc) || (nfc == xf))
            {
                fulc = nfc; ffulc = fnfc;
                nfc = x;    fnfc = fu;
            }
            else if ((fu <= ffulc) || (fulc == xf) || (fulc == nfc))
            {
                fulc = x;   ffulc = fu;
            }
        }
        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * fabs(xf) + xatol / 3.0;
        tol2 = 2.0 * tol1;
    }
    return xf;
}

fce_status_t fce_iir_auto_order(const fce_spec_t* sp, fce_auto_t* a)
{
    double wp1, wp2, ws1, ws2;   /* band edges [Hz] */
    double passb1 = 0.0, passb2 = 0.0, stopb1 = 0.0, stopb2 = 0.0; /* prewarped */
    double nat, gpass, gstop, gpass_db, gstop_db;
    uint32_t n = 0;
    double fs = sp->fs;

    a->clamped = 0;
    gpass_db = (sp->passband_ripple_db > 0.0) ? sp->passband_ripple_db : 3.0;
    gstop_db = sp->stopband_atten_db;
    if (!(gstop_db > gpass_db))
        return FCE_ERR_INVALID_SPEC;

    switch (sp->iir_type)
    {
    case FCE_IIR_LOWPASS:
        wp1 = sp->fc1; ws1 = sp->edge1_hz;
        if (!(ws1 > wp1))
            return FCE_ERR_INVALID_SPEC;
        passb1 = fce_prewarp(wp1, fs);
        stopb1 = fce_prewarp(ws1, fs);
        nat = stopb1 / passb1;
        a->design_fc1 = 0.0; a->design_fc2 = 0.0;
        break;
    case FCE_IIR_HIGHPASS:
        wp1 = sp->fc1; ws1 = sp->edge1_hz;
        if (!(ws1 < wp1))
            return FCE_ERR_INVALID_SPEC;
        passb1 = fce_prewarp(wp1, fs);
        stopb1 = fce_prewarp(ws1, fs);
        nat = passb1 / stopb1;
        a->design_fc1 = 0.0; a->design_fc2 = 0.0;
        break;
    case FCE_IIR_BANDPASS:
        wp1 = sp->fc1; wp2 = sp->fc2;
        ws1 = sp->edge1_hz; ws2 = sp->edge2_hz;
        if (!(ws1 < wp1 && wp2 < ws2))
            return FCE_ERR_INVALID_SPEC;
        passb1 = fce_prewarp(wp1, fs);
        passb2 = fce_prewarp(wp2, fs);
        stopb1 = fce_prewarp(ws1, fs);
        stopb2 = fce_prewarp(ws2, fs);
        {
            double n1 = fabs((stopb1 * stopb1 - passb1 * passb2) /
                             (stopb1 * (passb1 - passb2)));
            double n2 = fabs((stopb2 * stopb2 - passb1 * passb2) /
                             (stopb2 * (passb1 - passb2)));
            nat = (n1 < n2) ? n1 : n2;
        }
        a->design_fc1 = 0.0; a->design_fc2 = 0.0;
        break;
    case FCE_IIR_BANDSTOP:
        ws1 = sp->fc1; ws2 = sp->fc2;      /* stopband edges */
        wp1 = sp->edge1_hz; wp2 = sp->edge2_hz; /* passband edges */
        if (!(wp1 < ws1 && ws2 < wp2))
            return FCE_ERR_INVALID_SPEC;
        passb1 = fce_prewarp(wp1, fs);
        passb2 = fce_prewarp(wp2, fs);
        stopb1 = fce_prewarp(ws1, fs);
        stopb2 = fce_prewarp(ws2, fs);
        {
            /* scipy _find_nat_freq: optimize each passband edge so the
             * fractional order is minimal (the plain min(n1, n2) over
             * estimates the bandstop order), then derive nat from the
             * optimized edges. */
            fce_bs_obj_ctx_t ctx;
            double opt0, opt1;
            ctx.stopb0 = stopb1;
            ctx.stopb1 = stopb2;
            ctx.gpass = pow(10.0, 0.1 * gpass_db);
            ctx.gstop = pow(10.0, 0.1 * gstop_db);
            ctx.family = (int)sp->iir_family;

            /* both minimizations run against the ORIGINAL edges;
             * each keeps the other passband edge fixed */
            ctx.passb0 = passb1;
            ctx.passb1 = passb2;
            ctx.ind = 0;
            opt0 = fce_fminbound(&ctx, passb1, stopb1 - 1e-12);
            ctx.ind = 1;
            opt1 = fce_fminbound(&ctx, stopb2 + 1e-12, passb2);
            passb1 = opt0;
            passb2 = opt1;

            {
                double n1 = fabs((stopb1 * (passb1 - passb2)) /
                                 (stopb1 * stopb1 - passb1 * passb2));
                double n2 = fabs((stopb2 * (passb1 - passb2)) /
                                 (stopb2 * stopb2 - passb1 * passb2));
                nat = (n1 < n2) ? n1 : n2;
            }
        }
        a->design_fc1 = 0.0; a->design_fc2 = 0.0;
        break;
    default:
        return FCE_ERR_INVALID_SPEC;
    }

    if (!(nat > 1.0))
        return FCE_ERR_INVALID_SPEC;

    gpass = pow(10.0, 0.1 * gpass_db);
    gstop = pow(10.0, 0.1 * gstop_db);

    switch (sp->iir_family)
    {
    case FCE_IIR_BUTTERWORTH:
        n = (uint32_t)ceil(log10((gstop - 1.0) / (gpass - 1.0)) /
                           (2.0 * log10(nat)));
        break;
    case FCE_IIR_CHEBYSHEV1:
    case FCE_IIR_CHEBYSHEV2:
        n = (uint32_t)ceil(fce_acosh(sqrt((gstop - 1.0) / (gpass - 1.0))) /
                           fce_acosh(nat));
        break;
    case FCE_IIR_ELLIPTIC:
    {
        double arg1_sq = (gpass - 1.0) / (gstop - 1.0);
        double arg0 = 1.0 / nat;
        double d00 = fce_ellipk(arg0 * arg0);
        double d01 = fce_ellipkm1(arg0 * arg0);
        double d10 = fce_ellipk(arg1_sq);
        double d11 = fce_ellipkm1(arg1_sq);
        n = (uint32_t)ceil(d00 * d11 / (d01 * d10));
        break;
    }
    case FCE_IIR_BESSEL:
        return FCE_ERR_UNSUPPORTED; /* Bessel has no closed-form order rule */
    default:
        return FCE_ERR_INVALID_SPEC;
    }

    if (n < 1u)
        return FCE_ERR_INVALID_SPEC;
    if (n > FCE_MAX_AUTO_ORDER)
    {
        n = FCE_MAX_AUTO_ORDER;
        a->clamped = 1; /* reported via FCE_FLAG_ORDER_CLAMPED by the caller */
    }

    a->order = n;

    /* ---- design frequencies so the passband spec is met exactly ---- */
    {
        double w0 = 1.0;
        if (sp->iir_family == FCE_IIR_BUTTERWORTH)
            w0 = pow(gpass - 1.0, -1.0 / (2.0 * (double)n));

        switch (sp->iir_type)
        {
        case FCE_IIR_LOWPASS:
            if (sp->iir_family == FCE_IIR_CHEBYSHEV2)
            {
                a->design_fc1 = sp->edge1_hz; /* stopband edge */
            }
            else
            {
                double wc = passb1 * w0;
                a->design_fc1 = fs * atan(wc / (2.0 * fs)) / FCE_PI;
            }
            break;
        case FCE_IIR_HIGHPASS:
            if (sp->iir_family == FCE_IIR_CHEBYSHEV2)
            {
                a->design_fc1 = sp->edge1_hz;
            }
            else
            {
                double wc = passb1 / w0;
                a->design_fc1 = fs * atan(wc / (2.0 * fs)) / FCE_PI;
            }
            break;
        case FCE_IIR_BANDPASS:
        {
            double bw = passb2 - passb1;
            double w0b = (sp->iir_family == FCE_IIR_CHEBYSHEV2) ? 1.0 : w0;
            double wn_hi = sqrt(w0b * w0b * bw * bw / 4.0 + passb1 * passb2)
                           + w0b * bw * 0.5;
            double wn_lo = sqrt(w0b * w0b * bw * bw / 4.0 + passb1 * passb2)
                           - w0b * bw * 0.5;
            if (sp->iir_family == FCE_IIR_CHEBYSHEV2)
            {
                a->design_fc1 = sp->edge1_hz;
                a->design_fc2 = sp->edge2_hz;
            }
            else
            {
                a->design_fc1 = fs * atan(wn_lo / (2.0 * fs)) / FCE_PI;
                a->design_fc2 = fs * atan(wn_hi / (2.0 * fs)) / FCE_PI;
            }
            break;
        }
        case FCE_IIR_BANDSTOP:
        if (sp->iir_family == FCE_IIR_CHEBYSHEV2)
        {
            /* cheb2ord (filter_type=3): the -gpass response frequency is
             * re-derived via new_freq so that BOTH specs are met with
             * the integer-rounded order; not just the stopband edges */
            double v = fce_acosh(sqrt((gstop - 1.0) / (gpass - 1.0)));
            double new_freq = 1.0 / cosh(v / (double)n);
            double nat0 = (new_freq * 0.5 * (passb1 - passb2)) +
                sqrt(new_freq * new_freq * (passb2 - passb1) *
                     (passb2 - passb1) * 0.25 + passb1 * passb2);
            double nat1 = passb1 * passb2 / nat0;
            double lo = (nat0 < nat1) ? nat0 : nat1;
            double hi = (nat0 < nat1) ? nat1 : nat0;
            a->design_fc1 = fs * atan(lo / (2.0 * fs)) / FCE_PI;
            a->design_fc2 = fs * atan(hi / (2.0 * fs)) / FCE_PI;
        }
        else
        {
            /* butter: discr formula with the optimized edges;
             * cheb1/ellip: reduces to the optimized edges themselves */
            double bw = passb2 - passb1;
            double discr = sqrt(bw * bw + 4.0 * w0 * w0 * passb1 * passb2);
            double wn_hi = (bw + discr) / (2.0 * w0);
            double wn_lo = (discr - bw) / (2.0 * w0);
            a->design_fc1 = fs * atan(wn_lo / (2.0 * fs)) / FCE_PI;
            a->design_fc2 = fs * atan(wn_hi / (2.0 * fs)) / FCE_PI;
        }
        break;
        default:
            break;
        }
    }
    return FCE_OK;
}

/* ================================================================== */
/* main IIR design                                                     */
/* ================================================================== */

fce_status_t fce_iir_design(const fce_spec_t* sp, fce_result_t* r,
                            fce_layout_t* lay, void* base)
{
    fce_cplx_t* pp = (fce_cplx_t*)(void*)((char*)base + lay->off_proto_p);
    fce_cplx_t* pz = (fce_cplx_t*)(void*)((char*)base + lay->off_proto_z);
    fce_cplx_t* ap = (fce_cplx_t*)(void*)((char*)base + lay->off_ap);
    fce_cplx_t* az = (fce_cplx_t*)(void*)((char*)base + lay->off_az);
    fce_cplx_t* dp = (fce_cplx_t*)(void*)((char*)base + lay->off_dp);
    fce_cplx_t* dz = (fce_cplx_t*)(void*)((char*)base + lay->off_dz);
    double* sos = (double*)(void*)((char*)base + lay->off_sos);
    double* sec_gains = (double*)(void*)((char*)base + lay->off_sec_gains);
    uint32_t n = lay->order;
    uint32_t nz = 0, np = 0;
    double k = 1.0;
    fce_status_t st;
    uint32_t i;
    double w1 = 0.0, w2 = 0.0;
    double fs = sp->fs;

    /* ---- analog prototype ---- */
    switch (sp->iir_family)
    {
    case FCE_IIR_BUTTERWORTH:
        proto_butter(n, pp, &k);
        np = n; nz = 0;
        break;
    case FCE_IIR_CHEBYSHEV1:
        proto_cheb1(n, sp->passband_ripple_db, pp, &k);
        np = n; nz = 0;
        break;
    case FCE_IIR_CHEBYSHEV2:
        proto_cheb2(n, sp->stopband_atten_db, pz, &nz, pp, &k);
        np = n;
        break;
    case FCE_IIR_ELLIPTIC:
        st = proto_ellip(n, sp->passband_ripple_db, sp->stopband_atten_db,
                         pz, &nz, pp, &np, &k);
        if (st != FCE_OK)
            return st;
        break;
    case FCE_IIR_BESSEL:
#if FCE_ENABLE_IIR_BESSEL
        st = proto_bessel(n, pp, &k);
        if (st != FCE_OK)
            return st;
        np = n; nz = 0;
        break;
#else
        return FCE_ERR_UNSUPPORTED;
#endif
    default:
        return FCE_ERR_UNSUPPORTED;
    }

    /* ---- analog frequency transformation (frequencies already
     * prewarped into r->design_fc1/fc2 by the caller) ---- */
    if (sp->iir_type == FCE_IIR_BANDPASS || sp->iir_type == FCE_IIR_BANDSTOP)
    {
        w1 = fce_prewarp(r->design_fc1, fs);
        w2 = fce_prewarp(r->design_fc2, fs);
    }

    switch (sp->iir_type)
    {
    case FCE_IIR_LOWPASS:
        tr_lp2lp(pz, nz, pp, np, &k, fce_prewarp(r->design_fc1, fs));
        break;
    case FCE_IIR_HIGHPASS:
        tr_lp2hp(pz, &nz, pp, np, &k, fce_prewarp(r->design_fc1, fs));
        break;
    case FCE_IIR_BANDPASS:
    {
        double wo = sqrt(w1 * w2);
        double bw = w2 - w1;
        tr_lp2bp(pz, nz, pp, np, k, wo, bw, az, &nz, ap, &np, &k);
        break;
    }
    case FCE_IIR_BANDSTOP:
    {
        double wo = sqrt(w1 * w2);
        double bw = w2 - w1;
        tr_lp2bs(pz, nz, pp, np, k, wo, bw, az, &nz, ap, &np, &k);
        break;
    }
    default:
        return FCE_ERR_UNSUPPORTED;
    }

    /* ---- bilinear ---- */
    if (sp->iir_type == FCE_IIR_BANDPASS || sp->iir_type == FCE_IIR_BANDSTOP)
    {
        for (i = 0; i < np; i++) dp[i] = ap[i];
        for (i = 0; i < nz; i++) dz[i] = az[i];
    }
    else
    {
        for (i = 0; i < np; i++) dp[i] = pp[i];
        for (i = 0; i < nz; i++) dz[i] = pz[i];
    }
    tr_bilinear(dz, &nz, dp, np, &k, fs);
    if (!isfinite(k))
        return FCE_ERR_NUMERICAL; /* hopeless conditioning */

    /* ---- expose digital poles/zeros (no black box) ---- */
    r->iir_poles = (const double*)(const void*)dp;
    r->iir_zeros = (const double*)(const void*)dz;
    r->iir_npoles = (uint16_t)np;
    r->iir_nzeros = (uint16_t)nz;

    /* ---- zpk -> sos (nearest pairing) ---- */
    {
        uint32_t ns;
        uint32_t npad = np, nzpad = nz;
        uint8_t up[2 * FCE_MAX_IIR_ORDER + 2];
        uint8_t uz[2 * FCE_MAX_IIR_ORDER + 2];
        uint32_t np_rem, nz_rem;
        int32_t si;

        while (nzpad < npad) dz[nzpad++] = fce_cx(0.0, 0.0);
        while (npad < nzpad) dp[npad++] = fce_cx(0.0, 0.0);
        if ((npad & 1u) != 0u)
        {
            dp[npad++] = fce_cx(0.0, 0.0);
            dz[nzpad++] = fce_cx(0.0, 0.0);
        }
        ns = (npad + 1u) / 2u;
        /* the layout reserves `lay->order` sections (BP/BS double the
         * prototype order, so ns can reach order, not (order+1)/2) */
        if (ns > lay->order)
            return FCE_ERR_BUFFER_TOO_SMALL;

        memset(up, 0, sizeof(up));
        memset(uz, 0, sizeof(uz));
        np_rem = npad;
        nz_rem = nzpad;

        for (si = (int32_t)ns - 1; si >= 0; si--)
        {
            uint32_t p1i = npad;
            double worst = 1e300;
            fce_cplx_t p1;
            uint32_t nreal_p = 0, nreal_z = 0;
            uint32_t j;

            /* worst remaining pole: min |1 - |p|| */
            for (j = 0; j < npad; j++)
            {
                double d;
                if (up[j])
                    continue;
                d = fabs(1.0 - fce_cx_abs(dp[j]));
                if (p1i == npad || d < worst)
                {
                    p1i = j;
                    worst = d;
                }
            }
            if (p1i == npad)
                return FCE_ERR_NUMERICAL;
            p1 = dp[p1i];
            np_rem -= fce_consume(dp, up, npad, p1);

            for (j = 0; j < npad; j++)
                if (!up[j] && FCE_IS_REAL(dp[j]))
                    nreal_p++;
            for (j = 0; j < nzpad; j++)
                if (!uz[j] && FCE_IS_REAL(dz[j]))
                    nreal_z++;

            if (FCE_IS_REAL(p1) && nreal_p == 0)
            {
                /* last real pole: pair with nearest real zero (+ origin) */
                uint32_t z1i = fce_nearest_zero(dz, uz, nzpad, p1, 1);
                fce_cplx_t zs[2], ps[2];
                if (z1i < nzpad)
                    nz_rem -= fce_consume(dz, uz, nzpad, dz[z1i]);
                zs[0] = (z1i < nzpad) ? dz[z1i] : fce_cx(0.0, 0.0);
                zs[1] = fce_cx(0.0, 0.0);
                ps[0] = p1;
                ps[1] = fce_cx(0.0, 0.0);
                fce_section_make(zs, 2, ps, 2, sos + 5u * (uint32_t)si);
            }
            else if (np_rem + 1u == nz_rem &&
                     !FCE_IS_REAL(p1) && nreal_p == 1 && nreal_z == 1)
            {
                /* must pair with a complex zero to keep the real zero */
                uint32_t z1i = fce_nearest_zero(dz, uz, nzpad, p1, 2);
                fce_cplx_t zs[2], ps[2];
                if (z1i == nzpad)
                    return FCE_ERR_NUMERICAL;
                nz_rem -= fce_consume(dz, uz, nzpad, dz[z1i]);
                zs[0] = dz[z1i];
                zs[1] = fce_cx_conj(dz[z1i]);
                ps[0] = p1;
                ps[1] = fce_cx_conj(p1);
                fce_section_make(zs, 2, ps, 2, sos + 5u * (uint32_t)si);
            }
            else
            {
                fce_cplx_t p2;
                fce_cplx_t zs[2] = {{0.0, 0.0}, {0.0, 0.0}};
                fce_cplx_t ps[2] = {{0.0, 0.0}, {0.0, 0.0}};
                uint32_t nzsec = 0, npsec = 0;

                if (FCE_IS_REAL(p1))
                {
                    uint32_t p2i = npad;
                    double worst2 = 1e300;
                    for (j = 0; j < npad; j++)
                    {
                        double d;
                        if (up[j] || !FCE_IS_REAL(dp[j]))
                            continue;
                        d = fabs(1.0 - fce_cx_abs(dp[j]));
                        if (p2i == npad || d < worst2)
                        {
                            p2i = j;
                            worst2 = d;
                        }
                    }
                    if (p2i == npad)
                        return FCE_ERR_NUMERICAL;
                    p2 = dp[p2i];
                    np_rem -= fce_consume(dp, up, npad, p2);
                }
                else
                {
                    p2 = fce_cx_conj(p1);
                }

                {
                    uint32_t z1i = fce_nearest_zero(dz, uz, nzpad, p1, 0);
                    if (z1i == nzpad)
                    {
                        ps[npsec++] = p1;
                        ps[npsec++] = p2;
                        fce_section_make(zs, 0, ps, npsec,
                                         sos + 5u * (uint32_t)si);
                        continue;
                    }
                    {
                        fce_cplx_t z1 = dz[z1i];
                        nz_rem -= fce_consume(dz, uz, nzpad, z1);
                        if (!FCE_IS_REAL(z1))
                        {
                            zs[nzsec++] = z1;
                            zs[nzsec++] = fce_cx_conj(z1);
                        }
                        else
                        {
                            uint32_t z2i = fce_nearest_zero(dz, uz, nzpad, p1, 1);
                            zs[nzsec++] = z1;
                            if (z2i < nzpad)
                            {
                                zs[nzsec++] = dz[z2i];
                                nz_rem -= fce_consume(dz, uz, nzpad, dz[z2i]);
                            }
                        }
                        ps[npsec++] = p1;
                        ps[npsec++] = p2;
                        fce_section_make(zs, nzsec, ps, npsec,
                                         sos + 5u * (uint32_t)si);
                    }
                }
            }
        }
        lay->n_sections = ns;
        r->num_sections = (uint16_t)ns;
    }

    /* ---- section ordering + gain distribution ---- */
    {
        /* sections can reach FCE_MAX_IIR_ORDER for BP/BS designs */
        double peaks[FCE_MAX_IIR_ORDER];
        double key[FCE_MAX_IIR_ORDER];
        uint32_t ns = lay->n_sections;
        if (ns > FCE_MAX_IIR_ORDER)
            return FCE_ERR_BUFFER_TOO_SMALL;

        for (i = 0; i < ns; i++)
        {
            double c[5];
            uint32_t j;
            for (j = 0; j < 5; j++)
                c[j] = sos[5u * i + j];
            peaks[i] = fce_biquad_peak_gain(c, 256u);
        }

        if (sp->sos_order == FCE_SOS_ORDER_POLE_RADIUS_DESC ||
            sp->sos_order == FCE_SOS_ORDER_INTERNAL_GAIN)
        {
            for (i = 0; i < ns; i++)
            {
                if (sp->sos_order == FCE_SOS_ORDER_INTERNAL_GAIN)
                {
                    key[i] = peaks[i];
                }
                else
                {
                    fce_cplx_t p1, p2;
                    double r1, r2;
                    fce_biquad_poles(sos[5u * i + 3], sos[5u * i + 4], &p1, &p2);
                    r1 = fce_cx_abs(p1);
                    r2 = fce_cx_abs(p2);
                    key[i] = (r1 > r2) ? r1 : r2;
                }
            }
            for (i = 0; i < ns; i++)
            {
                uint32_t best = i, j;
                for (j = i + 1; j < ns; j++)
                {
                    int take = (sp->sos_order == FCE_SOS_ORDER_POLE_RADIUS_DESC)
                                   ? (key[j] > key[best])
                                   : (key[j] < key[best]);
                    if (take)
                        best = j;
                }
                if (best != i)
                {
                    double tmp[5];
                    double ktmp;
                    uint32_t jj;
                    for (jj = 0; jj < 5; jj++)
                    {
                        tmp[jj] = sos[5u * i + jj];
                        sos[5u * i + jj] = sos[5u * best + jj];
                        sos[5u * best + jj] = tmp[jj];
                    }
                    ktmp = key[i]; key[i] = key[best]; key[best] = ktmp;
                    ktmp = peaks[i]; peaks[i] = peaks[best]; peaks[best] = ktmp;
                }
            }
        }

        /* distribute total gain into first section (scipy convention) */
        for (i = 0; i < ns; i++)
        {
            uint32_t j;
            double g = (i == 0) ? k : 1.0;
            for (j = 0; j < 3; j++)
                sos[5u * i + j] *= g;
            sec_gains[i] = peaks[i] * g;
        }
        r->sos_order = sp->sos_order;
    }

    /* numerical safety net: never return FCE_OK with NaN/Inf
     * coefficients (only possible for hopelessly ill-conditioned specs);
     * report it honestly as FCE_ERR_NUMERICAL instead */
    for (i = 0; i < 5u * lay->n_sections; i++)
        if (!isfinite(sos[i]))
            return FCE_ERR_NUMERICAL;

    /* ---- float32 output (float64 always available) ---- */
    r->sos_f64 = sos;
    {
        float* sos32 = (float*)(void*)((char*)base + lay->off_sos32);
        uint32_t ns = lay->n_sections;
        for (i = 0; i < 5u * ns; i++)
            sos32[i] = (float)sos[i];
        r->sos_f32 = (sp->precision == FCE_PRECISION_FLOAT32) ? sos32 : NULL;
    }

    r->order = (uint16_t)n;
    r->section_gains = sec_gains;
    r->norm_factor = 1.0; /* total gain is embedded in the SOS */
    return FCE_OK;
}
