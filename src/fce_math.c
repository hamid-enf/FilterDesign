/*
 * fce_math.c - numerical core: elliptic integrals/functions, I0, sinc,
 * Kaiser parameter estimation, Durand-Kerner polynomial root finder.
 *
 * The elliptic-filter machinery follows Orfanidis' "Lecture Notes on
 * Elliptic Filter Design" and the SciPy reference implementation.
 */
#include "fce_internal.h"

/* ================================================================== */
/* basic helpers                                                       */
/* ================================================================== */

double fce_sinc(double x)
{
    if (x == 0.0)
        return 1.0;
    {
        double px = FCE_PI * x;
        return sin(px) / px;
    }
}

double fce_i0(double x)
{
    /* I0(x) = sum_k (x^2/4)^k / (k!)^2 ; converges quickly for x <= ~50 */
    double sum = 1.0;
    double term = 1.0;
    double x2 = 0.25 * x * x;
    uint32_t k = 1;
    while (k < 200u) {
        term *= x2 / ((double)k * (double)k);
        sum += term;
        if (term < 1e-18 * sum)
            break;
        k++;
    }
    return sum;
}

/* ================================================================== */
/* complete elliptic integral of the first kind: K(m), 0 <= m < 1      */
/* ================================================================== */

double fce_ellipk(double m)
{
    double a, b, t;
    if (m <= 0.0)
        return FCE_PI / 2.0;
    if (m >= 1.0)
        return 1.0e300; /* diverges; caller must avoid */

    a = 1.0;
    b = sqrt(1.0 - m);
    /* arithmetic-geometric mean */
    for (;;)
    {
        t = 0.5 * (a + b);
        b = sqrt(a * b);
        a = t;
        if (a - b < 1e-15 * a)
            break;
    }
    return FCE_PI / (2.0 * a);
}

double fce_ellipkm1(double m)
{
    /* K(1-m) = pi / (2 * AGM(1, sqrt(m)))
     * Computing it as ellipk(1-m) would lose ~4.7e-7 relative precision
     * for tiny m due to cancellation in (1 - m). */
    double a, b, t;
    if (m <= 0.0)
        return 1.0e300;
    if (m >= 1.0)
        return FCE_PI / 2.0;
    a = 1.0;
    b = sqrt(m);
    for (;;)
    {
        t = 0.5 * (a + b);
        b = sqrt(a * b);
        a = t;
        if (a - b < 1e-15 * a)
            break;
    }
    return FCE_PI / (2.0 * a);
}

/* ================================================================== */
/* Jacobi elliptic functions sn/cn/dn (real u, 0 <= m <= 1)            */
/* Descending Landen / AGM (Abramowitz & Stegun 17.6)                  */
/* ================================================================== */

void fce_ellipj(double u, double m, double* sn, double* cn, double* dn)
{
    double a[64];
    double c[64];
    double b, phi, t;
    int n = 0;
    int i;

    if (sn) *sn = 0.0;
    if (cn) *cn = 1.0;
    if (dn) *dn = 1.0;

    if (m <= 0.0)
    {
        if (sn) *sn = sin(u);
        if (cn) *cn = cos(u);
        if (dn) *dn = 1.0;
        return;
    }
    if (m >= 1.0)
    {
        double th = tanh(u);
        if (sn) *sn = th;
        if (cn) *cn = 1.0 / cosh(u);
        if (dn) *dn = 1.0 / cosh(u);
        return;
    }

    a[0] = 1.0;
    b = sqrt(1.0 - m);
    c[0] = sqrt(m);

    for (n = 1; n < 64; n++)
    {
        a[n] = 0.5 * (a[n - 1] + b);
        c[n] = 0.5 * (a[n - 1] - b);
        b = sqrt(a[n - 1] * b);
        if (c[n] < 1e-16 * a[n])
            break;
    }

    phi = ldexp(a[n], n) * u;

    for (i = n; i > 0; i--)
    {
        t = (c[i] / a[i]) * sin(phi);
        if (t > 1.0) t = 1.0;
        if (t < -1.0) t = -1.0;
        phi = 0.5 * (asin(t) + phi);
    }

    if (sn) *sn = sin(phi);
    if (cn) *cn = cos(phi);
    if (dn) *dn = sqrt(1.0 - m * sin(phi) * sin(phi));
}

/* ================================================================== */
/* complex elementary functions                                        */
/* ================================================================== */

fce_cplx_t fce_cx_sqrt(fce_cplx_t a)
{
    double r = hypot(a.re, a.im);
    double re = sqrt(0.5 * (r + a.re));
    double im = (a.im >= 0.0 ? 1.0 : -1.0) * sqrt(0.5 * (r - a.re));
    return fce_cx(re, im);
}

fce_cplx_t fce_cx_ln(fce_cplx_t a)
{
    return fce_cx(log(hypot(a.re, a.im)), atan2(a.im, a.re));
}

fce_cplx_t fce_cx_asin(fce_cplx_t a)
{
    /* asin(z) = -j ln(j z + sqrt(1 - z^2)) */
    fce_cplx_t one = fce_cx(1.0, 0.0);
    fce_cplx_t jz = fce_cx(-a.im, a.re);
    fce_cplx_t s = fce_cx_sqrt(fce_cx_sub(one, fce_cx_mul(a, a)));
    fce_cplx_t l = fce_cx_ln(fce_cx_add(jz, s));
    return fce_cx(l.im, -l.re);
}

fce_cplx_t fce_cx_atanh(fce_cplx_t a)
{
    /* atanh(z) = 0.5 (ln(1+z) - ln(1-z)) */
    fce_cplx_t one = fce_cx(1.0, 0.0);
    fce_cplx_t l1 = fce_cx_ln(fce_cx_add(one, a));
    fce_cplx_t l2 = fce_cx_ln(fce_cx_sub(one, a));
    return fce_cx(0.5 * (l1.re - l2.re), 0.5 * (l1.im - l2.im));
}

int fce_cx_isfinite(fce_cplx_t a)
{
    return isfinite(a.re) && isfinite(a.im);
}

/* ================================================================== */
/* inverse Jacobi sn (complex w), ascending Landen (Orfanidis eq. 56)  */
/* ================================================================== */

static double fce_cx_complement(double kx)
{
    /* sqrt(1 - kx^2), stable for small kx */
    return sqrt((1.0 - kx) * (1.0 + kx));
}

fce_cplx_t fce_arc_jac_sn(fce_cplx_t w, double m)
{
    double k;
    double ks[64];
    double kp;
    double K = 1.0;
    fce_cplx_t wns[64];
    fce_cplx_t u;
    int niter = 0;
    int i;

    k = sqrt(m);
    if (k > 1.0)
        return fce_cx(NAN, NAN);
    if (k == 1.0)
        return fce_cx_atanh(w); /* sn(z,1)=tanh(z) -> z=atanh(w) */

    /* ascending Landen: build the ks sequence */
    ks[0] = k;
    while (ks[niter] != 0.0)
    {
        kp = fce_cx_complement(ks[niter]);
        ks[niter + 1] = (1.0 - kp) / (1.0 + kp);
        niter++;
        if (niter >= 62)
            return fce_cx(NAN, NAN);
    }

    /* K(k) via the Landen product */
    for (i = 1; i <= niter; i++)
        K *= (1.0 + ks[i]);
    K *= FCE_PI / 2.0;

    /* ascending recurrence for w */
    wns[0] = w;
    for (i = 0; i < niter; i++)
    {
        fce_cplx_t num = fce_cx_scale(wns[i], 2.0);
        fce_cplx_t arg = fce_cx_scale(wns[i], ks[i]);
        fce_cplx_t cplx_c = fce_cx_sqrt(fce_cx_sub(fce_cx(1.0, 0.0),
                                                   fce_cx_mul(arg, arg)));
        fce_cplx_t den = fce_cx_mul(fce_cx(1.0 + ks[i + 1], 0.0),
                                    fce_cx(1.0 + cplx_c.re, cplx_c.im));
        wns[i + 1] = fce_cx_div(num, den);
        if (!fce_cx_isfinite(wns[i + 1]))
            return fce_cx(NAN, NAN);
    }

    /* u = 2/pi * asin(wns[-1])   (complex asin, as in SciPy) */
    u = fce_cx_scale(fce_cx_asin(wns[niter]), 2.0 / FCE_PI);
    return fce_cx_scale(u, K);
}

/* ================================================================== */
/* elliptic degree equation (Orfanidis eq. 49)                          */
/* ================================================================== */

double fce_ellipdeg(uint32_t n, double m1)
{
    double K1 = fce_ellipk(m1);
    double K1p = fce_ellipkm1(m1);
    double q1 = exp(-FCE_PI * K1p / K1);
    double q;
    double num = 0.0;
    double den = 0.0;
    uint32_t k;

    if (n == 0)
        return NAN;
    q = pow(q1, 1.0 / (double)n);

    for (k = 0; k <= 8; k++)
        num += pow(q, (double)(k * (k + 1)));
    den = 1.0;
    for (k = 1; k <= 9; k++)
        den += 2.0 * pow(q, (double)(k * k));

    return 16.0 * q * pow(num / den, 4.0);
}

/* ================================================================== */
/* Aberth-Ehrlich simultaneous root finding                             */
/* ================================================================== */

static double fce_polyval(double x, const double* c, uint32_t n)
{
    /* c[0] + c[1] x + ... + c[n-1] x^(n-1) */
    double s = 0.0;
    uint32_t i;
    for (i = n; i-- > 0;)
        s = s * x + c[i];
    return s;
}

fce_status_t fce_roots_durand_kerner(const double* coeff, uint32_t deg,
                                     fce_cplx_t* roots)
{
    return fce_roots_durand_kerner_start(coeff, deg, roots, NULL);
}

fce_status_t fce_roots_durand_kerner_start(const double* coeff, uint32_t deg,
                                           fce_cplx_t* roots,
                                           const fce_cplx_t* start)
{
    double a[FCE_MAX_IIR_ORDER + 1];
    fce_cplx_t r[FCE_MAX_IIR_ORDER];
    uint32_t i, iter, j;
    int converged;

    if (deg == 0 || deg > FCE_MAX_IIR_ORDER)
        return FCE_ERR_INVALID_ARGUMENT;

    /* normalize by leading coefficient */
    {
        double lead = coeff[deg];
        if (!(lead > 0.0) && !(lead < 0.0))
            return FCE_ERR_NUMERICAL;
        for (i = 0; i <= deg; i++)
            a[i] = coeff[i] / lead;
    }

    if (start != NULL)
    {
        for (i = 0; i < deg; i++)
            r[i] = start[i];
    }
    else
    {
        /* default starts: circle around the origin, radius scaled by the
         * coefficient magnitudes */
        double scale = 0.0;
        for (i = 0; i <= deg; i++)
            scale += fabs(a[i]);
        scale /= (double)(deg + 1);
        if (!(scale > 0.0))
            return FCE_ERR_NUMERICAL;
        for (i = 0; i < deg; i++)
        {
            double ang = 2.0 * FCE_PI * (double)i / (double)deg;
            double rad = 0.7 + 0.3 * (double)((i * 7u) % 5u) / 4.0;
            r[i] = fce_cx(rad * cos(ang) * scale, rad * sin(ang) * scale);
        }
    }

    for (iter = 0; iter < 200u; iter++)
    {
        converged = 1;
        for (i = 0; i < deg; i++)
        {
            fce_cplx_t pv, dp;
            fce_cplx_t ssum = fce_cx(0.0, 0.0);
            fce_cplx_t w, den, step;

            /* p(r_i) and p'(r_i) via simultaneous Horner */
            pv = fce_cx(1.0, 0.0);
            dp = fce_cx(0.0, 0.0);
            for (j = deg; j-- > 0u;)
            {
                dp = fce_cx_add(fce_cx_mul(dp, r[i]), pv);
                pv = fce_cx_add(fce_cx_mul(pv, r[i]), fce_cx(a[j], 0.0));
            }
            for (j = 0; j < deg; j++)
            {
                if (j != i)
                    ssum = fce_cx_add(ssum,
                                      fce_cx_div(fce_cx(1.0, 0.0),
                                                 fce_cx_sub(r[i], r[j])));
            }
            /* Aberth-Ehrlich update: r -= w / (1 - w * sum 1/(r-rj)) */
            w = fce_cx_div(pv, dp);
            den = fce_cx_sub(fce_cx(1.0, 0.0), fce_cx_mul(w, ssum));
            if (fce_cx_abs(den) < 1e-300)
                return FCE_ERR_NUMERICAL;
            step = fce_cx_div(w, den);
            r[i] = fce_cx_sub(r[i], step);
            if (fce_cx_abs(step) > 1e-12 * (1.0 + fce_cx_abs(r[i])))
                converged = 0;
        }
        if (converged)
            break;
    }

    if (!converged)
        return FCE_ERR_NUMERICAL;

    for (i = 0; i < deg; i++)
        roots[i] = r[i];
    return FCE_OK;
}

/* ================================================================== */
/* Bessel prototype poles                                              */
/* ================================================================== */

/*
 * Approximate zeros of the ordinary Bessel polynomial y_n (Campos &
 * Calderon 2011, arXiv:1105.0957) - used as starting points, as in SciPy.
 */
void fce_campos_zeros(uint32_t n, fce_cplx_t* z)
{
    double s, r, b0, b1, b2, b3, a1, a2;
    uint32_t k;

    if (n == 1)
    {
        z[0] = fce_cx(-1.0, 0.0);
        return;
    }

    s = fce_polyval((double)n, (double[]){0.0, 0.0, 2.0, 0.0, -3.0, 1.0}, 6);
    b3 = fce_polyval((double)n, (double[]){16.0, -8.0}, 2) / s;
    b2 = fce_polyval((double)n, (double[]){-24.0, -12.0, 12.0}, 3) / s;
    b1 = fce_polyval((double)n, (double[]){8.0, 24.0, -12.0, -2.0}, 4) / s;
    b0 = fce_polyval((double)n, (double[]){0.0, -6.0, 0.0, 5.0, -1.0}, 5) / s;

    r = fce_polyval((double)n, (double[]){0.0, 0.0, 2.0, 1.0}, 4);
    a1 = fce_polyval((double)n, (double[]){-6.0, -6.0}, 2) / r;
    a2 = 6.0 / r;

    for (k = 1; k <= n; k++)
    {
        double kk = (double)k;
        z[k - 1] = fce_cx(a2 * kk * kk + a1 * kk,
                          ((b3 * kk + b2) * kk + b1) * kk + b0);
    }
}

/*
 * Poles of the Bessel filter (roots of the reverse Bessel polynomial
 * theta_n), 'delay' normalization (DC gain 1, group delay ~1 at DC).
 */
fce_status_t fce_bessel_poles(uint32_t n, fce_cplx_t* poles)
{
    double a[FCE_MAX_IIR_ORDER + 1];
    fce_cplx_t start[FCE_MAX_IIR_ORDER];
    fce_cplx_t yz[FCE_MAX_IIR_ORDER];
    uint32_t k;
    fce_status_t st;
    double sum = 0.0;

    if (n == 0 || n > FCE_MAX_IIR_ORDER)
        return FCE_ERR_INVALID_ARGUMENT;

    /* theta_n(s) = sum_k a_k s^k,  a_k = (2n-k)! / (2^(n-k) k! (n-k)!) */
    for (k = 0; k <= n; k++)
    {
        double num = 1.0;
        uint32_t i;
        double den = 1.0;
        for (i = 0; i < (2u * n - k); i++)
            num *= (double)(2u * n - k - i);
        den = ldexp(1.0, (int)(n - k));
        for (i = 2; i <= k; i++)
            den *= (double)i;
        for (i = 2; i <= (n - k); i++)
            den *= (double)i;
        a[k] = num / den;
    }

    fce_campos_zeros(n, yz);
    for (k = 0; k < n; k++)
        start[k] = fce_cx_div(fce_cx(1.0, 0.0), yz[k]); /* poles ~ 1/zeros */

    st = fce_roots_durand_kerner_start(a, n, poles, start);
    if (st != FCE_OK)
        return st;

    /* sanity: sum of roots == -a_{n-1}/a_n */
    for (k = 0; k < n; k++)
        sum += poles[k].re;
    if (fabs(sum + a[n - 1] / a[n]) > 1e-8 * (1.0 + fabs(sum)))
        return FCE_ERR_NUMERICAL;

    return FCE_OK;
}

/* ================================================================== */
/* biquad helpers                                                      */
/* ================================================================== */

/* H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2) at z^-1=e^-jw */
fce_cplx_t fce_eval_biquad(const double* c, double w)
{
    double cw = cos(w), sw = sin(w);
    double c2w = cos(2.0 * w), s2w = sin(2.0 * w);
    fce_cplx_t num = fce_cx(c[0] + c[1] * cw + c[2] * c2w,
                            -(c[1] * sw + c[2] * s2w));
    fce_cplx_t den = fce_cx(1.0 + c[3] * cw + c[4] * c2w,
                            -(c[3] * sw + c[4] * s2w));
    return fce_cx_div(num, den);
}

/* poles of z^2 + a1 z + a2 (roots of the denominator) */
void fce_biquad_poles(double a1, double a2, fce_cplx_t* p1, fce_cplx_t* p2)
{
    double disc = a1 * a1 - 4.0 * a2;
    if (disc >= 0.0)
    {
        double s = sqrt(disc);
        *p1 = fce_cx(0.5 * (-a1 + s), 0.0);
        *p2 = fce_cx(0.5 * (-a1 - s), 0.0);
    }
    else
    {
        double s = sqrt(-disc);
        *p1 = fce_cx(-0.5 * a1, 0.5 * s);
        *p2 = fce_cx(-0.5 * a1, -0.5 * s);
    }
}

/* peak |H| of a biquad over [0, pi] on a `grid`-point scan */
double fce_biquad_peak_gain(const double* c, uint32_t grid)
{
    double best = 0.0;
    uint32_t i;
    for (i = 0; i <= grid; i++)
    {
        double w = FCE_PI * (double)i / (double)grid;
        double g = fce_cx_abs(fce_eval_biquad(c, w));
        if (g > best)
            best = g;
    }
    return best;
}

double fce_biquad_group_delay(const double* c, double w)
{
    /* GD = -d/dw arg(H); per biquad via (B' B* - ...) analytic formula */
    double cw = cos(w), sw = sin(w);
    double c2w = cos(2.0 * w), s2w = sin(2.0 * w);
    double br = c[0] + c[1] * cw + c[2] * c2w;
    double bi = -(c[1] * sw + c[2] * s2w);
    double bdr = -(c[1] * sw + 2.0 * c[2] * s2w);
    double bdi = -(c[1] * cw + 2.0 * c[2] * c2w);
    double ar = 1.0 + c[3] * cw + c[4] * c2w;
    double ai = -(c[3] * sw + c[4] * s2w);
    double adr = -(c[3] * sw + 2.0 * c[4] * s2w);
    double adi = -(c[3] * cw + 2.0 * c[4] * c2w);
    double bm2 = br * br + bi * bi;
    double am2 = ar * ar + ai * ai;
    double gd_b = (bm2 > 0.0) ? (br * bdi - bi * bdr) / bm2 : 0.0;
    double gd_a = (am2 > 0.0) ? (ar * adi - ai * adr) / am2 : 0.0;
    return gd_b - gd_a; /* samples */
}

/* ================================================================== */
/* Kaiser parameter estimation (Oppenheim & Schafer pp. 475-476)       */
/* ================================================================== */

double fce_kaiser_beta(double atten_db)
{
    double a = fabs(atten_db);
    if (a > 50.0)
        return 0.1102 * (a - 8.7);
    if (a > 21.0)
        return 0.5842 * pow(a - 21.0, 0.4) + 0.07886 * (a - 21.0);
    return 0.0;
}

uint32_t fce_kaiser_taps(double atten_db, double transition_hz, double fs)
{
    double a = fabs(atten_db);
    double dw = 2.0 * FCE_PI * transition_hz / fs; /* rad/sample */
    double n;
    uint32_t taps;
    if (dw <= 0.0 || fs <= 0.0)
        return 0;
    n = (a - 7.95) / (2.285 * dw) + 1.0;
    if (n < 1.0)
        n = 1.0;
    taps = (uint32_t)ceil(n);
    /* keep odd so the design is Type-I/II symmetric */
    if ((taps & 1u) == 0u)
        taps += 1u;
    return taps;
}
