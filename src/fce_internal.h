/*
 * fce_internal.h - private helpers shared between translation units.
 * Not installed; not part of the public API.
 */
#ifndef FCE_INTERNAL_H
#define FCE_INTERNAL_H

#include "filtercoeff.h"

#include <math.h>
#include <string.h>

/* PI constant (M_PI is not portable under strict -std=c99) */
#ifndef FCE_PI
#define FCE_PI 3.14159265358979323846264338327950288
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* complex double                                                      */
/* ------------------------------------------------------------------ */

typedef struct fce_cplx
{
    double re;
    double im;
} fce_cplx_t;

static inline fce_cplx_t fce_cx(double re, double im)
{
    fce_cplx_t z;
    z.re = re;
    z.im = im;
    return z;
}

static inline fce_cplx_t fce_cx_add(fce_cplx_t a, fce_cplx_t b)
{
    return fce_cx(a.re + b.re, a.im + b.im);
}

static inline fce_cplx_t fce_cx_sub(fce_cplx_t a, fce_cplx_t b)
{
    return fce_cx(a.re - b.re, a.im - b.im);
}

static inline fce_cplx_t fce_cx_mul(fce_cplx_t a, fce_cplx_t b)
{
    return fce_cx(a.re * b.re - a.im * b.im,
                  a.re * b.im + a.im * b.re);
}

static inline fce_cplx_t fce_cx_div(fce_cplx_t a, fce_cplx_t b)
{
    double d = b.re * b.re + b.im * b.im;
    return fce_cx((a.re * b.re + a.im * b.im) / d,
                  (a.im * b.re - a.re * b.im) / d);
}

static inline fce_cplx_t fce_cx_scale(fce_cplx_t a, double s)
{
    return fce_cx(a.re * s, a.im * s);
}

static inline fce_cplx_t fce_cx_conj(fce_cplx_t a)
{
    return fce_cx(a.re, -a.im);
}

static inline double fce_cx_abs(fce_cplx_t a)
{
    return hypot(a.re, a.im);
}

static inline double fce_cx_abs_sq(fce_cplx_t a)
{
    return a.re * a.re + a.im * a.im;
}

/* complex sqrt / ln / asin / atanh (principal branches) */
fce_cplx_t fce_cx_sqrt(fce_cplx_t a);
fce_cplx_t fce_cx_ln(fce_cplx_t a);
fce_cplx_t fce_cx_asin(fce_cplx_t a);
fce_cplx_t fce_cx_atanh(fce_cplx_t a);
int fce_cx_isfinite(fce_cplx_t a);

/* ------------------------------------------------------------------ */
/* math helpers                                                        */
/* ------------------------------------------------------------------ */

/* sinc(x) = sin(pi*x)/(pi*x), sinc(0) = 1  (numpy convention) */
double fce_sinc(double x);

/* modified Bessel I0 (series) */
double fce_i0(double x);

/* complete elliptic integral K(m), 0 <= m < 1 (AGM) */
double fce_ellipk(double m);
/* K(1-m) */
double fce_ellipkm1(double m);

/* Jacobi elliptic functions sn/cn/dn for real u, 0<=m<=1 (AGM/Landen) */
void fce_ellipj(double u, double m, double* sn, double* cn, double* dn);

/* inverse Jacobi sn: solve w = sn(z, m) for complex w (ascending Landen) */
fce_cplx_t fce_arc_jac_sn(fce_cplx_t w, double m);

/* solve the elliptic degree equation: n*K(m)/K(1-m) = K(m1)/K(1-m1) */
double fce_ellipdeg(uint32_t n, double m1);

/* Durand-Kerner roots of a real-coefficient polynomial.
 * coeff[k] = coefficient of x^k, k = 0..deg. Roots written to roots[].
 * Returns FCE_OK / FCE_ERR_NUMERICAL. */
fce_status_t fce_roots_durand_kerner(const double* coeff, uint32_t deg,
                                     fce_cplx_t* roots);
fce_status_t fce_roots_durand_kerner_start(const double* coeff, uint32_t deg,
                                           fce_cplx_t* roots,
                                           const fce_cplx_t* start);

/* Bessel prototype poles (delay normalization), 'mag' renormalization done
 * by the caller. */
fce_status_t fce_bessel_poles(uint32_t n, fce_cplx_t* poles);

/* ------------------------------------------------------------------ */
/* workspace layout                                                    */
/* ------------------------------------------------------------------ */

typedef struct fce_layout
{
    /* FIR */
    size_t off_fir_work;    /* double[2*N]: ideal + windowed             */
    size_t off_fir_h;       /* double[N] final float64 coefficients      */
    size_t off_fir_h32;     /* float[N]                                  */
    size_t off_fir_q15;     /* int16[N]                                  */
    size_t off_fir_q31;     /* int32[N]                                  */
    size_t off_fir_scales;  /* double[N] coefficient-wise scales         */
    /* IIR */
    size_t off_proto_p;     /* double[2*O]   analog prototype poles      */
    size_t off_proto_z;     /* double[2*O]   analog prototype zeros      */
    size_t off_ap;          /* double[4*O]   analog poles (transformed)  */
    size_t off_az;          /* double[4*O]   analog zeros                */
    size_t off_dp;          /* double[4*O]   digital poles               */
    size_t off_dz;          /* double[4*O]   digital zeros               */
    size_t off_sos;         /* double[5*S]   SOS float64                 */
    size_t off_sos32;       /* float[5*S]                                */
    size_t off_sos_q15;     /* int16[5*S]                                */
    size_t off_sos_q31;     /* int32[5*S]                                */
    size_t off_sec_scales;  /* double[S]                                 */
    size_t off_sec_gains;   /* double[S]  peak gain per section          */
    size_t off_coeff_scales;/* double[5*S]                               */
    size_t off_scratch;     /* scratch for pairing (double[4*O])         */
    size_t total;
    uint32_t n_taps;        /* effective FIR taps (0 = none)             */
    uint32_t order;         /* effective IIR order (0 = none)            */
    uint32_t n_sections;    /* (order+1)/2                               */
} fce_layout_t;

/* Compute the layout for a spec. Returns false if the spec can never fit. */
bool fce_layout_compute(const fce_spec_t* spec, fce_layout_t* lay);

/* Alignment helper: next offset aligned to 8 bytes. */
static inline size_t fce_align8(size_t off)
{
    return (off + 7u) & ~(size_t)7u;
}

/* Typed pointer helpers. */
static inline double* fce_lp_dbl(void* base, size_t off)
{
    return (double*)(void*)((char*)base + off);
}

static inline float* fce_lp_flt(void* base, size_t off)
{
    return (float*)(void*)((char*)base + off);
}

/* ------------------------------------------------------------------ */
/* shared numeric helpers                                              */
/* ------------------------------------------------------------------ */

/* window value at index n of an N-point window (shared with tests) */
double fce_window_value(fce_window_t win, uint32_t n, uint32_t N,
                        double kaiser_beta, double tukey_alpha);
/* transition half-width of a windowed design in Hz (main-lobe edge) */
double fce_window_transition_half_hz(fce_window_t win, double atten_db,
                                     double kaiser_beta,
                                     uint32_t num_taps, double fs);

/* Kaiser parameters from attenuation (Oppenheim & Schafer). */
double fce_kaiser_beta(double atten_db);
/* Kaiser tap estimate: transition_hz, fs, atten_db. Returns >= 1. */
uint32_t fce_kaiser_taps(double atten_db, double transition_hz, double fs);

/* prewarp: analog angular frequency for a digital frequency f [Hz] */
static inline double fce_prewarp(double f, double fs)
{
    return 2.0 * fs * tan(FCE_PI * f / fs);
}

/* Evaluate a biquad H(z) = (b0+b1 z^-1+b2 z^-2)/(1+a1 z^-1+a2 z^-2)
 * at z^-1 = e^{-jw}. */
fce_cplx_t fce_eval_biquad(const double* c /* b0,b1,b2,a1,a2 */, double w);

/* peak |H| of a biquad over [0, pi] using an adaptive scan */
double fce_biquad_peak_gain(const double* c, uint32_t grid);

/* poles of z^2 + a1*z + a2 (our convention) */
void fce_biquad_poles(double a1, double a2, fce_cplx_t* p1, fce_cplx_t* p2);

/* analytic group delay (samples) of a biquad at w: gd_a - gd_b */
double fce_biquad_group_delay(const double* c, double w);

/* band scans (validation internals, shared with fce_generate) */
void fce_scan_fir_band(const double* h, uint32_t n, double fs,
                       double f_lo, double f_hi,
                       double* min_db, double* max_db);
void fce_scan_sos_band(const double* sos, uint32_t ns, double fs,
                       double f_lo, double f_hi,
                       double* min_db, double* max_db);

/* ------------------------------------------------------------------ */
/* design module entry points (shared with fce_generate.c)             */
/* ------------------------------------------------------------------ */

fce_status_t fce_fir_design(const fce_spec_t* sp, fce_result_t* r,
                            fce_layout_t* lay, void* base);
fce_status_t fce_iir_design(const fce_spec_t* sp, fce_result_t* r,
                            fce_layout_t* lay, void* base);

typedef struct fce_auto
{
    uint32_t order;
    int      clamped;   /* 1 when the order was clamped to FCE_MAX_AUTO_ORDER */
    double design_fc1; /* Hz (LP/HP: cutoff; BP/BS: band edges) */
    double design_fc2;
} fce_auto_t;

fce_status_t fce_iir_auto_order(const fce_spec_t* sp, fce_auto_t* a);

/* ------------------------------------------------------------------ */
/* quantization (shared with fce_generate.c)                           */
/* ------------------------------------------------------------------ */

#if FCE_ENABLE_FIXED_POINT
fce_status_t fce_quant_fir(const double* h, uint32_t n,
                           fce_qformat_t qf, fce_scale_strategy_t strat,
                           int16_t* q15, int32_t* q31,
                           double* scale_scratch,
                           double* scale, double* coeff_scales,
                           uint8_t* int_bits,
                           double* max_abs, double* rms, double* max_rel,
                           uint32_t* overflow_cnt);
fce_status_t fce_quant_sos(const double* sos, uint32_t n_sections,
                           fce_qformat_t qf, fce_scale_strategy_t strat,
                           int16_t* q15, int32_t* q31,
                           double* scale_scratch,
                           double* scale, double* sec_scales,
                           double* coeff_scales, uint8_t* int_bits,
                           double* max_abs, double* rms, double* max_rel,
                           uint32_t* overflow_cnt);
#endif

/* ------------------------------------------------------------------ */
/* validation internals                                                */
/* ------------------------------------------------------------------ */

#if FCE_ENABLE_VALIDATION
void fce_validate_iir_measures(const fce_spec_t* sp, fce_result_t* r,
                               const double* sos, uint32_t ns);
void fce_validate_fir_measures(const fce_spec_t* sp, fce_result_t* r,
                               const double* h, uint32_t n);
double fce_validate_quant_response(const double* ref, const double* test,
                                   uint32_t n, int is_sos, double fs);
#endif

#ifdef __cplusplus
}
#endif

#endif /* FCE_INTERNAL_H */
