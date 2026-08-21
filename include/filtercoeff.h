/*
 * filtercoeff.h - Digital Filter Coefficient Generator.
 *
 * A professional, embedded-friendly C99 library that turns a filter
 * specification into validated, ready-to-use coefficients:
 *
 *     Specification -> Calculation -> Validation -> Float/Fixed-Point -> Export
 *
 * This library ONLY generates coefficients. It does not run filters.
 * Pair it with your existing FIR/IIR runtime library.
 *
 * Key properties:
 *   - C99, no dynamic allocation, no C++, no exceptions
 *   - internal math in float64 (double) for maximum accuracy
 *   - optional final float32 / Q15 / Q31 output
 *   - dependency: <stdint.h> <stddef.h> <stdbool.h> <math.h>
 *   - works on STM32H7, ESP32-C2, PC, ...
 *
 * Sign convention (IIR, SOS):
 *
 *     y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
 *
 * SOS array layout per section: { b0, b1, b2, a1, a2 }  (a0 == 1 implicit)
 */

#ifndef FILTERCOEFF_H
#define FILTERCOEFF_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include "filtercoeff_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================== */
/* Version                                                            */
/* ================================================================== */

#define FCE_VERSION_MAJOR 1
#define FCE_VERSION_MINOR 0
#define FCE_VERSION_PATCH 0
#define FCE_VERSION_STRING "1.0.0"

/* ================================================================== */
/* Error codes                                                        */
/* ================================================================== */

typedef enum fce_status
{
    FCE_OK = 0,                    /* success                                */
    FCE_ERR_INVALID_ARGUMENT,      /* NULL pointer / bad enum / bad combo    */
    FCE_ERR_INVALID_SPEC,          /* spec values are not usable             */
    FCE_ERR_UNSUPPORTED,           /* feature disabled or not implemented    */
    FCE_ERR_NUMERICAL,             /* internal math did not converge / NaN   */
    FCE_ERR_UNSTABLE,              /* designed filter is unstable            */
    FCE_ERR_OVERFLOW,              /* coefficient overflow (fixed point)     */
    FCE_ERR_QUANTIZATION,          /* quantization produced unusable result  */
    FCE_ERR_BUFFER_TOO_SMALL,      /* caller workspace/export buffer too big */
    FCE_ERR_NOT_AVAILABLE          /* requested data was not produced        */
} fce_status_t;

/* Human readable description of a status code. */
const char* fce_status_str(fce_status_t status);

/* ================================================================== */
/* Result flags (bitmask, warnings are NOT fatal)                      */
/* ================================================================== */

enum fce_flag
{
    FCE_FLAG_NONE                 = 0,
    FCE_FLAG_NUMERICAL_WARNING    = (1u << 0), /* precision/numerical risk   */
    FCE_FLAG_COEFFICIENT_OVERFLOW = (1u << 1), /* fixed-point overflow seen  */
    FCE_FLAG_UNSTABLE             = (1u << 2), /* float design unstable      */
    FCE_FLAG_QUANTIZATION_UNSTABLE= (1u << 3), /* quantized IIR unstable     */
    FCE_FLAG_QUANTIZATION_WARNING = (1u << 4), /* large quantization error   */
    FCE_FLAG_ORDER_CLAMPED        = (1u << 5), /* auto order/taps clamped    */
    FCE_FLAG_SPEC_MARGINAL        = (1u << 6), /* spec near numeric limits   */
    FCE_FLAG_SYMMETRY_WARNING     = (1u << 7)  /* FIR symmetry limitation    */
};

/* ================================================================== */
/* Filter kinds / types                                                */
/* ================================================================== */

typedef enum fce_kind
{
    FCE_KIND_FIR = 0,
    FCE_KIND_IIR
} fce_kind_t;

typedef enum fce_fir_type
{
    FCE_FIR_LOWPASS = 0,
    FCE_FIR_HIGHPASS,
    FCE_FIR_BANDPASS,
    FCE_FIR_BANDSTOP,
    FCE_FIR_HILBERT,        /* 90-degree phase shifter (Type III)     */
    FCE_FIR_DIFFERENTIATOR  /* discrete differentiator (Type IV)      */
} fce_fir_type_t;

typedef enum fce_iir_type
{
    FCE_IIR_LOWPASS = 0,
    FCE_IIR_HIGHPASS,
    FCE_IIR_BANDPASS,
    FCE_IIR_BANDSTOP
} fce_iir_type_t;

typedef enum fce_iir_family
{
    FCE_IIR_BUTTERWORTH = 0,
    FCE_IIR_CHEBYSHEV1,     /* passband ripple                         */
    FCE_IIR_CHEBYSHEV2,     /* stopband ripple (inverse Chebyshev)     */
    FCE_IIR_ELLIPTIC,       /* Cauer: ripple in both bands             */
    FCE_IIR_BESSEL          /* maximally flat group delay              */
} fce_iir_family_t;

/* FIR linear-phase symmetry types. */
typedef enum fce_fir_symmetry
{
    FCE_SYMMETRY_NONE = 0,  /* non-symmetric (not produced by this lib) */
    FCE_SYMMETRY_I,         /* odd  taps, symmetric                     */
    FCE_SYMMETRY_II,        /* even taps, symmetric (Nyquist null)      */
    FCE_SYMMETRY_III,       /* odd  taps, anti-symmetric (Hilbert)      */
    FCE_SYMMETRY_IV         /* even taps, anti-symmetric (differentiator)*/
} fce_fir_symmetry_t;

/* ================================================================== */
/* Design options                                                      */
/* ================================================================== */

typedef enum fce_window
{
    FCE_WIN_RECTANGULAR = 0,
    FCE_WIN_HANN,
    FCE_WIN_HAMMING,
    FCE_WIN_BLACKMAN,
    FCE_WIN_KAISER,
    FCE_WIN_BLACKMAN_HARRIS,
    FCE_WIN_BARTLETT,
    FCE_WIN_TUKEY
} fce_window_t;

typedef enum fce_precision
{
    FCE_PRECISION_FLOAT32 = 0,  /* internal float64, output float32 */
    FCE_PRECISION_FLOAT64       /* internal float64, output float64 */
} fce_precision_t;

/* Coefficient normalization strategy (what "unity gain" means). */
typedef enum fce_norm
{
    FCE_NORM_AUTO = 0,      /* LP: DC, HP: Nyquist, BP: passband peak,   */
                            /* BS: DC, Hilbert: in-band peak,            */
                            /* differentiator: Nyquist                   */
    FCE_NORM_DC,            /* |H(0)|        = 1                         */
    FCE_NORM_NYQUIST,       /* |H(pi*fs)|    = 1                         */
    FCE_NORM_PASSBAND_PEAK, /* max |H| inside passband = 1               */
    FCE_NORM_NONE           /* no gain scaling                           */
} fce_norm_t;

typedef enum fce_qformat
{
    FCE_QFORMAT_NONE = 0,   /* float only                                */
    FCE_QFORMAT_Q15,        /* 15 fractional bits (int16)                */
    FCE_QFORMAT_Q31         /* 31 fractional bits (int32)                */
} fce_qformat_t;

/* Fixed-point scaling strategy. */
typedef enum fce_scale_strategy
{
    FCE_SCALE_SYMMETRIC = 0,       /* one scale for the whole coefficient set  */
    FCE_SCALE_SECTION_WISE,        /* per-SOS-section scale (IIR)              */
    FCE_SCALE_COEFFICIENT_WISE     /* per-coefficient scale (storage only)     */
} fce_scale_strategy_t;

/* SOS section ordering strategy (IIR). */
typedef enum fce_sos_order
{
    FCE_SOS_ORDER_DEFAULT = 0,     /* pole radius ascending ("worst last"),    */
                                   /* matches scipy zpk2sos default            */
    FCE_SOS_ORDER_POLE_RADIUS_ASC, /* same as DEFAULT, explicit                */
    FCE_SOS_ORDER_POLE_RADIUS_DESC,/* high-Q sections first                    */
    FCE_SOS_ORDER_INTERNAL_GAIN    /* by measured section peak gain            */
} fce_sos_order_t;

/* Validation subset selector (bitmask). 0 = FCE_VALIDATE_DEFAULT. */
enum fce_validate
{
    FCE_VALIDATE_NONE      = 0,
    FCE_VALIDATE_STABILITY = (1u << 0), /* pole radii, stability margin        */
    FCE_VALIDATE_RESPONSE  = (1u << 1), /* response scans (gains, ripple, ...) */
    FCE_VALIDATE_SPEC      = (1u << 2), /* spec vs measured attenuation        */
    FCE_VALIDATE_QUANT     = (1u << 3), /* quantization error + stability      */
    FCE_VALIDATE_DEFAULT   = (1u << 0) | (1u << 1) | (1u << 2) | (1u << 3)
};

/* ================================================================== */
/* Specification                                                       */
/* ================================================================== */

typedef struct fce_spec
{
    fce_kind_t kind;                 /* FCE_KIND_FIR or FCE_KIND_IIR */

    /* --- common --- */
    double fs;                       /* sample rate [Hz], > 0           */
    double fc1;                      /* cutoff / first edge [Hz]        */
    double fc2;                      /* second edge [Hz] (BP/BS only)   */
    fce_precision_t precision;       /* output precision                */
    fce_norm_t normalization;        /* gain normalization strategy     */
    fce_qformat_t qformat;           /* Q15 / Q31 / none                */
    fce_scale_strategy_t scale_strategy;
    uint32_t validate;               /* bitmask of fce_validate, 0=default */

    /* --- FIR --- */
    fce_fir_type_t fir_type;
    uint16_t num_taps;               /* 0 = automatic (Kaiser)          */
    fce_window_t window;
    double kaiser_beta;              /* 0 = automatic from attenuation  */
    double transition_hz;            /* for Kaiser auto-taps            */

    /* --- IIR --- */
    fce_iir_family_t iir_family;
    fce_iir_type_t iir_type;
    uint16_t order;                  /* 0 = automatic                   */
    double passband_ripple_db;       /* gpass: cheby1/ellip (+auto)     */
    double stopband_atten_db;        /* gstop: cheby2/ellip (+auto)     */
    double edge1_hz;                 /* opposite-band edge 1 (auto-order only) */
    double edge2_hz;                 /* opposite-band edge 2 (auto-order only) */
    fce_sos_order_t sos_order;
} fce_spec_t;

/* ================================================================== */
/* Result                                                              */
/* ================================================================== */

typedef struct fce_result
{
    /* overall */
    fce_status_t status;             /* FCE_OK or first error            */
    uint32_t flags;                  /* fce_flag bitmask                 */
    fce_kind_t kind;

    /* design parameters actually used */
    fce_fir_type_t fir_type;
    fce_iir_family_t iir_family;
    fce_iir_type_t iir_type;
    fce_window_t window;
    double kaiser_beta;              /* effective Kaiser beta (0 if N/A) */
    double fs, fc1, fc2;             /* user frequencies                 */
    double design_fc1, design_fc2;   /* frequencies actually used        */
    double passband_ripple_db, stopband_atten_db;
    double edge1_hz, edge2_hz;       /* opposite-band edges (auto-order) */
    double transition_hz;            /* transition band (Kaiser auto)    */
    uint16_t num_taps;               /* FIR tap count                    */
    uint16_t order;                  /* IIR prototype order              */
    uint16_t num_sections;           /* IIR SOS section count            */
    fce_fir_symmetry_t symmetry;     /* FIR symmetry type                */
    fce_precision_t precision;
    fce_norm_t normalization;
    double norm_factor;              /* gain multiplier applied          */
    fce_qformat_t qformat;
    fce_scale_strategy_t scale_strategy;
    fce_sos_order_t sos_order;

    /* --- float coefficients (point into the workspace) --- */
    const float*  h_f32;             /* FIR, length num_taps             */
    const double* h_f64;             /* FIR, length num_taps             */
    const float*  sos_f32;           /* IIR, length 5*num_sections       */
    const double* sos_f64;           /* IIR, length 5*num_sections       */

    /* --- FIR internals (transparency) --- */
    const double* fir_ideal;         /* ideal (pre-window) response      */
    const double* fir_window;        /* window values                    */

    /* --- fixed-point coefficients --- */
    const int16_t* q15;              /* FIR: num_taps; IIR: 5*num_sections */
    const int32_t* q31;
    double scale;                    /* symmetric scale (q = round(c*scale)) */
    const double* section_scales;    /* section-wise, length num_sections  */
    const double* coeff_scales;      /* coefficient-wise                   */
    uint8_t q_int_bits;              /* effective integer bits (0 = Q0.F)  */

    /* --- fixed-point error metrics (coefficient domain) --- */
    double q_max_abs_error;
    double q_rms_error;
    double q_max_rel_error;

    /* --- validation --- */
    double stability_margin;         /* min(1 - |pole|) over sections    */
    double max_pole_radius;
    double dc_gain_db;
    double nyquist_gain_db;
    double passband_ripple_measured_db;
    double stopband_atten_measured_db;
    double cutoff_measured_hz;       /* -3 dB (or -gpass) crossing       */
    double quant_response_max_error_db; /* float vs quantized response   */
    double quant_max_pole_radius;    /* quantized IIR pole radius        */
    double quant_stability_margin;   /* quantized IIR stability margin   */

    /* --- IIR internals (transparency / no-black-box) --- */
    const double* section_gains;     /* peak gain per section, length S  */
    uint16_t sos_pairs_valid;        /* sections kept after validation   */
    const double* iir_poles;         /* digital poles (re,im pairs)      */
    const double* iir_zeros;         /* digital zeros (re,im pairs)      */
    uint16_t iir_npoles;
    uint16_t iir_nzeros;

    /* workspace bookkeeping */
    size_t workspace_size;           /* bytes actually used              */
} fce_result_t;

/* ================================================================== */
/* Workspace                                                           */
/* ================================================================== */

/*
 * Caller-provided memory. No dynamic allocation anywhere in the library.
 * The coefficient arrays in fce_result_t point INTO this workspace, so it
 * must stay alive as long as the coefficients are used.
 */
typedef struct fce_workspace
{
    void* data;
    size_t size;
} fce_workspace_t;

/*
 * Number of bytes needed to design the given spec.
 * Returns 0 if the spec is invalid.
 */
size_t fce_workspace_required(const fce_spec_t* spec);

/*
 * Recommended workspace for an arbitrary/unknown spec
 * (sized for the maximum supported taps/order).
 */
size_t fce_workspace_required_max(void);

/* ================================================================== */
/* Main API                                                            */
/* ================================================================== */

/*
 * Design a filter from `spec`.
 *
 *   - fills `result` (metadata + coefficient pointers into ws->data)
 *   - runs validation according to spec->validate
 *   - performs fixed-point conversion if spec->qformat != FCE_QFORMAT_NONE
 *
 * Returns FCE_OK on success. Detailed diagnostics live in result->flags.
 */
fce_status_t fce_generate(const fce_spec_t* spec,
                          fce_result_t* result,
                          fce_workspace_t* ws);

/*
 * Fill a spec with defaults (all zeros; kind unset).
 * Use the fce_spec_*_... constructors below for convenience.
 */
void fce_spec_defaults(fce_spec_t* spec);

/* ================================================================== */
/* Convenience spec constructors (simple API)                          */
/* ================================================================== */

static inline void fce_spec_fir(fce_spec_t* sp,
                                fce_fir_type_t type,
                                double fs,
                                double fc1,
                                double fc2,
                                uint16_t num_taps,
                                fce_window_t win,
                                double kaiser_beta,
                                fce_precision_t prec)
{
    fce_spec_defaults(sp);
    sp->kind       = FCE_KIND_FIR;
    sp->fir_type   = type;
    sp->fs         = fs;
    sp->fc1        = fc1;
    sp->fc2        = fc2;
    sp->num_taps   = num_taps;
    sp->window     = win;
    sp->kaiser_beta = kaiser_beta;
    sp->precision  = prec;
}

static inline void fce_spec_iir(fce_spec_t* sp,
                                fce_iir_family_t family,
                                fce_iir_type_t type,
                                double fs,
                                double fc1,
                                double fc2,
                                uint16_t order,
                                double ripple_db,
                                double atten_db,
                                fce_precision_t prec)
{
    fce_spec_defaults(sp);
    sp->kind      = FCE_KIND_IIR;
    sp->iir_family = family;
    sp->iir_type  = type;
    sp->fs        = fs;
    sp->fc1       = fc1;
    sp->fc2       = fc2;
    sp->order     = order;
    sp->passband_ripple_db = ripple_db;
    sp->stopband_atten_db  = atten_db;
    sp->precision = prec;
}

/* ================================================================== */
/* Validation engine (also available standalone)                        */
/* ================================================================== */

#if FCE_ENABLE_VALIDATION

/* A single point of a frequency-response scan. */
typedef struct fce_response_point
{
    double f_hz;        /* frequency [Hz]                     */
    double mag;         /* linear magnitude |H|               */
    double mag_db;      /* 20*log10(mag)                      */
    double phase_deg;   /* phase [degrees]                    */
    double group_delay; /* group delay [samples]              */
} fce_response_point_t;

/* Callback receiving each scan point. Return false to abort the scan. */
typedef bool (*fce_response_cb)(void* ctx, const fce_response_point_t* pt);

/*
 * Scan the frequency response of an FIR filter (h, n taps) or an IIR SOS
 * filter (sos, 5 coefficients per section) on a grid of n_points points
 * from f_start to f_stop [Hz].
 */
fce_status_t fce_response_fir(const double* h, uint16_t n,
                              double fs, uint32_t n_points,
                              double f_start, double f_stop,
                              fce_response_cb cb, void* ctx);

fce_status_t fce_response_sos(const double* sos, uint16_t n_sections,
                              double fs, uint32_t n_points,
                              double f_start, double f_stop,
                              fce_response_cb cb, void* ctx);

/*
 * Stability check of an SOS filter. Returns FCE_OK if every pole radius
 * is strictly below 1. `max_radius` and `margin` are always written.
 */
fce_status_t fce_stability_sos(const double* sos, uint16_t n_sections,
                               double* max_radius, double* margin);

/*
 * Error metrics between two coefficient sets (float reference vs
 * quantized, or generated vs ideal): max abs, RMS, max relative.
 * n = number of coefficients.
 */
void fce_coeff_error(const double* ref, const double* test, uint32_t n,
                     double* max_abs, double* rms, double* max_rel);

#endif /* FCE_ENABLE_VALIDATION */

/* ================================================================== */
/* Simulation (validation aid ONLY - not a filter runtime)             */
/* ================================================================== */

#if FCE_ENABLE_SIMULATION

/*
 * Minimal direct-form FIR application and transposed direct-form-II SOS
 * application. These exist purely to verify generated coefficients on a
 * test signal; they are NOT a production runtime.
 */
void fce_sim_fir(const double* h, uint16_t n,
                 const double* x, double* y, uint32_t nx);
void fce_sim_sos(const double* sos, uint16_t n_sections,
                 const double* x, double* y, uint32_t nx,
                 double* state /* 2*n_sections, zeroed by caller */);

/* Test signal generators. Returns number of samples written. */
typedef enum fce_sim_signal
{
    FCE_SIM_SINE = 0,
    FCE_SIM_MULTITONE,
    FCE_SIM_IMPULSE,
    FCE_SIM_STEP,
    FCE_SIM_WHITE_NOISE,
    FCE_SIM_CHIRP,
    FCE_SIM_DC_PLUS_NOISE
} fce_sim_signal_t;

/*
 * Generate a test signal. `a` and `b` are generator arguments:
 *   SINE:          a = frequency [Hz],        b = amplitude
 *   MULTITONE:     a = f1 [Hz],               b = f2 [Hz]  (three unity tones: f1, f2, 2*f1+f2)
 *   IMPULSE:       a = sample index of impulse, b = amplitude
 *   STEP:          a = step start sample,     b = amplitude
 *   WHITE_NOISE:   a = seed,                  b = amplitude (uniform [-b, b])
 *   CHIRP:         a = f_start [Hz],          b = f_stop [Hz]
 *   DC_PLUS_NOISE: a = DC level,              b = noise amplitude
 */
uint32_t fce_sim_signal(fce_sim_signal_t kind, double a, double b,
                        double fs, double* out, uint32_t n);

#endif /* FCE_ENABLE_SIMULATION */

/* ================================================================== */
/* Export                                                              */
/* ================================================================== */

#if FCE_ENABLE_EXPORT

/* Writer abstraction: the export functions stream text through `write`. */
typedef size_t (*fce_write_fn)(void* ctx, const char* data, size_t len);

typedef struct fce_writer
{
    fce_write_fn write;
    void* ctx;
} fce_writer_t;

/* Memory-backed writer (embedded friendly). */
typedef struct fce_mem_writer
{
    char*  buf;
    size_t size;
    size_t pos;
    int    truncated;
} fce_mem_writer_t;

void fce_writer_mem_init(fce_writer_t* w, fce_mem_writer_t* mw,
                         char* buf, size_t size);
/* FILE*-backed writer (host; needs FCE_ENABLE_EXPORT_STDIO). */
fce_status_t fce_writer_file_init(fce_writer_t* w, void* file);

/* Options shared by all exporters. */
typedef struct fce_export_opts
{
    const char* name;      /* C identifier for the array / object name     */
    uint32_t    precision; /* decimal digits (0 = auto)                    */
    bool        include_quantized; /* include Q15/Q31 arrays if available  */
} fce_export_opts_t;

/* C header/source arrays. */
fce_status_t fce_export_c_fir(const fce_result_t* r, const fce_export_opts_t* o,
                              fce_writer_t* w);
fce_status_t fce_export_c_sos(const fce_result_t* r, const fce_export_opts_t* o,
                              fce_writer_t* w);

/* CSV (one row per coefficient; FIR or SOS). */
fce_status_t fce_export_csv(const fce_result_t* r, fce_writer_t* w);

/* JSON document with metadata + arrays. */
fce_status_t fce_export_json(const fce_result_t* r, fce_writer_t* w);

/* Human-readable Markdown report. */
fce_status_t fce_export_report(const fce_result_t* r, fce_writer_t* w);

#endif /* FCE_ENABLE_EXPORT */

/* ================================================================== */
/* Simple-API compatibility aliases                                    */
/* ================================================================== */

typedef fce_spec_t     FilterCoeffSpec;
typedef fce_result_t   FilterCoeffResult;
typedef fce_workspace_t FilterCoeffWorkspace;
typedef fce_status_t   FilterCoeffStatus;

static inline fce_status_t FilterCoeff_Generate(const FilterCoeffSpec* spec,
                                                FilterCoeffResult* result,
                                                FilterCoeffWorkspace* ws)
{
    return fce_generate(spec, result, ws);
}

#ifdef __cplusplus
}
#endif

#endif /* FILTERCOEFF_H */
