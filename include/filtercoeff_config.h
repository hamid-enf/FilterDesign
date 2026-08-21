/*
 * filtercoeff_config.h - Feature / size configuration for FilterCoeff.
 *
 * Every macro below has a safe default. Define any of them (e.g. in your
 * build system with -DFCE_ENABLE_FIR=0, or by editing this file) to shrink
 * the library to the features your project actually needs.
 *
 * The core never allocates memory dynamically. All sizes are bounded by the
 * FCE_MAX_* limits below; the required workspace is computed by
 * fce_workspace_required().
 */

#ifndef FILTERCOEFF_CONFIG_H
#define FILTERCOEFF_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* Feature switches                                                    */
/* ------------------------------------------------------------------ */

/* FIR coefficient generation (windowed-sinc + Kaiser + ...) */
#ifndef FCE_ENABLE_FIR
#define FCE_ENABLE_FIR 1
#endif

/* IIR coefficient generation (Butterworth/Chebyshev/Elliptic/Bessel) */
#ifndef FCE_ENABLE_IIR
#define FCE_ENABLE_IIR 1
#endif

/* Bessel family (requires a polynomial root finder; disable to save code) */
#ifndef FCE_ENABLE_IIR_BESSEL
#define FCE_ENABLE_IIR_BESSEL 1
#endif

/* Fixed-point conversion (Q15 / Q31) */
#ifndef FCE_ENABLE_FIXED_POINT
#define FCE_ENABLE_FIXED_POINT 1
#endif

#ifndef FCE_ENABLE_Q15
#define FCE_ENABLE_Q15 1
#endif

#ifndef FCE_ENABLE_Q31
#define FCE_ENABLE_Q31 1
#endif

/* Validation engine (stability, frequency response, error metrics) */
#ifndef FCE_ENABLE_VALIDATION
#define FCE_ENABLE_VALIDATION 1
#endif

/* Light-weight signal simulation used only to sanity-check coefficients */
#ifndef FCE_ENABLE_SIMULATION
#define FCE_ENABLE_SIMULATION 1
#endif

/* Exporters: C arrays, CSV, JSON, Markdown report */
#ifndef FCE_ENABLE_EXPORT
#define FCE_ENABLE_EXPORT 1
#endif

/* FILE*-based export writer (needs <stdio.h>; host only, optional) */
#ifndef FCE_ENABLE_EXPORT_STDIO
#define FCE_ENABLE_EXPORT_STDIO 1
#endif

/* ------------------------------------------------------------------ */
/* Size limits (workspace sizing)                                      */
/* ------------------------------------------------------------------ */

/* Maximum number of FIR taps supported by the workspace layout. */
#ifndef FCE_MAX_FIR_TAPS
#define FCE_MAX_FIR_TAPS 2048
#endif

/* Maximum IIR order (analog prototype order). */
#ifndef FCE_MAX_IIR_ORDER
#define FCE_MAX_IIR_ORDER 32
#endif

/* Maximum number of SOS sections: (order+1)/2. */
#ifndef FCE_MAX_SECTIONS
#define FCE_MAX_SECTIONS ((FCE_MAX_IIR_ORDER + 1) / 2)
#endif

/* Automatic-order guard: reject specs that would need an order above this. */
#ifndef FCE_MAX_AUTO_ORDER
#define FCE_MAX_AUTO_ORDER 24
#endif

/* Frequency-response grid used by validation / normalization scans. */
#ifndef FCE_VALIDATE_GRID_POINTS
#define FCE_VALIDATE_GRID_POINTS 512
#endif

/* ------------------------------------------------------------------ */
/* Numeric guards                                                      */
/* ------------------------------------------------------------------ */

/* Minimum pole-radius margin (1 - |p|) considered "safe". */
#ifndef FCE_STABILITY_MARGIN_MIN
#define FCE_STABILITY_MARGIN_MIN 1e-9
#endif

/* Relative epsilon for conjugate-pair detection in zpk2sos. */
#ifndef FCE_CPLX_TOL
#define FCE_CPLX_TOL 1e-9
#endif

#ifdef __cplusplus
}
#endif

#endif /* FILTERCOEFF_CONFIG_H */
