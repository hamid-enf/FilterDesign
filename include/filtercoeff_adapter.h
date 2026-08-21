/*
 * filtercoeff_adapter.h - optional adapter layer (NOT part of the core).
 *
 * The core is library-agnostic. These tiny inline helpers only repack the
 * coefficient arrays into the layout expected by a few common runtimes,
 * so you never have to convert formats by hand.
 *
 *   - scipy / MATLAB sos:  [b0, b1, b2, a0, a1, a2]   (a0 = 1)
 *   - CMSIS-DSP biquad:    b0, b1, b2, a1, a2        (same as ours)
 *   - generic FIR:         plain coefficient array    (same as ours)
 */
#ifndef FILTERCOEFF_ADAPTER_H
#define FILTERCOEFF_ADAPTER_H

#include "filtercoeff.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* scipy / MATLAB style SOS: add the explicit a0 = 1 column.           */
/* dst must hold 6 * num_sections doubles/floats.                      */
/* ------------------------------------------------------------------ */
static inline void fce_adapter_sos_to_scipy(const fce_result_t* r,
                                            double* dst)
{
    uint32_t s, j;
    for (s = 0; s < r->num_sections; s++)
    {
        for (j = 0; j < 5u; j++)
            dst[6u * s + (j < 3u ? j : j + 1u)] = r->sos_f64[5u * s + j];
        dst[6u * s + 3u] = 1.0; /* a0 */
    }
}

static inline void fce_adapter_sos_to_scipy_f32(const fce_result_t* r,
                                                float* dst)
{
    uint32_t s, j;
    for (s = 0; s < r->num_sections; s++)
    {
        for (j = 0; j < 5u; j++)
            dst[6u * s + (j < 3u ? j : j + 1u)] = r->sos_f32[5u * s + j];
        dst[6u * s + 3u] = 1.0f;
    }
}

/* ------------------------------------------------------------------ */
/* Q15/Q31 SOS with explicit a0:  dst[6*s + 3] = 2^F (a0 in the same   */
/* fixed-point scale).                                                 */
/* ------------------------------------------------------------------ */
static inline void fce_adapter_sos_q15_to_6(const fce_result_t* r,
                                            int16_t* dst)
{
    uint32_t s, j;
    for (s = 0; s < r->num_sections; s++)
    {
        for (j = 0; j < 5u; j++)
            dst[6u * s + (j < 3u ? j : j + 1u)] = r->q15[5u * s + j];
        dst[6u * s + 3u] = 32767; /* a0 = 1.0 in Q15 */
    }
}

/* ------------------------------------------------------------------ */
/* per-section fixed-point scale as a power-of-two shift (Q format):   */
/* returns the shift `sh` such that  c_q15 * 2^-sh ~= c / sec_scale.   */
/* (Useful for runtimes that apply per-section gains by shifting.)     */
/* ------------------------------------------------------------------ */
static inline int fce_adapter_sec_shift_bits(const fce_result_t* r,
                                             uint32_t section,
                                             double* actual_scale)
{
    double s = (r->scale_strategy == FCE_SCALE_SECTION_WISE)
                   ? r->section_scales[section]
                   : r->scale;
    int sh = (int)floor(log2(s) + 0.5);
    if (sh < 0)
        sh = 0;
    if (actual_scale)
        *actual_scale = ldexp(1.0, sh);
    return sh;
}

#ifdef __cplusplus
}
#endif

#endif /* FILTERCOEFF_ADAPTER_H */
