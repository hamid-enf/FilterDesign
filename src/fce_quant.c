/*
 * fce_quant.c - fixed-point conversion (Q15 / Q31).
 *
 * This is NOT a plain cast. Pipeline per strategy:
 *
 *   float coefficients -> determine scale -> range analysis
 *        -> quantization (round-to-nearest) -> saturation/overflow check
 *        -> quantization error metrics
 *
 * Strategies:
 *   FCE_SCALE_SYMMETRIC       one scale for the whole coefficient set
 *   FCE_SCALE_SECTION_WISE    one scale per SOS section
 *   FCE_SCALE_COEFFICIENT_WISE one scale per coefficient (storage only;
 *                              every coefficient saturates the format)
 *
 * Quantization:  q = round(c * scale),  scale = (2^F - 1) / max|c|
 * Reconstruction: c_tilde = q / scale
 * The reported "effective Q format" is Qm.F with m = ceil(log2(max|c|))
 * integer bits (0 when max|c| <= 1).
 */
#include "fce_internal.h"

#if FCE_ENABLE_FIXED_POINT

static int64_t fce_round_half(double x)
{
    return (int64_t)(x >= 0.0 ? floor(x + 0.5) : ceil(x - 0.5));
}

typedef struct fce_quant_cfg
{
    int      frac_bits;   /* 15 or 31 */
    int64_t  qmax;        /* 2^F - 1 */
    double   qmax_d;
    uint32_t sec_len;     /* coefficients per section (n or 5) */
    uint32_t n_sections;
} fce_quant_cfg_t;

static void fce_quant_apply(const double* c, uint32_t n,
                            const fce_quant_cfg_t* cfg,
                            const double* scales, /* per coefficient */
                            int16_t* q15, int32_t* q31,
                            double* max_abs, double* rms, double* max_rel,
                            uint32_t* overflow_cnt)
{
    uint32_t i;
    double sse = 0.0;
    double mx = 0.0;
    double mr = 0.0;

    if (overflow_cnt)
        *overflow_cnt = 0;

    for (i = 0; i < n; i++)
    {
        double v = c[i];
        double s = scales[i];
        int64_t q;
        double err;

        if (s > 0.0)
            q = fce_round_half(v * s);
        else
            q = 0;

        /* saturation / overflow check: count it HERE - the stored value
         * is clamped, so inspecting the array later can never see it */
        if (q > cfg->qmax)
        {
            q = cfg->qmax;
            if (overflow_cnt)
                (*overflow_cnt)++;
        }
        else if (q < -cfg->qmax)
        {
            q = -cfg->qmax;
            if (overflow_cnt)
                (*overflow_cnt)++;
        }

        if (q15)
            q15[i] = (int16_t)q;
        if (q31)
            q31[i] = (int32_t)q;

        /* quantization error in coefficient units */
        {
            double vt = (s > 0.0) ? (double)q / s : 0.0;
            err = fabs(v - vt);
            if (err > mx)
                mx = err;
            sse += err * err;
            if (v != 0.0)
            {
                double rel = err / fabs(v);
                if (rel > mr)
                    mr = rel;
            }
        }
    }

    *max_abs = mx;
    *rms = (n > 0u) ? sqrt(sse / (double)n) : 0.0;
    *max_rel = mr;
}

/*
 * Quantize `n` coefficients. `sec_len` = coefficients per scaling unit
 * (n for FIR / symmetric; 5 for SOS section-wise; 1 for coefficient-wise).
 * Fills the caller's arrays and error metrics.
 */
static fce_status_t fce_quant_core(const double* c, uint32_t n,
                                   fce_qformat_t qf,
                                   fce_scale_strategy_t strat,
                                   uint32_t sec_len,
                                   int16_t* q15, int32_t* q31,
                                   double* scales /* scratch, >= n */,
                                   double* sym_scale,
                                   double* sec_scales,
                                   double* coeff_scales,
                                   uint8_t* int_bits,
                                   double* max_abs, double* rms,
                                   double* max_rel,
                                   uint32_t* overflow_cnt)
{
    fce_quant_cfg_t cfg;
    uint32_t n_sec = (sec_len > 0u) ? (n + sec_len - 1u) / sec_len : 1u;
    uint32_t i, sec;

    if (n == 0u || c == NULL || scales == NULL)
        return FCE_ERR_INVALID_ARGUMENT;

    cfg.frac_bits = (qf == FCE_QFORMAT_Q15) ? 15 : 31;
    cfg.qmax = ((int64_t)1 << cfg.frac_bits) - 1;
    cfg.qmax_d = (double)cfg.qmax;
    cfg.sec_len = sec_len;
    cfg.n_sections = n_sec;

    *overflow_cnt = 0;

    switch (strat)
    {
    case FCE_SCALE_SYMMETRIC:
    {
        double mx = 0.0;
        double scale;
        for (i = 0; i < n; i++)
            if (fabs(c[i]) > mx)
                mx = fabs(c[i]);
        if (mx <= 0.0)
            scale = 1.0;
        else
            scale = cfg.qmax_d / mx;
        for (i = 0; i < n; i++)
            scales[i] = scale;
        if (sym_scale)
            *sym_scale = scale;
        if (int_bits)
        {
            if (mx > 1.0)
                *int_bits = (uint8_t)ceil(log2(mx));
            else
                *int_bits = 0;
        }
        if (sec_scales)
            for (sec = 0; sec < n_sec; sec++)
                sec_scales[sec] = scale;
        break;
    }

    case FCE_SCALE_SECTION_WISE:
    {
        for (sec = 0; sec < n_sec; sec++)
        {
            double mx = 0.0;
            double scale;
            uint32_t base = sec * sec_len;
            uint32_t cnt = (base + sec_len <= n) ? sec_len : (n - base);
            for (i = 0; i < cnt; i++)
                if (fabs(c[base + i]) > mx)
                    mx = fabs(c[base + i]);
            if (mx <= 0.0)
                scale = 1.0;
            else
                scale = cfg.qmax_d / mx;
            for (i = 0; i < cnt; i++)
                scales[base + i] = scale;
            if (sec_scales)
                sec_scales[sec] = scale;
        }
        if (sym_scale)
            *sym_scale = 0.0; /* not applicable */
        if (int_bits)
            *int_bits = 0;
        break;
    }

    case FCE_SCALE_COEFFICIENT_WISE:
    {
        double mx = 0.0;
        for (i = 0; i < n; i++)
        {
            double a = fabs(c[i]);
            if (a > mx)
                mx = a;
            if (a > 0.0)
                scales[i] = cfg.qmax_d / a;
            else
                scales[i] = 1.0;
            if (coeff_scales)
                coeff_scales[i] = scales[i];
        }
        if (sym_scale)
            *sym_scale = 0.0;
        if (int_bits)
        {
            if (mx > 1.0)
                *int_bits = (uint8_t)ceil(log2(mx));
            else
                *int_bits = 0;
        }
        break;
    }

    default:
        return FCE_ERR_INVALID_ARGUMENT;
    }

    fce_quant_apply(c, n, &cfg, scales,
                    qf == FCE_QFORMAT_Q15 ? q15 : NULL,
                    qf == FCE_QFORMAT_Q31 ? q31 : NULL,
                    max_abs, rms, max_rel, overflow_cnt);

    return FCE_OK;
}

fce_status_t fce_quant_fir(const double* h, uint32_t n,
                           fce_qformat_t qf, fce_scale_strategy_t strat,
                           int16_t* q15, int32_t* q31,
                           double* scale_scratch /* >= n */,
                           double* scale, double* coeff_scales,
                           uint8_t* int_bits,
                           double* max_abs, double* rms, double* max_rel,
                           uint32_t* overflow_cnt)
{
    return fce_quant_core(h, n, qf, strat, n, q15, q31, scale_scratch,
                          scale, NULL, coeff_scales, int_bits, max_abs,
                          rms, max_rel, overflow_cnt);
}

fce_status_t fce_quant_sos(const double* sos, uint32_t n_sections,
                           fce_qformat_t qf, fce_scale_strategy_t strat,
                           int16_t* q15, int32_t* q31,
                           double* scale_scratch /* >= 5*n_sections */,
                           double* scale, double* sec_scales,
                           double* coeff_scales, uint8_t* int_bits,
                           double* max_abs, double* rms, double* max_rel,
                           uint32_t* overflow_cnt)
{
    fce_status_t st;
    uint32_t sec_len = (strat == FCE_SCALE_SECTION_WISE) ? 5u
                       : (strat == FCE_SCALE_COEFFICIENT_WISE) ? 1u
                       : (5u * n_sections);
    st = fce_quant_core(sos, 5u * n_sections, qf, strat, sec_len,
                        q15, q31, scale_scratch, scale, sec_scales,
                        coeff_scales, int_bits, max_abs, rms, max_rel,
                        overflow_cnt);
    return st;
}

#endif /* FCE_ENABLE_FIXED_POINT */
