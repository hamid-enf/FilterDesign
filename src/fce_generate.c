/*
 * fce_generate.c - public API orchestration.
 *
 *   spec validation -> workspace layout -> auto taps/order
 *        -> FIR/IIR design -> quantization -> validation -> result
 *
 * The workspace is caller-provided; the library never allocates memory.
 */
#include "fce_internal.h"

/* ================================================================== */
/* misc                                                                */
/* ================================================================== */

const char* fce_status_str(fce_status_t status)
{
    switch (status)
    {
    case FCE_OK: return "OK";
    case FCE_ERR_INVALID_ARGUMENT: return "invalid argument";
    case FCE_ERR_INVALID_SPEC: return "invalid specification";
    case FCE_ERR_UNSUPPORTED: return "unsupported feature";
    case FCE_ERR_NUMERICAL: return "numerical error";
    case FCE_ERR_UNSTABLE: return "unstable filter";
    case FCE_ERR_OVERFLOW: return "coefficient overflow";
    case FCE_ERR_QUANTIZATION: return "quantization error";
    case FCE_ERR_BUFFER_TOO_SMALL: return "buffer too small";
    case FCE_ERR_NOT_AVAILABLE: return "data not available";
    default: return "unknown";
    }
}

void fce_spec_defaults(fce_spec_t* spec)
{
    if (spec == NULL)
        return;
    memset(spec, 0, sizeof(*spec));
    spec->precision = FCE_PRECISION_FLOAT32;
    spec->normalization = FCE_NORM_AUTO;
    spec->scale_strategy = FCE_SCALE_SYMMETRIC;
    spec->sos_order = FCE_SOS_ORDER_DEFAULT;
    spec->validate = FCE_VALIDATE_DEFAULT;
    spec->window = FCE_WIN_HAMMING;
}

/* ================================================================== */
/* workspace layout                                                    */
/* ================================================================== */

bool fce_layout_compute(const fce_spec_t* spec, fce_layout_t* lay)
{
    size_t off = 0;
    uint32_t n_taps = 0, order = 0, sections;

    memset(lay, 0, sizeof(*lay));

    if (spec->kind == FCE_KIND_FIR)
    {
        n_taps = spec->num_taps;
        if (n_taps == 0u)
        {
            /* auto (Kaiser): estimate an upper bound for sizing */
            uint32_t est = 0;
            if (spec->window == FCE_WIN_KAISER &&
                spec->stopband_atten_db > 0.0 &&
                spec->transition_hz > 0.0 && spec->fs > 0.0)
                est = fce_kaiser_taps(spec->stopband_atten_db,
                                      spec->transition_hz, spec->fs);
            n_taps = (est > 0u) ? est : FCE_MAX_FIR_TAPS;
        }
        if (n_taps == 0u || n_taps > FCE_MAX_FIR_TAPS)
            return false;
        lay->n_taps = n_taps;

        off = 0;
        lay->off_fir_work  = off; off += 2u * (size_t)n_taps * sizeof(double);
        lay->off_fir_h     = fce_align8(off); off = lay->off_fir_h + (size_t)n_taps * sizeof(double);
        lay->off_fir_h32   = fce_align8(off); off = lay->off_fir_h32 + (size_t)n_taps * sizeof(float);
        lay->off_fir_q15   = fce_align8(off); off = lay->off_fir_q15 + (size_t)n_taps * sizeof(int16_t);
        lay->off_fir_q31   = fce_align8(off); off = lay->off_fir_q31 + (size_t)n_taps * sizeof(int32_t);
        lay->off_fir_scales= fce_align8(off); off = lay->off_fir_scales + (size_t)n_taps * sizeof(double);
        lay->total = fce_align8(off);
    }
    else if (spec->kind == FCE_KIND_IIR)
    {
        order = spec->order;
        if (order == 0u)
            order = FCE_MAX_AUTO_ORDER; /* upper bound for sizing */
        if (order == 0u || order > FCE_MAX_IIR_ORDER)
            return false;
        sections = order; /* BP/BS double the number of sections */
        lay->order = order;
        lay->n_sections = sections;

        off = 0;
        lay->off_proto_p  = off; off += 2u * (size_t)order * sizeof(double);
        lay->off_proto_z  = fce_align8(off); off = lay->off_proto_z + 2u * (size_t)order * sizeof(double);
        lay->off_ap       = fce_align8(off); off = lay->off_ap + 4u * (size_t)order * sizeof(double);
        lay->off_az       = fce_align8(off); off = lay->off_az + 4u * (size_t)order * sizeof(double);
        lay->off_dp       = fce_align8(off); off = lay->off_dp + 4u * (size_t)order * sizeof(double);
        lay->off_dz       = fce_align8(off); off = lay->off_dz + 4u * (size_t)order * sizeof(double);
        lay->off_sos      = fce_align8(off); off = lay->off_sos + 5u * (size_t)sections * sizeof(double);
        lay->off_sos32    = fce_align8(off); off = lay->off_sos32 + 5u * (size_t)sections * sizeof(float);
        lay->off_sos_q15  = fce_align8(off); off = lay->off_sos_q15 + 5u * (size_t)sections * sizeof(int16_t);
        lay->off_sos_q31  = fce_align8(off); off = lay->off_sos_q31 + 5u * (size_t)sections * sizeof(int32_t);
        lay->off_sec_scales = fce_align8(off); off = lay->off_sec_scales + (size_t)sections * sizeof(double);
        lay->off_sec_gains  = fce_align8(off); off = lay->off_sec_gains + (size_t)sections * sizeof(double);
        lay->off_coeff_scales = fce_align8(off); off = lay->off_coeff_scales + 5u * (size_t)sections * sizeof(double);
        lay->off_scratch  = fce_align8(off); off = lay->off_scratch + 8u * (size_t)order * sizeof(double);
        lay->total = fce_align8(off);
    }
    else
    {
        return false;
    }
    return true;
}

size_t fce_workspace_required(const fce_spec_t* spec)
{
    fce_layout_t lay;
    if (spec == NULL)
        return 0u;
    if (!fce_layout_compute(spec, &lay))
        return 0u;
    return lay.total;
}

size_t fce_workspace_required_max(void)
{
    fce_spec_t sp;
    fce_layout_t lay;
    size_t a = 0, b = 0;
    memset(&sp, 0, sizeof(sp));
    sp.kind = FCE_KIND_FIR;
    sp.num_taps = FCE_MAX_FIR_TAPS;
    if (fce_layout_compute(&sp, &lay))
        a = lay.total;
    memset(&sp, 0, sizeof(sp));
    sp.kind = FCE_KIND_IIR;
    sp.order = FCE_MAX_IIR_ORDER;
    if (fce_layout_compute(&sp, &lay))
        b = lay.total;
    return (a > b) ? a : b;
}

/* ================================================================== */
/* quantization pipeline                                               */
/* ================================================================== */

#if FCE_ENABLE_FIXED_POINT

static fce_status_t fce_pipeline_quant(fce_result_t* r, fce_layout_t* lay,
                                       void* base)
{
    fce_status_t st;
    uint32_t overflow = 0;
    double max_abs = 0.0, rms = 0.0, max_rel = 0.0;

    if (r->qformat == FCE_QFORMAT_NONE)
        return FCE_OK;

    if (r->kind == FCE_KIND_FIR)
    {
        int16_t* q15 = (int16_t*)(void*)((char*)base + lay->off_fir_q15);
        int32_t* q31 = (int32_t*)(void*)((char*)base + lay->off_fir_q31);
        double* cscales = (double*)(void*)((char*)base + lay->off_fir_scales);
        st = fce_quant_fir(r->h_f64, r->num_taps, r->qformat,
                           r->scale_strategy, q15, q31, cscales,
                           &r->scale, cscales,
                           &r->q_int_bits, &max_abs, &rms, &max_rel,
                           &overflow);
        if (st != FCE_OK)
            return st;
        if (r->qformat == FCE_QFORMAT_Q15)
            r->q15 = q15;
        else
            r->q31 = q31;
        if (r->scale_strategy == FCE_SCALE_COEFFICIENT_WISE)
            r->coeff_scales = cscales;
    }
    else
    {
        int16_t* q15 = (int16_t*)(void*)((char*)base + lay->off_sos_q15);
        int32_t* q31 = (int32_t*)(void*)((char*)base + lay->off_sos_q31);
        double* sscales = (double*)(void*)((char*)base + lay->off_sec_scales);
        double* cscales = (double*)(void*)((char*)base + lay->off_coeff_scales);
        st = fce_quant_sos(r->sos_f64, r->num_sections, r->qformat,
                           r->scale_strategy, q15, q31, cscales,
                           &r->scale, sscales,
                           cscales, &r->q_int_bits, &max_abs, &rms, &max_rel,
                           &overflow);
        if (st != FCE_OK)
            return st;
        if (r->qformat == FCE_QFORMAT_Q15)
            r->q15 = q15;
        else
            r->q31 = q31;
        if (r->scale_strategy == FCE_SCALE_SECTION_WISE)
            r->section_scales = sscales;
        if (r->scale_strategy == FCE_SCALE_COEFFICIENT_WISE)
            r->coeff_scales = cscales;
    }

    r->q_max_abs_error = max_abs;
    r->q_rms_error = rms;
    r->q_max_rel_error = max_rel;
    if (overflow > 0u)
    {
        r->flags |= FCE_FLAG_COEFFICIENT_OVERFLOW;
        return FCE_ERR_OVERFLOW;
    }

#if FCE_ENABLE_VALIDATION
    /* response error float vs quantized + stability of quantized IIR */
    {
        double* recon = (double*)(void*)((char*)base + lay->off_scratch);
        uint32_t n;
        uint32_t i;

        if (r->kind == FCE_KIND_FIR)
        {
            n = r->num_taps;
            for (i = 0; i < n; i++)
            {
                double q = (r->q15 != NULL) ? (double)r->q15[i]
                                            : (double)r->q31[i];
                double s = (r->scale_strategy == FCE_SCALE_COEFFICIENT_WISE)
                               ? r->coeff_scales[i] : r->scale;
                recon[i] = (s > 0.0) ? q / s : 0.0;
            }
            r->quant_response_max_error_db =
                fce_validate_quant_response(r->h_f64, recon, n, 0, r->fs);
        }
        else
        {
            uint32_t ns = r->num_sections;
            n = 5u * ns;
            for (i = 0; i < n; i++)
            {
                double q = (r->q15 != NULL) ? (double)r->q15[i]
                                            : (double)r->q31[i];
                double s;
                if (r->scale_strategy == FCE_SCALE_COEFFICIENT_WISE)
                    s = r->coeff_scales[i];
                else if (r->scale_strategy == FCE_SCALE_SECTION_WISE)
                    s = r->section_scales[i / 5u];
                else
                    s = r->scale;
                recon[i] = (s > 0.0) ? q / s : 0.0;
            }
            r->quant_response_max_error_db =
                fce_validate_quant_response(r->sos_f64, recon, ns, 1, r->fs);

            /* stability of the quantized filter */
            {
                double mr = 0.0;
                for (i = 0; i < ns; i++)
                {
                    double a1 = recon[5u * i + 3];
                    double a2 = recon[5u * i + 4];
                    fce_cplx_t p1, p2;
                    double r1, r2;
                    fce_biquad_poles(a1, a2, &p1, &p2);
                    r1 = fce_cx_abs(p1);
                    r2 = fce_cx_abs(p2);
                    if (r1 > mr) mr = r1;
                    if (r2 > mr) mr = r2;
                }
                r->quant_max_pole_radius = mr;
                r->quant_stability_margin = 1.0 - mr;
                if (mr >= 1.0)
                {
                    r->flags |= FCE_FLAG_QUANTIZATION_UNSTABLE;
                    return FCE_ERR_QUANTIZATION;
                }
            }
        }
    }
#endif /* FCE_ENABLE_VALIDATION */

    return FCE_OK;
}

#endif /* FCE_ENABLE_FIXED_POINT */

/* ================================================================== */
/* main entry                                                          */
/* ================================================================== */

fce_status_t fce_generate(const fce_spec_t* spec,
                          fce_result_t* result,
                          fce_workspace_t* ws)
{
    fce_layout_t lay;
    fce_status_t st = FCE_OK;
    fce_status_t hard = FCE_OK;
    fce_spec_t sp;
    void* base;

    if (spec == NULL || result == NULL || ws == NULL)
        return FCE_ERR_INVALID_ARGUMENT;

    memset(result, 0, sizeof(*result));
    sp = *spec; /* local copy */

    if (ws->data == NULL || ws->size == 0u)
    {
        result->status = FCE_ERR_INVALID_ARGUMENT;
        return result->status;
    }
    base = ws->data;

    /* ---------- spec validation ---------- */
    if (sp.kind != FCE_KIND_FIR && sp.kind != FCE_KIND_IIR)
        hard = FCE_ERR_INVALID_SPEC;
    if (sp.fs <= 0.0)
        hard = FCE_ERR_INVALID_SPEC;
    if (sp.precision != FCE_PRECISION_FLOAT32 &&
        sp.precision != FCE_PRECISION_FLOAT64)
        hard = FCE_ERR_INVALID_SPEC;
    if (sp.normalization < FCE_NORM_AUTO || sp.normalization > FCE_NORM_NONE)
        hard = FCE_ERR_INVALID_SPEC;
#if FCE_ENABLE_FIXED_POINT
    if (sp.qformat != FCE_QFORMAT_NONE &&
        sp.qformat != FCE_QFORMAT_Q15 && sp.qformat != FCE_QFORMAT_Q31)
        hard = FCE_ERR_INVALID_SPEC;
#else
    if (sp.qformat != FCE_QFORMAT_NONE)
        hard = FCE_ERR_UNSUPPORTED;
#endif

    if (sp.kind == FCE_KIND_FIR)
    {
#if FCE_ENABLE_FIR
        if (sp.fir_type < FCE_FIR_LOWPASS || sp.fir_type > FCE_FIR_DIFFERENTIATOR)
            hard = FCE_ERR_INVALID_SPEC;
        /* Hilbert / differentiator may be full-band (fc1 = 0) */
        if (sp.fir_type == FCE_FIR_HILBERT ||
            sp.fir_type == FCE_FIR_DIFFERENTIATOR)
        {
            if (!(sp.fc1 >= 0.0) || !(sp.fc1 < 0.5 * sp.fs))
                hard = FCE_ERR_INVALID_SPEC;
        }
        else if (!(sp.fc1 > 0.0) || !(sp.fc1 < 0.5 * sp.fs))
        {
            hard = FCE_ERR_INVALID_SPEC;
        }
        if ((sp.fir_type == FCE_FIR_BANDPASS ||
             sp.fir_type == FCE_FIR_BANDSTOP) &&
            (!(sp.fc2 > sp.fc1) || !(sp.fc2 < 0.5 * sp.fs)))
            hard = FCE_ERR_INVALID_SPEC;
        if (sp.num_taps == 0u)
        {
            if (sp.window != FCE_WIN_KAISER ||
                !(sp.stopband_atten_db > 0.0) ||
                !(sp.transition_hz > 0.0) ||
                !(sp.transition_hz < 0.5 * sp.fs))
                hard = FCE_ERR_INVALID_SPEC;
        }
        else if (sp.num_taps > FCE_MAX_FIR_TAPS)
        {
            hard = FCE_ERR_BUFFER_TOO_SMALL;
        }
        if (sp.num_taps == 1u)
            hard = FCE_ERR_INVALID_SPEC;
        if (sp.window < FCE_WIN_RECTANGULAR || sp.window > FCE_WIN_TUKEY)
            hard = FCE_ERR_INVALID_SPEC;
#else
        hard = FCE_ERR_UNSUPPORTED;
#endif
    }
    else if (sp.kind == FCE_KIND_IIR)
    {
#if FCE_ENABLE_IIR
        if (sp.iir_family < FCE_IIR_BUTTERWORTH || sp.iir_family > FCE_IIR_BESSEL)
            hard = FCE_ERR_INVALID_SPEC;
        if (sp.iir_type < FCE_IIR_LOWPASS || sp.iir_type > FCE_IIR_BANDSTOP)
            hard = FCE_ERR_INVALID_SPEC;
        if (!(sp.fc1 > 0.0) || !(sp.fc1 < 0.5 * sp.fs))
            hard = FCE_ERR_INVALID_SPEC;
        if ((sp.iir_type == FCE_IIR_BANDPASS ||
             sp.iir_type == FCE_IIR_BANDSTOP) &&
            (!(sp.fc2 > sp.fc1) || !(sp.fc2 < 0.5 * sp.fs)))
            hard = FCE_ERR_INVALID_SPEC;
        if (sp.iir_family == FCE_IIR_CHEBYSHEV1 ||
            sp.iir_family == FCE_IIR_ELLIPTIC)
        {
            if (!(sp.passband_ripple_db > 0.0))
                hard = FCE_ERR_INVALID_SPEC;
        }
        if (sp.iir_family == FCE_IIR_CHEBYSHEV2 ||
            sp.iir_family == FCE_IIR_ELLIPTIC)
        {
            if (!(sp.stopband_atten_db > 0.0))
                hard = FCE_ERR_INVALID_SPEC;
        }
        if (sp.order > FCE_MAX_IIR_ORDER)
            hard = FCE_ERR_BUFFER_TOO_SMALL;
#else
        hard = FCE_ERR_UNSUPPORTED;
#endif
    }

    if (hard != FCE_OK)
    {
        result->status = hard;
        return hard;
    }

    /* ---------- workspace ---------- */
    if (!fce_layout_compute(&sp, &lay))
    {
        result->status = FCE_ERR_BUFFER_TOO_SMALL;
        return result->status;
    }
    if (ws->size < lay.total)
    {
        result->status = FCE_ERR_BUFFER_TOO_SMALL;
        result->workspace_size = lay.total;
        return result->status;
    }
    result->workspace_size = lay.total;

    /* ---------- shared metadata ---------- */
    result->kind = sp.kind;
    result->fs = sp.fs;
    result->fc1 = sp.fc1;
    result->fc2 = sp.fc2;
    result->precision = sp.precision;
    result->normalization = sp.normalization;
    result->qformat = sp.qformat;
    result->scale_strategy = sp.scale_strategy;
    result->sos_order = sp.sos_order;

    /* ---------- design ---------- */
    if (sp.kind == FCE_KIND_FIR)
    {
#if FCE_ENABLE_FIR
        uint32_t taps = sp.num_taps;
        int auto_taps = (taps == 0u);
        if (auto_taps)
        {
            taps = fce_kaiser_taps(sp.stopband_atten_db, sp.transition_hz,
                                   sp.fs);
            if (taps > FCE_MAX_FIR_TAPS)
            {
                taps = FCE_MAX_FIR_TAPS;
                result->flags |= FCE_FLAG_ORDER_CLAMPED;
            }
            sp.num_taps = (uint16_t)taps;
            lay.n_taps = taps;
        }
        result->fir_type = sp.fir_type;
        result->window = sp.window;
        result->transition_hz = sp.transition_hz;
        result->stopband_atten_db = sp.stopband_atten_db;
        result->design_fc1 = sp.fc1;
        result->design_fc2 = sp.fc2;

        st = fce_fir_design(&sp, result, &lay, base);
        if (st == FCE_OK)
            result->num_taps = (uint16_t)taps;

        if (st == FCE_OK && auto_taps)
        {
            /* Kaiser auto-taps spec check: stopband attenuation (LP) */
            if (sp.fir_type == FCE_FIR_LOWPASS)
            {
                double f0 = sp.fc1 + sp.transition_hz;
                double mn, mx;
                if (f0 < 0.5 * sp.fs)
                {
                    fce_scan_fir_band(result->h_f64, result->num_taps,
                                      sp.fs, f0, 0.5 * sp.fs, &mn, &mx);
                    if (-mx < sp.stopband_atten_db - 1.0)
                        result->flags |= FCE_FLAG_SPEC_MARGINAL;
                }
            }
        }
#else
        st = FCE_ERR_UNSUPPORTED;
#endif
    }
    else
    {
#if FCE_ENABLE_IIR
        uint32_t order = sp.order;
        int auto_order = (order == 0u);
        result->iir_family = sp.iir_family;
        result->iir_type = sp.iir_type;
        result->passband_ripple_db = sp.passband_ripple_db;
        result->stopband_atten_db = sp.stopband_atten_db;
        result->edge1_hz = sp.edge1_hz;
        result->edge2_hz = sp.edge2_hz;
        result->transition_hz = sp.transition_hz;

        if (auto_order)
        {
            fce_auto_t auto_info;
            st = fce_iir_auto_order(&sp, &auto_info);
            if (st == FCE_OK)
            {
                if (auto_info.order > FCE_MAX_AUTO_ORDER)
                {
                    result->flags |= FCE_FLAG_ORDER_CLAMPED;
                }
                order = auto_info.order;
                sp.order = (uint16_t)order;
                lay.order = order;
                lay.n_sections = order;
                result->design_fc1 = auto_info.design_fc1;
                result->design_fc2 = auto_info.design_fc2;
            }
        }
        else
        {
            result->design_fc1 = sp.fc1;
            result->design_fc2 = sp.fc2;
        }

        if (st == FCE_OK)
            st = fce_iir_design(&sp, result, &lay, base);

        if (st == FCE_OK && auto_order && sp.stopband_atten_db > 0.0)
        {
            /* spec check: attenuation at the stopband edge */
            double f_edge = sp.edge1_hz;
            if (f_edge > 0.0 && f_edge < 0.5 * sp.fs)
            {
                double w = 2.0 * FCE_PI * f_edge / sp.fs;
                double g = 1.0;
                uint32_t s;
                for (s = 0; s < result->num_sections; s++)
                {
                    fce_cplx_t h = fce_eval_biquad(result->sos_f64 + 5u * s, w);
                    g *= fce_cx_abs(h);
                }
                if (g > 0.0)
                {
                    double att = -20.0 * log10(g);
                    if (att < sp.stopband_atten_db - 1.0)
                        result->flags |= FCE_FLAG_SPEC_MARGINAL;
                }
            }
        }
#else
        st = FCE_ERR_UNSUPPORTED;
#endif
    }

    if (st != FCE_OK)
    {
        result->status = st;
        return st;
    }

    /* ---------- fixed point ---------- */
#if FCE_ENABLE_FIXED_POINT
    if (result->qformat != FCE_QFORMAT_NONE)
    {
        st = fce_pipeline_quant(result, &lay, base);
        if (st != FCE_OK)
        {
            result->status = st;
            return st;
        }
    }
#endif

    /* ---------- validation ---------- */
#if FCE_ENABLE_VALIDATION
    if (sp.validate != FCE_VALIDATE_NONE)
    {
        if (result->kind == FCE_KIND_IIR)
        {
            double mr = 0.0, margin = 0.0;
            fce_status_t vst = fce_stability_sos(result->sos_f64,
                                                 result->num_sections,
                                                 &mr, &margin);
            result->max_pole_radius = mr;
            result->stability_margin = margin;
            if (vst != FCE_OK)
            {
                result->flags |= FCE_FLAG_UNSTABLE;
                result->status = FCE_ERR_UNSTABLE;
                return result->status;
            }
            if (margin < FCE_STABILITY_MARGIN_MIN)
                result->flags |= FCE_FLAG_SPEC_MARGINAL;

            fce_validate_iir_measures(&sp, result, result->sos_f64,
                                      result->num_sections);
        }
        else
        {
            fce_validate_fir_measures(&sp, result, result->h_f64,
                                      result->num_taps);
        }
    }
#endif /* FCE_ENABLE_VALIDATION */

    result->status = FCE_OK;
    return FCE_OK;
}
