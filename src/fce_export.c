/*
 * fce_export.c - coefficient exporters: C arrays, CSV, JSON, Markdown report.
 *
 * All exporters stream text through a fce_writer_t so they work on any
 * target (memory buffer on MCU, FILE* on host, UART, ...).
 * Requires a working snprintf (available on every supported toolchain).
 */
#include "fce_internal.h"

#if FCE_ENABLE_EXPORT

#include <stdio.h>
#include <stdarg.h>

/* forward declarations (used by the int-array emitters below) */
static int fce_w_printf(fce_writer_t* w, const char* fmt, ...);
static int fce_w_str(fce_writer_t* w, const char* s);

/* dB values are clamped at 20log10(1e-300) = -6000 dB when the gain is
 * at/below the float64 floor (e.g. a filter with a zero at Nyquist).
 * Rendering that raw looks like a broken number; show a floor note. */
static int fce_w_db_row(fce_writer_t* w, const char* label, double db)
{
    if (db <= -6000.0 + 1e-3)
        return fce_w_printf(w, "| %s | < -600 dB |\n", label);
    return fce_w_printf(w, "| %s | %.4f dB |\n", label, db);
}

/* ------------------------------------------------------------------ */
/* writers                                                             */
/* ------------------------------------------------------------------ */

static int fce_emit_int_array16(fce_writer_t* w, const char* name,
                                const int16_t* v, uint32_t n)
{
    uint32_t i;
    if (!fce_w_printf(w, "static const int16_t %s[%u] = {\n", name,
                      (unsigned)n))
        return 0;
    for (i = 0; i < n; i++)
        if (!fce_w_printf(w, "    %d%s\n", (int)v[i],
                          (i + 1u < n) ? "," : ""))
            return 0;
    return fce_w_str(w, "};\n");
}

static int fce_emit_int_array32(fce_writer_t* w, const char* name,
                                const int32_t* v, uint32_t n)
{
    uint32_t i;
    if (!fce_w_printf(w, "static const int32_t %s[%u] = {\n", name,
                      (unsigned)n))
        return 0;
    for (i = 0; i < n; i++)
        if (!fce_w_printf(w, "    %ld%s\n", (long)v[i],
                          (i + 1u < n) ? "," : ""))
            return 0;
    return fce_w_str(w, "};\n");
}

size_t fce_mem_writer_write(void* ctx, const char* data, size_t len)
{
    fce_mem_writer_t* mw = (fce_mem_writer_t*)ctx;
    size_t room = (mw->pos < mw->size) ? (mw->size - mw->pos) : 0u;
    size_t n;
    /* keep one byte for the NUL terminator so the buffer is always a
     * valid C string (strstr-safe) */
    if (room > 0u)
        room -= 1u;
    n = (len < room) ? len : room;
    if (n > 0u)
    {
        memcpy(mw->buf + mw->pos, data, n);
        mw->pos += n;
    }
    if (mw->pos < mw->size)
        mw->buf[mw->pos] = '\0';
    if (n < len)
        mw->truncated = 1;
    return n;
}

void fce_writer_mem_init(fce_writer_t* w, fce_mem_writer_t* mw,
                         char* buf, size_t size)
{
    if (w == NULL || mw == NULL)
        return;
    mw->buf = buf;
    mw->size = size;
    mw->pos = 0u;
    mw->truncated = 0;
    w->write = fce_mem_writer_write;
    w->ctx = mw;
}

/* internal: formatted write; returns false on truncation */
static int fce_w_printf(fce_writer_t* w, const char* fmt, ...)
{
    char tmp[160];
    va_list ap;
    int n;
    va_start(ap, fmt);
    n = vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    if (n < 0)
        return 0;
    if ((size_t)n >= sizeof(tmp))
        n = (int)sizeof(tmp) - 1;
    return w->write(w->ctx, tmp, (size_t)n) == (size_t)n;
}

static int fce_w_str(fce_writer_t* w, const char* s)
{
    size_t n = strlen(s);
    return w->write(w->ctx, s, n) == n;
}

/* ------------------------------------------------------------------ */
/* shared helpers                                                      */
/* ------------------------------------------------------------------ */

static const char* fce_kind_str(fce_kind_t k)
{
    return (k == FCE_KIND_FIR) ? "FIR" : "IIR";
}

static const char* fce_fir_type_str(fce_fir_type_t t)
{
    switch (t)
    {
    case FCE_FIR_LOWPASS: return "Lowpass";
    case FCE_FIR_HIGHPASS: return "Highpass";
    case FCE_FIR_BANDPASS: return "Bandpass";
    case FCE_FIR_BANDSTOP: return "Bandstop";
    case FCE_FIR_HILBERT: return "Hilbert";
    case FCE_FIR_DIFFERENTIATOR: return "Differentiator";
    default: return "?";
    }
}

static const char* fce_iir_family_str(fce_iir_family_t f)
{
    switch (f)
    {
    case FCE_IIR_BUTTERWORTH: return "Butterworth";
    case FCE_IIR_CHEBYSHEV1: return "Chebyshev I";
    case FCE_IIR_CHEBYSHEV2: return "Chebyshev II";
    case FCE_IIR_ELLIPTIC: return "Elliptic";
    case FCE_IIR_BESSEL: return "Bessel";
    default: return "?";
    }
}

static const char* fce_iir_type_str(fce_iir_type_t t)
{
    switch (t)
    {
    case FCE_IIR_LOWPASS: return "Lowpass";
    case FCE_IIR_HIGHPASS: return "Highpass";
    case FCE_IIR_BANDPASS: return "Bandpass";
    case FCE_IIR_BANDSTOP: return "Bandstop";
    default: return "?";
    }
}

static const char* fce_window_str(fce_window_t win)
{
    switch (win)
    {
    case FCE_WIN_RECTANGULAR: return "Rectangular";
    case FCE_WIN_HANN: return "Hann";
    case FCE_WIN_HAMMING: return "Hamming";
    case FCE_WIN_BLACKMAN: return "Blackman";
    case FCE_WIN_KAISER: return "Kaiser";
    case FCE_WIN_BLACKMAN_HARRIS: return "Blackman-Harris";
    case FCE_WIN_BARTLETT: return "Bartlett";
    case FCE_WIN_TUKEY: return "Tukey";
    default: return "?";
    }
}

static const char* fce_norm_str(fce_norm_t n)
{
    switch (n)
    {
    case FCE_NORM_AUTO: return "Auto";
    case FCE_NORM_DC: return "Unity DC gain";
    case FCE_NORM_NYQUIST: return "Unity Nyquist gain";
    case FCE_NORM_PASSBAND_PEAK: return "Unity passband peak";
    case FCE_NORM_NONE: return "None";
    default: return "?";
    }
}

static const char* fce_qformat_str(fce_qformat_t q)
{
    switch (q)
    {
    case FCE_QFORMAT_NONE: return "none";
    case FCE_QFORMAT_Q15: return "Q15";
    case FCE_QFORMAT_Q31: return "Q31";
    default: return "?";
    }
}

static int fce_export_header(fce_writer_t* w, const fce_result_t* r)
{
    if (!fce_w_str(w, "/* Generated by FilterCoeff v" FCE_VERSION_STRING " */\n"))
        return 0;
    if (!fce_w_printf(w, "/* Filter: %s %s", fce_kind_str(r->kind),
                      (r->kind == FCE_KIND_FIR) ? fce_fir_type_str(r->fir_type)
                                                : fce_iir_type_str(r->iir_type)))
        return 0;
    if (r->kind == FCE_KIND_IIR)
    {
        if (!fce_w_printf(w, " (%s, order %u, sections %u)",
                          fce_iir_family_str(r->iir_family),
                          (unsigned)r->order, (unsigned)r->num_sections))
            return 0;
    }
    if (!fce_w_printf(w, " */\n"))
        return 0;
    if (!fce_w_printf(w, "/* Fs = %.1f Hz", r->fs))
        return 0;
    if (r->kind == FCE_KIND_FIR)
    {
        if (!fce_w_printf(w, ", Fc = %.4f Hz, Taps = %u, Window = %s",
                          r->fc1, (unsigned)r->num_taps,
                          fce_window_str(r->window)))
            return 0;
    }
    else
    {
        if (!fce_w_printf(w, ", Fc = %.4f Hz", r->fc1))
            return 0;
    }
    if (!fce_w_str(w, " */\n"))
        return 0;
    return 1;
}

static int fce_emit_double_array(fce_writer_t* w, const char* name,
                                 const char* type, const double* v,
                                 uint32_t n, uint32_t digits, int is_float)
{
    uint32_t i;
    if (!fce_w_printf(w, "static const %s %s[%u] = {\n", type, name,
                      (unsigned)n))
        return 0;
    for (i = 0; i < n; i++)
    {
        if (is_float)
        {
            if (!fce_w_printf(w, "    %.*gf%s\n", (int)digits, v[i],
                              (i + 1u < n) ? "," : ""))
                return 0;
        }
        else
        {
            if (!fce_w_printf(w, "    %.*g\n", (int)digits, v[i]))
                return 0;
        }
    }
    return fce_w_str(w, "};\n");
}

/* ------------------------------------------------------------------ */
/* C export                                                            */
/* ------------------------------------------------------------------ */

fce_status_t fce_export_c_fir(const fce_result_t* r, const fce_export_opts_t* o,
                              fce_writer_t* w)
{
    fce_export_opts_t def;
    const double* h;
    uint32_t n, digits;
    int ok;

    if (r == NULL || w == NULL || w->write == NULL)
        return FCE_ERR_INVALID_ARGUMENT;
    if (r->kind != FCE_KIND_FIR)
        return FCE_ERR_INVALID_SPEC;
    if (o == NULL)
    {
        memset(&def, 0, sizeof(def));
        def.name = "coeffs";
        def.include_quantized = 1;
        o = &def;
    }
    h = r->h_f64;
    n = r->num_taps;
    if (h == NULL || n == 0u)
        return FCE_ERR_NOT_AVAILABLE;
    digits = (o->precision > 0u) ? o->precision : 17u;

    ok = fce_export_header(w, r);
    if (r->precision == FCE_PRECISION_FLOAT32)
    {
        /* emit float32 array from the float copy when available */
        if (r->h_f32 != NULL)
        {
            uint32_t i;
            ok = ok && fce_w_printf(w,
                "static const float32_t %s[%u] = {\n", o->name, (unsigned)n);
            for (i = 0; i < n; i++)
                ok = ok && fce_w_printf(w, "    %.9gf%s\n",
                                        (double)r->h_f32[i],
                                        (i + 1u < n) ? "," : "");
            ok = ok && fce_w_str(w, "};\n");
        }
        else
        {
            ok = ok && fce_emit_double_array(w, o->name, "float32_t", h, n,
                                             digits > 9u ? 9u : digits, 1);
        }
    }
    else
    {
        ok = ok && fce_emit_double_array(w, o->name, "float64_t", h, n,
                                         digits, 0);
    }

    if (ok && o->include_quantized && r->q15 != NULL)
        ok = ok && fce_emit_int_array16(w, "coeffs_q15", r->q15, n);
    if (ok && o->include_quantized && r->q31 != NULL)
        ok = ok && fce_emit_int_array32(w, "coeffs_q31", r->q31, n);

    return ok ? FCE_OK : FCE_ERR_BUFFER_TOO_SMALL;
}

fce_status_t fce_export_c_sos(const fce_result_t* r, const fce_export_opts_t* o,
                              fce_writer_t* w)
{
    fce_export_opts_t def;
    const double* sos;
    uint32_t ns, digits;
    uint32_t i;
    int ok;

    if (r == NULL || w == NULL || w->write == NULL)
        return FCE_ERR_INVALID_ARGUMENT;
    if (r->kind != FCE_KIND_IIR)
        return FCE_ERR_INVALID_SPEC;
    if (o == NULL)
    {
        memset(&def, 0, sizeof(def));
        def.name = "sos";
        def.include_quantized = 1;
        o = &def;
    }
    sos = r->sos_f64;
    ns = r->num_sections;
    if (sos == NULL || ns == 0u)
        return FCE_ERR_NOT_AVAILABLE;
    digits = (o->precision > 0u) ? o->precision : 17u;

    ok = fce_export_header(w, r);
    ok = ok && fce_w_str(w, "/* SOS layout per section: { b0, b1, b2, a1, a2 } */\n");
    ok = ok && fce_w_str(w, "/* y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]");
    ok = ok && fce_w_str(w, " - a1*y[n-1] - a2*y[n-2] */\n");

    if (r->precision == FCE_PRECISION_FLOAT32)
    {
        ok = ok && fce_w_printf(w, "static const float32_t %s[][5] = {\n",
                                o->name);
        for (i = 0; i < ns; i++)
        {
            uint32_t j;
            ok = ok && fce_w_str(w, "    { ");
            for (j = 0; j < 5u; j++)
            {
                double v = (r->sos_f32 != NULL) ? (double)r->sos_f32[5u * i + j]
                                                : sos[5u * i + j];
                ok = ok && fce_w_printf(w, "%.9gf%s", v,
                                        (j + 1u < 5u) ? ", " : "");
            }
            ok = ok && fce_w_printf(w, " }%s\n", (i + 1u < ns) ? "," : "");
        }
        ok = ok && fce_w_str(w, "};\n");
    }
    else
    {
        ok = ok && fce_w_printf(w, "static const float64_t %s[][5] = {\n",
                                o->name);
        for (i = 0; i < ns; i++)
        {
            uint32_t j;
            ok = ok && fce_w_str(w, "    { ");
            for (j = 0; j < 5u; j++)
                ok = ok && fce_w_printf(w, "%.*g%s", (int)digits,
                                        sos[5u * i + j],
                                        (j + 1u < 5u) ? ", " : "");
            ok = ok && fce_w_printf(w, " }%s\n", (i + 1u < ns) ? "," : "");
        }
        ok = ok && fce_w_str(w, "};\n");
    }

    if (ok && o->include_quantized && r->q15 != NULL)
    {
        ok = ok && fce_w_printf(w, "static const int16_t %s_q15[][5] = {\n",
                                o->name);
        for (i = 0; i < ns; i++)
        {
            uint32_t j;
            ok = ok && fce_w_str(w, "    { ");
            for (j = 0; j < 5u; j++)
                ok = ok && fce_w_printf(w, "%d%s", (int)r->q15[5u * i + j],
                                        (j + 1u < 5u) ? ", " : "");
            ok = ok && fce_w_printf(w, " }%s\n", (i + 1u < ns) ? "," : "");
        }
        ok = ok && fce_w_str(w, "};\n");
    }
    if (ok && o->include_quantized && r->q31 != NULL)
    {
        ok = ok && fce_w_printf(w, "static const int32_t %s_q31[][5] = {\n",
                                o->name);
        for (i = 0; i < ns; i++)
        {
            uint32_t j;
            ok = ok && fce_w_str(w, "    { ");
            for (j = 0; j < 5u; j++)
                ok = ok && fce_w_printf(w, "%ld%s", (long)r->q31[5u * i + j],
                                        (j + 1u < 5u) ? ", " : "");
            ok = ok && fce_w_printf(w, " }%s\n", (i + 1u < ns) ? "," : "");
        }
        ok = ok && fce_w_str(w, "};\n");
    }

    return ok ? FCE_OK : FCE_ERR_BUFFER_TOO_SMALL;
}

/* ------------------------------------------------------------------ */
/* CSV                                                                 */
/* ------------------------------------------------------------------ */

fce_status_t fce_export_csv(const fce_result_t* r, fce_writer_t* w)
{
    int ok;
    uint32_t i, rows, cols;
    if (r == NULL || w == NULL || w->write == NULL)
        return FCE_ERR_INVALID_ARGUMENT;

    if (r->kind == FCE_KIND_FIR)
    {
        const double* h = r->h_f64;
        rows = r->num_taps;
        cols = 1u;
        if (h == NULL || rows == 0u)
            return FCE_ERR_NOT_AVAILABLE;
        ok = fce_w_str(w, "index,coefficient");
        if (r->q15) ok = ok && fce_w_str(w, ",q15");
        if (r->q31) ok = ok && fce_w_str(w, ",q31");
        ok = ok && fce_w_str(w, "\n");
        for (i = 0; i < rows; i++)
        {
            ok = ok && fce_w_printf(w, "%u,%.17g", (unsigned)i, h[i]);
            if (r->q15) ok = ok && fce_w_printf(w, ",%d", (int)r->q15[i]);
            if (r->q31) ok = ok && fce_w_printf(w, ",%ld", (long)r->q31[i]);
            ok = ok && fce_w_str(w, "\n");
        }
    }
    else
    {
        const double* sos = r->sos_f64;
        rows = r->num_sections;
        cols = 5u;
        if (sos == NULL || rows == 0u)
            return FCE_ERR_NOT_AVAILABLE;
        ok = fce_w_str(w, "section,b0,b1,b2,a1,a2");
        if (r->q15) ok = ok && fce_w_str(w, ",q15_b0,q15_b1,q15_b2,q15_a1,q15_a2");
        if (r->q31) ok = ok && fce_w_str(w, ",q31_b0,q31_b1,q31_b2,q31_a1,q31_a2");
        ok = ok && fce_w_str(w, "\n");
        for (i = 0; i < rows; i++)
        {
            uint32_t j;
            ok = ok && fce_w_printf(w, "%u", (unsigned)i);
            for (j = 0; j < 5u; j++)
                ok = ok && fce_w_printf(w, ",%.17g", sos[5u * i + j]);
            if (r->q15)
                for (j = 0; j < 5u; j++)
                    ok = ok && fce_w_printf(w, ",%d",
                                            (int)r->q15[5u * i + j]);
            if (r->q31)
                for (j = 0; j < 5u; j++)
                    ok = ok && fce_w_printf(w, ",%ld",
                                            (long)r->q31[5u * i + j]);
            ok = ok && fce_w_str(w, "\n");
        }
    }
    (void)cols;
    return ok ? FCE_OK : FCE_ERR_BUFFER_TOO_SMALL;
}

/* ------------------------------------------------------------------ */
/* JSON                                                                */
/* ------------------------------------------------------------------ */

static int fce_json_double_array(fce_writer_t* w, const char* key,
                                 const double* v, uint32_t n)
{
    uint32_t i;
    int ok = fce_w_printf(w, "  \"%s\": [", key);
    for (i = 0; i < n; i++)
        ok = ok && fce_w_printf(w, "%.17g%s", v[i], (i + 1u < n) ? ", " : "");
    ok = ok && fce_w_str(w, "]");
    return ok;
}

/* float variant: values convert on the fly (no big stack buffer) */
static int fce_json_float_array(fce_writer_t* w, const char* key,
                                const float* v, uint32_t n)
{
    uint32_t i;
    int ok = fce_w_printf(w, "  \"%s\": [", key);
    for (i = 0; i < n; i++)
        ok = ok && fce_w_printf(w, "%.9g%s", (double)v[i],
                                (i + 1u < n) ? ", " : "");
    ok = ok && fce_w_str(w, "]");
    return ok;
}

fce_status_t fce_export_json(const fce_result_t* r, fce_writer_t* w)
{
    int ok;
    uint32_t i;
    if (r == NULL || w == NULL || w->write == NULL)
        return FCE_ERR_INVALID_ARGUMENT;
    /* require the primary array; without it the document would be
     * malformed (the metadata block ends with a pending comma) */
    if (r->kind == FCE_KIND_FIR && (r->h_f64 == NULL || r->num_taps == 0u))
        return FCE_ERR_NOT_AVAILABLE;
    if (r->kind == FCE_KIND_IIR && (r->sos_f64 == NULL || r->num_sections == 0u))
        return FCE_ERR_NOT_AVAILABLE;

    ok = fce_w_str(w, "{\n");
    ok = ok && fce_w_printf(w, "  \"library\": \"FilterCoeff\",\n");
    ok = ok && fce_w_printf(w, "  \"version\": \"" FCE_VERSION_STRING "\",\n");
    ok = ok && fce_w_printf(w, "  \"kind\": \"%s\",\n", fce_kind_str(r->kind));
    ok = ok && fce_w_printf(w, "  \"fs\": %.17g,\n", r->fs);
    ok = ok && fce_w_printf(w, "  \"fc1\": %.17g,\n", r->fc1);
    if (r->kind == FCE_KIND_FIR)
    {
        ok = ok && fce_w_printf(w, "  \"type\": \"%s\",\n",
                                fce_fir_type_str(r->fir_type));
        ok = ok && fce_w_printf(w, "  \"num_taps\": %u,\n",
                                (unsigned)r->num_taps);
        ok = ok && fce_w_printf(w, "  \"window\": \"%s\",\n",
                                fce_window_str(r->window));
        ok = ok && fce_w_printf(w, "  \"kaiser_beta\": %.17g,\n",
                                r->kaiser_beta);
    }
    else
    {
        ok = ok && fce_w_printf(w, "  \"type\": \"%s\",\n",
                                fce_iir_type_str(r->iir_type));
        ok = ok && fce_w_printf(w, "  \"family\": \"%s\",\n",
                                fce_iir_family_str(r->iir_family));
        ok = ok && fce_w_printf(w, "  \"order\": %u,\n", (unsigned)r->order);
        ok = ok && fce_w_printf(w, "  \"num_sections\": %u,\n",
                                (unsigned)r->num_sections);
        ok = ok && fce_w_printf(w, "  \"design_fc1\": %.17g,\n",
                                r->design_fc1);
        ok = ok && fce_w_printf(w, "  \"design_fc2\": %.17g,\n",
                                r->design_fc2);
    }
    ok = ok && fce_w_printf(w, "  \"precision\": \"%s\",\n",
                            (r->precision == FCE_PRECISION_FLOAT32)
                                ? "float32" : "float64");
    ok = ok && fce_w_printf(w, "  \"normalization\": \"%s\",\n",
                            fce_norm_str(r->normalization));
    ok = ok && fce_w_printf(w, "  \"norm_factor\": %.17g,\n", r->norm_factor);
    ok = ok && fce_w_printf(w, "  \"qformat\": \"%s\",\n",
                            fce_qformat_str(r->qformat));
    ok = ok && fce_w_printf(w, "  \"flags\": %u,\n", (unsigned)r->flags);

    if (r->kind == FCE_KIND_FIR)
    {
        if (r->h_f64)
            ok = ok && fce_json_double_array(w, "coefficients_f64", r->h_f64,
                                             r->num_taps);
        if (r->h_f32)
        {
            ok = ok && fce_w_str(w, ",\n");
            ok = ok && fce_json_float_array(w, "coefficients_f32", r->h_f32,
                                            r->num_taps);
        }
        if (r->q15)
        {
            ok = ok && fce_w_str(w, ",\n  \"q15\": [");
            for (i = 0; i < r->num_taps; i++)
                ok = ok && fce_w_printf(w, "%d%s", (int)r->q15[i],
                                        (i + 1u < r->num_taps) ? ", " : "");
            ok = ok && fce_w_str(w, "]");
        }
        if (r->q31)
        {
            ok = ok && fce_w_str(w, ",\n  \"q31\": [");
            for (i = 0; i < r->num_taps; i++)
                ok = ok && fce_w_printf(w, "%ld%s", (long)r->q31[i],
                                        (i + 1u < r->num_taps) ? ", " : "");
            ok = ok && fce_w_str(w, "]");
        }
    }
    else
    {
        if (r->sos_f64)
        {
            ok = ok && fce_w_str(w, "  \"sos_f64\": [\n");
            for (i = 0; i < r->num_sections; i++)
            {
                uint32_t j;
                ok = ok && fce_w_str(w, "    [");
                for (j = 0; j < 5u; j++)
                    ok = ok && fce_w_printf(w, "%.17g%s",
                                            r->sos_f64[5u * i + j],
                                            (j + 1u < 5u) ? ", " : "");
                ok = ok && fce_w_printf(w, "]%s\n",
                                        (i + 1u < r->num_sections) ? "," : "");
            }
            ok = ok && fce_w_str(w, "  ]");
        }
        if (r->q15)
        {
            ok = ok && fce_w_str(w, ",\n  \"sos_q15\": [\n");
            for (i = 0; i < r->num_sections; i++)
            {
                uint32_t j;
                ok = ok && fce_w_str(w, "    [");
                for (j = 0; j < 5u; j++)
                    ok = ok && fce_w_printf(w, "%d%s",
                                            (int)r->q15[5u * i + j],
                                            (j + 1u < 5u) ? ", " : "");
                ok = ok && fce_w_printf(w, "]%s\n",
                                        (i + 1u < r->num_sections) ? "," : "");
            }
            ok = ok && fce_w_str(w, "  ]");
        }
        if (r->q31)
        {
            ok = ok && fce_w_str(w, ",\n  \"sos_q31\": [\n");
            for (i = 0; i < r->num_sections; i++)
            {
                uint32_t j;
                ok = ok && fce_w_str(w, "    [");
                for (j = 0; j < 5u; j++)
                    ok = ok && fce_w_printf(w, "%ld%s",
                                            (long)r->q31[5u * i + j],
                                            (j + 1u < 5u) ? ", " : "");
                ok = ok && fce_w_printf(w, "]%s\n",
                                        (i + 1u < r->num_sections) ? "," : "");
            }
            ok = ok && fce_w_str(w, "  ]");
        }
        ok = ok && fce_w_printf(w, ",\n  \"scale\": %.17g", r->scale);
        if (r->section_scales)
        {
            /* leading comma of the previous key lives in the array
             * helper's caller: without it (or with a stray trailing
             * comma on "scale") the document is malformed JSON */
            ok = ok && fce_w_str(w, ",\n");
            ok = ok && fce_json_double_array(w, "section_scales",
                                             r->section_scales,
                                             r->num_sections);
        }
    }

    ok = ok && fce_w_str(w, "\n}\n");
    return ok ? FCE_OK : FCE_ERR_BUFFER_TOO_SMALL;
}

/* ------------------------------------------------------------------ */
/* Markdown report                                                     */
/* ------------------------------------------------------------------ */

fce_status_t fce_export_report(const fce_result_t* r, fce_writer_t* w)
{
    int ok;
    const char* status;
    if (r == NULL || w == NULL || w->write == NULL)
        return FCE_ERR_INVALID_ARGUMENT;

    ok = fce_w_str(w, "# FILTER COEFFICIENT REPORT\n\n");
    ok = ok && fce_w_printf(w, "**Library:** FilterCoeff v" FCE_VERSION_STRING "\n\n");
    ok = ok && fce_w_printf(w, "**Status:** %s\n\n", fce_status_str(r->status));
    ok = ok && fce_w_printf(w, "**Flags:** 0x%04X\n\n", (unsigned)r->flags);

    ok = ok && fce_w_str(w, "## Design\n\n");
    ok = ok && fce_w_printf(w, "| Parameter | Value |\n|---|---|\n");
    ok = ok && fce_w_printf(w, "| Type | %s %s |\n", fce_kind_str(r->kind),
                            (r->kind == FCE_KIND_FIR)
                                ? fce_fir_type_str(r->fir_type)
                                : fce_iir_type_str(r->iir_type));
    if (r->kind == FCE_KIND_IIR)
        ok = ok && fce_w_printf(w, "| Family | %s |\n",
                                fce_iir_family_str(r->iir_family));
    ok = ok && fce_w_printf(w, "| Sample rate | %.4f Hz |\n", r->fs);
    ok = ok && fce_w_printf(w, "| Cutoff | %.4f Hz |\n", r->fc1);
    if (r->kind == FCE_KIND_IIR && r->num_sections > 0u)
        ok = ok && fce_w_printf(w, "| Order / Sections | %u / %u |\n",
                                (unsigned)r->order,
                                (unsigned)r->num_sections);
    if (r->kind == FCE_KIND_FIR)
    {
        ok = ok && fce_w_printf(w, "| Taps | %u |\n", (unsigned)r->num_taps);
        ok = ok && fce_w_printf(w, "| Window | %s |\n",
                                fce_window_str(r->window));
        ok = ok && fce_w_printf(w, "| Kaiser beta | %.4f |\n", r->kaiser_beta);
        ok = ok && fce_w_printf(w, "| Symmetry | Type %d |\n",
                                (int)r->symmetry);
    }
    ok = ok && fce_w_printf(w, "| Precision | %s |\n",
                            (r->precision == FCE_PRECISION_FLOAT32)
                                ? "float32 (internal float64)"
                                : "float64");
    ok = ok && fce_w_printf(w, "| Normalization | %s |\n",
                            fce_norm_str(r->normalization));
    ok = ok && fce_w_printf(w, "| Norm factor | %.9g |\n", r->norm_factor);

    if (r->qformat != FCE_QFORMAT_NONE)
    {
        ok = ok && fce_w_str(w, "\n## Fixed Point\n\n");
        ok = ok && fce_w_printf(w, "| Parameter | Value |\n|---|---|\n");
        ok = ok && fce_w_printf(w, "| Format | %s |\n",
                                fce_qformat_str(r->qformat));
        ok = ok && fce_w_printf(w, "| Scale | %.9g |\n", r->scale);
        ok = ok && fce_w_printf(w, "| Integer bits | %u |\n",
                                (unsigned)r->q_int_bits);
        ok = ok && fce_w_printf(w, "| Max abs error | %.6g |\n",
                                r->q_max_abs_error);
        ok = ok && fce_w_printf(w, "| RMS error | %.6g |\n", r->q_rms_error);
        ok = ok && fce_w_printf(w, "| Max rel error | %.6g |\n",
                                r->q_max_rel_error);
        if (r->quant_response_max_error_db > 0.0)
            ok = ok && fce_w_printf(w, "| Max response error | %.6g dB |\n",
                                    r->quant_response_max_error_db);
    }

    ok = ok && fce_w_str(w, "\n## Validation\n\n");
    ok = ok && fce_w_printf(w, "| Check | Value |\n|---|---|\n");
    if (r->kind == FCE_KIND_IIR)
    {
        ok = ok && fce_w_printf(w, "| Max pole radius | %.9g |\n",
                                r->max_pole_radius);
        ok = ok && fce_w_printf(w, "| Stability margin | %.9g |\n",
                                r->stability_margin);
    }
    ok = ok && fce_w_db_row(w, "DC gain", r->dc_gain_db);
    ok = ok && fce_w_db_row(w, "Nyquist gain", r->nyquist_gain_db);
    ok = ok && fce_w_printf(w, "| Passband ripple | %.6f dB |\n",
                            r->passband_ripple_measured_db);
    ok = ok && fce_w_printf(w, "| Stopband attenuation | %.4f dB |\n",
                            r->stopband_atten_measured_db);
    if (r->cutoff_measured_hz > 0.0)
        ok = ok && fce_w_printf(w, "| Measured cutoff | %.4f Hz |\n",
                                r->cutoff_measured_hz);

    status = (r->status == FCE_OK) ? "PASS" : "FAIL";
    ok = ok && fce_w_printf(w, "\n**Validation: %s**\n", status);
    return ok ? FCE_OK : FCE_ERR_BUFFER_TOO_SMALL;
}

#endif /* FCE_ENABLE_EXPORT */
