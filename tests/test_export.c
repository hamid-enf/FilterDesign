/* test_export.c - C/CSV/JSON/report exporters. */
#include "filtercoeff.h"
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static uint8_t ws_mem[1 << 16];

static void test_export_c_fir(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;
    fce_export_opts_t o;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 33;
    sp.window = FCE_WIN_HAMMING;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    memset(&o, 0, sizeof(o));
    o.name = "coeffs";
    o.include_quantized = 1;
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_c_fir(&r, &o, &w));
    TEST_ASSERT(mw.truncated == 0);
    TEST_ASSERT(strstr(buf, "static const float64_t coeffs[33]") != NULL);
    TEST_ASSERT(strstr(buf, "static const int16_t coeffs_q15[33]") != NULL);
    TEST_ASSERT(strstr(buf, "FilterCoeff") != NULL);

    /* small buffer -> FCE_ERR_BUFFER_TOO_SMALL */
    fce_writer_mem_init(&w, &mw, buf, 64);
    TEST_ASSERT(fce_export_c_fir(&r, &o, &w) == FCE_ERR_BUFFER_TOO_SMALL);
}

static void test_export_c_sos(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_c_sos(&r, NULL, &w));
    TEST_ASSERT(strstr(buf, "static const float64_t sos[][5]") != NULL);
    TEST_ASSERT(strstr(buf, "a1*y[n-1]") != NULL);
    /* count rows: two sections */
    {
        int rows = 0;
        char* p = buf;
        while ((p = strstr(p, "    { ")) != NULL)
        {
            rows++;
            p += 6;
        }
        TEST_ASSERT(rows == 2);
    }
}

static void test_export_csv_json_report(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    char buf[32768];
    fce_mem_writer_t mw;
    fce_writer_t w;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    /* CSV */
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_csv(&r, &w));
    TEST_ASSERT(strstr(buf, "section,b0,b1,b2,a1,a2") != NULL);
    TEST_ASSERT(strstr(buf, "0,") != NULL);

    /* JSON */
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_json(&r, &w));
    TEST_ASSERT(strstr(buf, "\"library\": \"FilterCoeff\"") != NULL);
    TEST_ASSERT(strstr(buf, "\"sos_f64\"") != NULL);
    TEST_ASSERT(strstr(buf, "\"family\": \"Elliptic\"") != NULL);

    /* report */
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_report(&r, &w));
    TEST_ASSERT(strstr(buf, "FILTER COEFFICIENT REPORT") != NULL);
    TEST_ASSERT(strstr(buf, "Validation: PASS") != NULL);
    TEST_ASSERT(strstr(buf, "Elliptic") != NULL);
}

/* Cheap structural JSON check: no dangling comma before a closing
 * brace/bracket, balanced braces/brackets (values here never embed
 * braces inside strings).  Regression: the IIR export used to emit
 * "scale": 0,\n} for unquantized designs -- invalid JSON. */
static int json_sane(const char* s)
{
    int braces = 0, brackets = 0;
    const char* p;
    for (p = s; *p; p++)
    {
        if (*p == ',')
        {
            const char* q = p + 1;
            while (*q == ' ' || *q == '\n' || *q == '\t')
                q++;
            if (*q == '}' || *q == ']')
                return 0;
        }
        else if (*p == '{') braces++;
        else if (*p == '}') braces--;
        else if (*p == '[') brackets++;
        else if (*p == ']') brackets--;
    }
    return braces == 0 && brackets == 0;
}

static void test_export_json_wellformed(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    char buf[32768];
    fce_mem_writer_t mw;
    fce_writer_t w;

    /* IIR, no quantization: no section_scales array is emitted. */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_json(&r, &w));
    TEST_ASSERT(strstr(buf, "\"scale\": 0") != NULL);
    TEST_ASSERT(strstr(buf, "section_scales") == NULL);
    TEST_ASSERT(json_sane(buf));

    /* IIR with quantization: section_scales present, still well formed. */
    sp.qformat = FCE_QFORMAT_Q31;
    sp.scale_strategy = FCE_SCALE_SECTION_WISE;
    TEST_OK(fce_generate(&sp, &r, &ws));
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_json(&r, &w));
    TEST_ASSERT(strstr(buf, "section_scales") != NULL);
    TEST_ASSERT(json_sane(buf));

    /* FIR float32: well formed too. */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 33;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT32;
    TEST_OK(fce_generate(&sp, &r, &ws));
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_json(&r, &w));
    TEST_ASSERT(json_sane(buf));
}

static void test_export_float32(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 33;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT32;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.h_f32 != NULL);
    TEST_ASSERT(r.h_f64 != NULL); /* doubles always available */

    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    TEST_OK(fce_export_c_fir(&r, NULL, &w));
    TEST_ASSERT(strstr(buf, "static const float32_t coeffs[33]") != NULL);
}

int main(void)
{
    test_export_c_fir();
    test_export_c_sos();
    test_export_csv_json_report();
    test_export_json_wellformed();
    test_export_float32();
    printf("test_export: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
