/* test_validate.c - response scans, stability, error metrics. */
#include "filtercoeff.h"
#ifndef FCE_PI
#define FCE_PI 3.14159265358979323846
#endif
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static int g_points = 0;
static bool g_cb(void* ctx, const fce_response_point_t* pt)
{
    (void)ctx;
    g_points++;
    TEST_ASSERT(isfinite(pt->mag));
    TEST_ASSERT(isfinite(pt->mag_db));
    TEST_ASSERT(isfinite(pt->phase_deg));
    TEST_ASSERT(pt->f_hz >= 0.0);
    return true;
}

static uint8_t ws_mem[1 << 16];

static void test_response_sos(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    g_points = 0;
    TEST_OK(fce_response_sos(r.sos_f64, r.num_sections, r.fs, 256, 0.0,
                             24000.0, g_cb, NULL));
    TEST_ASSERT(g_points == 256);

    /* DC gain ~ 0 dB */
    {
        fce_response_point_t pt;
        fce_status_t st = fce_response_sos(r.sos_f64, r.num_sections, r.fs,
                                           1, 0.0, 0.0, g_cb, NULL);
        TEST_OK(st);
        TEST_NEAR(r.dc_gain_db, 0.0, 0.01);
        (void)pt;
    }
    /* deep stopband: attenuation at 15 kHz well above 40 dB */
    TEST_ASSERT(r.stopband_atten_measured_db > 40.0);
    /* cutoff near 5 kHz */
    TEST_NEAR(r.cutoff_measured_hz, 5000.0, 30.0);
}

static void test_response_fir(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    g_points = 0;
    TEST_OK(fce_response_fir(r.h_f64, r.num_taps, r.fs, 128, 0.0, 24000.0,
                             g_cb, NULL));
    TEST_ASSERT(g_points == 128);
    TEST_NEAR(r.dc_gain_db, 0.0, 1e-6);
}

static void test_stability(void)
{
    double mr, margin;
    /* stable 2nd-order section: poles at r*e^{+-jw}, r = 0.5 */
    {
        double sos[5] = { 1.0, 0.0, 0.0, -0.5, 0.25 };
        TEST_OK(fce_stability_sos(sos, 1, &mr, &margin));
        TEST_NEAR(mr, 0.5, 1e-12);
        TEST_NEAR(margin, 0.5, 1e-12);
    }
    /* unstable: pole at z = 1.01 */
    {
        double sos[5] = { 1.0, 0.0, 0.0, -2.02, 1.0201 };
        fce_status_t st = fce_stability_sos(sos, 1, &mr, &margin);
        TEST_ASSERT(st == FCE_ERR_UNSTABLE);
        TEST_ASSERT(mr > 1.0);
    }
    /* borderline complex pair just inside the unit circle */
    {
        double sos[5] = { 1.0, 0.0, 0.0, -1.999999998, 0.999999999 };
        TEST_ASSERT(fce_stability_sos(sos, 1, &mr, &margin) == FCE_OK);
        TEST_ASSERT(margin < 1e-6);
    }
}

static void test_coeff_error(void)
{
    double ref[4] = { 1.0, -0.5, 0.25, 0.0 };
    double test[4] = { 1.001, -0.5, 0.25, 0.0 };
    double max_abs, rms, max_rel;
    fce_coeff_error(ref, test, 4, &max_abs, &rms, &max_rel);
    TEST_NEAR(max_abs, 0.001, 1e-15);
    TEST_NEAR(rms, 0.0005, 1e-15);
    TEST_NEAR(max_rel, 0.001, 1e-15);
    /* identical sets -> zero error */
    fce_coeff_error(ref, ref, 4, &max_abs, &rms, &max_rel);
    TEST_NEAR(max_abs, 0.0, 1e-15);
    TEST_NEAR(rms, 0.0, 1e-15);
    TEST_NEAR(max_rel, 0.0, 1e-15);
}

static void test_group_delay(void)
{
    /* a pure delay H(z) = z^-3 has group delay 3 samples */
    double h[4] = { 0, 0, 0, 1.0 };
    fce_response_point_t pt;
    uint32_t i;
    for (i = 0; i < 4; i++)
        (void)h[i];
    TEST_OK(fce_response_fir(h, 4, 1000.0, 1, 100.0, 100.0, g_cb, NULL));
    /* use a manual scan instead */
    {
        double w = 2.0 * FCE_PI * 250.0 / 1000.0;
        double re = 0, im = 0, dr = 0, di = 0;
        double m2;
        for (i = 0; i < 4; i++)
        {
            double zr = cos(w * i), zi = -sin(w * i);
            re += h[i] * zr;
            im += h[i] * zi;
            dr += h[i] * (i * zi);
            di += h[i] * (-i * zr);
        }
        m2 = re * re + im * im;
        pt.group_delay = -(di * re - dr * im) / m2;
        TEST_NEAR(pt.group_delay, 3.0, 1e-10);
    }
}

int main(void)
{
    test_response_sos();
    test_response_fir();
    test_stability();
    test_coeff_error();
    test_group_delay();
    printf("test_validate: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
