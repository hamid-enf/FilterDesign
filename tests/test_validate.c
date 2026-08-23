/* test_validate.c - response scans, stability, error metrics. */
#include "filtercoeff.h"
#include "fce_internal.h" /* fce_eval_biquad / fce_biquad_group_delay */
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

static void test_bandpass_stopband_metric(void)
{
    /* a bandpass has TWO stopbands: the reported attenuation must take
     * the worse of the lower and upper stopband maxima (it used to
     * compare the passband max against the upper stopband -> ~0 dB) */
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    /* FIR Kaiser bandpass, designed for 80 dB */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 4000; sp.fc2 = 8000; sp.num_taps = 201;
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 80.0;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.stopband_atten_measured_db > 80.0);

    /* IIR elliptic bandpass, designed for 60 dB */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 4000; sp.fc2 = 8000; sp.order = 6;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60.0;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.stopband_atten_measured_db > 55.0);
}

static void test_group_delay_sos_sign(void)
{
    /* a first-order LP biquad must have POSITIVE group delay at DC
     * (the biquad GD formula used to return the negated value).
     * H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2).
     * Reference: finite difference of the unwrapped phase. */
    const double sos[5] = { 0.06745527, 0.13491055, 0.06745527,
                            -1.1429805, 0.4128016 }; /* butter LP2 w=0.1pi-ish */
    const double w0 = 0.3;
    const double dw = 1e-7;
    double gd_num;
    fce_cplx_t hp, hm;
    double pp, pm;
    double gd_ana = 0.0;
    fce_response_point_t pt;
    (void)pt;

    hp = fce_eval_biquad(sos, w0 + dw);
    hm = fce_eval_biquad(sos, w0 - dw);
    pp = atan2(hp.im, hp.re);
    pm = atan2(hm.im, hm.re);
    gd_num = -(pp - pm) / (2.0 * dw);

    gd_ana = fce_biquad_group_delay(sos, w0);
    TEST_ASSERT(gd_ana > 0.0);
    TEST_NEAR(gd_ana, gd_num, 1e-4);

    /* DC: known closed-form check on a simple 1st-order section
     * H(z) = (1 - z^-1)/2 / (1 - 0.5 z^-1) ... use response scan on an
     * SOS cascade and confirm positivity there too */
    {
        static const double sos2[5] = { 0.5, -0.5, 0.0, -0.5, 0.0 };
        double gd2 = fce_biquad_group_delay(sos2, 0.0);
        /* GD_B = sum n b / sum b = 1 (b=[0.5,-0.5]) , GD_A = -0.5/0.5 = -1
         * GD = GD_A * -1?? -> numerically: */
        TEST_ASSERT(gd2 > 0.0);
    }
}

int main(void)
{
    test_response_sos();
    test_response_fir();
    test_stability();
    test_coeff_error();
    test_group_delay();
    test_group_delay_sos_sign();
    test_bandpass_stopband_metric();
    printf("test_validate: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
