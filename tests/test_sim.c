/* test_sim.c - simulation module (validation aid) sanity checks. */
#include "filtercoeff.h"
#ifndef FCE_PI
#define FCE_PI 3.14159265358979323846
#endif
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static void test_signals(void)
{
    double buf[256];
    uint32_t n;

    n = fce_sim_signal(FCE_SIM_SINE, 1000.0, 1.0, 8000.0, buf, 256);
    TEST_ASSERT(n == 256);
    TEST_NEAR(buf[0], 0.0, 1e-12);
    TEST_NEAR(buf[2], 1.0, 1e-6);  /* 1000 Hz at 8 kHz: pi/2 at sample 2 */
    TEST_NEAR(buf[4], 0.0, 1e-6);  /* pi at sample 4 */

    n = fce_sim_signal(FCE_SIM_IMPULSE, 10.0, 2.0, 1000.0, buf, 256);
    TEST_ASSERT(n == 256);
    TEST_NEAR(buf[10], 2.0, 1e-15);
    TEST_NEAR(buf[9], 0.0, 1e-15);

    n = fce_sim_signal(FCE_SIM_STEP, 5.0, 1.0, 1000.0, buf, 256);
    TEST_ASSERT(n == 256);
    TEST_NEAR(buf[4], 0.0, 1e-15);
    TEST_NEAR(buf[5], 1.0, 1e-15);

    n = fce_sim_signal(FCE_SIM_WHITE_NOISE, 1.0, 1.0, 1000.0, buf, 256);
    TEST_ASSERT(n == 256);
    {
        double mx = 0.0;
        uint32_t i;
        for (i = 0; i < 256; i++)
            if (fabs(buf[i]) > mx)
                mx = fabs(buf[i]);
        TEST_ASSERT(mx <= 1.0);
        TEST_ASSERT(mx > 0.9); /* statistically certain for 256 samples */
    }
}

static void test_fir_impulse_response(void)
{
    /* impulse response of a filter must equal its coefficients */
    fce_spec_t sp;
    fce_result_t r;
    uint8_t mem[1 << 16];
    fce_workspace_t ws = { mem, sizeof(mem) };
    double x[128] = { 0 };
    double y[128];
    uint32_t i;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 64;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    x[0] = 1.0;
    fce_sim_fir(r.h_f64, r.num_taps, x, y, 128);
    for (i = 0; i < r.num_taps; i++)
        TEST_NEAR(y[i], r.h_f64[i], 1e-15);
    for (i = r.num_taps; i < 128; i++)
        TEST_NEAR(y[i], 0.0, 1e-15);
}

static void test_sine_amplitude_matches_response(void)
{
    /* sine at 1 kHz through a 5 kHz LP: amplitude = |H(1k)| */
    fce_spec_t sp;
    fce_result_t r;
    uint8_t mem[1 << 16];
    fce_workspace_t ws = { mem, sizeof(mem) };
    double x[4096];
    double y[4096];
    double state[8] = { 0 };
    double amp = 0.0;
    uint32_t i;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));

    fce_sim_signal(FCE_SIM_SINE, 1000.0, 1.0, 48000.0, x, 4096);
    fce_sim_sos(r.sos_f64, r.num_sections, x, y, 4096, state);

    /* steady-state amplitude via RMS over the last 20 periods
     * (20 periods at 1 kHz in 48 kHz = 960 samples) */
    {
        double p = 0.0;
        for (i = 4096 - 960; i < 4096; i++)
            p += y[i] * y[i];
        amp = sqrt(2.0 * p / 960.0);
    }

    /* |H| at 1 kHz from the response scan */
    {
        double w = 2.0 * FCE_PI * 1000.0 / 48000.0;
        double g = 1.0;
        uint32_t s;
        for (s = 0; s < r.num_sections; s++)
        {
            double cw = cos(w), sw = sin(w);
            double c2w = cos(2.0 * w), s2w = sin(2.0 * w);
            double br = r.sos_f64[5 * s] + r.sos_f64[5 * s + 1] * cw +
                        r.sos_f64[5 * s + 2] * c2w;
            double bi = -(r.sos_f64[5 * s + 1] * sw +
                          r.sos_f64[5 * s + 2] * s2w);
            double ar = 1.0 + r.sos_f64[5 * s + 3] * cw +
                        r.sos_f64[5 * s + 4] * c2w;
            double ai = -(r.sos_f64[5 * s + 3] * sw +
                          r.sos_f64[5 * s + 4] * s2w);
            g *= sqrt(br * br + bi * bi) / sqrt(ar * ar + ai * ai);
        }
        TEST_NEAR(amp, g, 1e-3);
    }
}

static void test_multitone_energy(void)
{
    double buf[512];
    /* fs = 16 kHz so the third tone (2*f1+f2 = 4 kHz) is below Nyquist */
    uint32_t n = fce_sim_signal(FCE_SIM_MULTITONE, 1000.0, 2000.0, 16000.0,
                                buf, 512);
    TEST_ASSERT(n == 512);
    /* three tones at 1k, 2k, 4k: power ~ 1.5 * amplitude^2 */
    {
        double p = 0.0;
        uint32_t i;
        for (i = 0; i < 512; i++)
            p += buf[i] * buf[i];
        p /= 512.0;
        TEST_NEAR(p, 1.5, 0.05);
    }
}

int main(void)
{
    test_signals();
    test_fir_impulse_response();
    test_sine_amplitude_matches_response();
    test_multitone_energy();
    printf("test_sim: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
