/* test_fir.c - FIR design: correctness, symmetry, normalization, windows. */
#include "filtercoeff.h"
#ifndef FCE_PI
#define FCE_PI 3.14159265358979323846
#endif
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static uint8_t ws_mem[1 << 16];

static fce_status_t design(fce_spec_t* sp, fce_result_t* r)
{
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    return fce_generate(sp, r, &ws);
}

static void test_lowpass_basics(void)
{
    fce_spec_t sp;
    fce_result_t r;
    uint32_t i;
    double dc = 0.0;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_taps == 101);
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_I);

    /* symmetry */
    for (i = 0; i < 101; i++)
        TEST_NEAR(r.h_f64[i], r.h_f64[100u - i], 1e-14);

    /* DC gain = 1 (normalized) */
    for (i = 0; i < 101; i++)
        dc += r.h_f64[i];
    TEST_NEAR(dc, 1.0, 1e-12);

    /* coefficients finite */
    for (i = 0; i < 101; i++)
        TEST_ASSERT(isfinite(r.h_f64[i]));

    /* ideal response available */
    TEST_ASSERT(r.fir_ideal != NULL);
    TEST_ASSERT(r.fir_window != NULL);

    /* Hamming window in result matches the window function */
    TEST_NEAR(r.fir_window[0], 0.08, 1e-14);
    TEST_NEAR(r.fir_window[50], 1.0, 1e-14);
}

static void test_highpass_nyquist(void)
{
    fce_spec_t sp;
    fce_result_t r;
    uint32_t i;
    double g = 0.0;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_HIGHPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    for (i = 0; i < 101; i++)
        g += r.h_f64[i] * ((i & 1u) ? -1.0 : 1.0);
    TEST_NEAR(g, 1.0, 1e-12); /* unity at Nyquist */
    /* DC gain ~ 0 (window leakage only) */
    {
        double dc = 0.0;
        for (i = 0; i < 101; i++)
            dc += r.h_f64[i];
        TEST_NEAR(dc, 0.0, 1e-3);
    }
}

static void test_bandpass_peak(void)
{
    fce_spec_t sp;
    fce_result_t r;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 6000; sp.num_taps = 121;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.norm_factor > 0.9 && r.norm_factor < 1.1);
}

static void test_kaiser_auto(void)
{
    fce_spec_t sp;
    fce_result_t r;

    /* the mission example: Fs=48k, Fc=5k, 80 dB, 1 kHz transition */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 80; sp.transition_hz = 1000;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_taps == 243);
    TEST_NEAR(r.kaiser_beta, 7.8573, 1e-3);
    TEST_ASSERT((r.flags & FCE_FLAG_SPEC_MARGINAL) == 0);

    /* achieved stopband attenuation must meet the spec */
    TEST_ASSERT(r.stopband_atten_measured_db > 79.0);
}

static void test_hilbert_differentiator(void)
{
    fce_spec_t sp;
    fce_result_t r;
    uint32_t i;

    /* Hilbert: odd taps, Type III antisymmetric */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_HILBERT;
    sp.fs = 48000; sp.fc1 = 0; sp.fc2 = 0; sp.num_taps = 65;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_III);
    for (i = 0; i < 65; i++)
        TEST_NEAR(r.h_f64[i], -r.h_f64[64u - i], 1e-12);

    /* differentiator: even taps, Type IV antisymmetric */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_DIFFERENTIATOR;
    sp.fs = 48000; sp.fc1 = 0; sp.fc2 = 0; sp.num_taps = 64;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_IV);
    for (i = 0; i < 64; i++)
        TEST_NEAR(r.h_f64[i], -r.h_f64[63u - i], 1e-12);

    /* differentiator frequency response ~ |H| = w (linear ramp) */
    {
        double f = 12000.0;
        double w = 2.0 * FCE_PI * f / 48000.0;
        double re = 0.0, im = 0.0;
        for (i = 0; i < 64; i++)
        {
            re += r.h_f64[i] * cos(w * i);
            im -= r.h_f64[i] * sin(w * i);
        }
        TEST_NEAR(sqrt(re * re + im * im), 0.5, 0.02); /* f/nyq (windowed) */
    }
}

static void test_symmetry_warning(void)
{
    /* even taps + highpass -> Type II has a Nyquist null: warn */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_HIGHPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 64;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT((r.flags & FCE_FLAG_SYMMETRY_WARNING) != 0);
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_II);
}

static void test_fir_known_coefficients(void)
{
    /* 3-tap windowed-sinc LP at fc = fs/4 (Nyquist/2):
     * ideal h = [0.5*sinc(-0.5), 0.5, 0.5*sinc(0.5)] = [0.31831, 0.5, 0.31831]
     * normalized to unity DC gain */
    {
        double s_ = sin(-FCE_PI * 0.5) / (-FCE_PI * 0.5);
        double h0 = 0.5 * s_;
        double sum = 2.0 * h0 + 0.5;
        TEST_NEAR(h0 / sum, 0.28004957675577868, 1e-12);
    }
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 1000; sp.fc1 = 250; sp.num_taps = 3;
    sp.window = FCE_WIN_RECTANGULAR;
    sp.normalization = FCE_NORM_DC;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_NEAR(r.h_f64[0], 0.28004957675577868, 1e-14);
    TEST_NEAR(r.h_f64[1], 0.43990084648844263, 1e-14);
    TEST_NEAR(r.h_f64[2], 0.28004957675577868, 1e-14);
}

int main(void)
{
    test_lowpass_basics();
    test_highpass_nyquist();
    test_bandpass_peak();
    test_kaiser_auto();
    test_hilbert_differentiator();
    test_symmetry_warning();
    test_fir_known_coefficients();
    printf("test_fir: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
