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

static void test_antisymmetric_parity(void)
{
    /* anti-symmetric responses: odd taps -> Type III, even -> Type IV
     * (it used to be hardcoded III for Hilbert / IV for differentiator) */
    fce_spec_t sp;
    fce_result_t r;

    fce_spec_fir(&sp, FCE_FIR_HILBERT, 48000, 0, 0, 100,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT64);
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_IV);

    fce_spec_fir(&sp, FCE_FIR_DIFFERENTIATOR, 48000, 0, 0, 51,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT64);
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_III);
}

/* |H(f)| of the designed FIR (unquantized float64 coefficients). */
static double fir_gain_at(const fce_result_t* r, double f)
{
    double re = 0.0, im = 0.0;
    uint32_t i;
    for (i = 0; i < r->num_taps; i++)
    {
        double ph = -2.0 * FCE_PI * f * (double)i / 48000.0;
        re += r->h_f64[i] * cos(ph);
        im += r->h_f64[i] * sin(ph);
    }
    return sqrt(re * re + im * im);
}

static void test_even_tap_hp_bs(void)
{
    fce_spec_t sp;
    fce_result_t r;
    uint32_t i;
    double f, peak = 0.0;

    /* Regression: an even-tap highpass must be a true highpass.  The
     * band-limited delta (sinc(m), equal to 1 at m=0) used to be added
     * only when m == 0.0 exactly, which never happens for even lengths
     * (half-integer m) -- the result was a negated LOWPASS with
     * H(0) ~= 1 instead of ~= 0. */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_HIGHPASS;
    sp.fs = 48000; sp.fc1 = 12000; sp.num_taps = 64;
    sp.window = FCE_WIN_HAMMING;
    sp.normalization = FCE_NORM_NONE;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_taps == 64u);
    TEST_ASSERT(r.symmetry == FCE_SYMMETRY_II);
    TEST_ASSERT(fir_gain_at(&r, 0.0) < 1e-3);          /* blocks DC        */
    TEST_NEAR(fir_gain_at(&r, 12000.0), 0.5, 0.05);    /* -6 dB at cutoff  */
    TEST_NEAR(fir_gain_at(&r, 0.40 * 48000.0), 1.0, 0.05); /* passband     */
    TEST_ASSERT(r.flags & FCE_FLAG_SYMMETRY_WARNING);  /* Type II null     */
    for (i = 0; i < 64; i++)
        TEST_NEAR(r.h_f64[i], r.h_f64[63u - i], 1e-14);

    /* AUTO normalization must not divide by the Type-II Nyquist null;
     * it falls back to the passband peak and warns. */
    sp.normalization = FCE_NORM_AUTO;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.flags & FCE_FLAG_SYMMETRY_WARNING);
    TEST_ASSERT(r.normalization == FCE_NORM_PASSBAND_PEAK);
    for (f = 0.26 * 48000.0; f < 0.495 * 48000.0; f += 200.0)
    {
        double g = fir_gain_at(&r, f);
        if (g > peak)
            peak = g;
    }
    TEST_NEAR(peak, 1.0, 0.02); /* passband peak normalized to unity */

    /* Regression: an even-tap bandstop must pass DC (~1), not block it. */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_BANDSTOP;
    sp.fs = 48000; sp.fc1 = 6000; sp.fc2 = 18000; sp.num_taps = 64;
    sp.window = FCE_WIN_HAMMING;
    sp.normalization = FCE_NORM_NONE;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_NEAR(fir_gain_at(&r, 0.0), 1.0, 0.02);        /* passes DC        */
    TEST_ASSERT(fir_gain_at(&r, 12000.0) < 0.01);      /* stopband dip     */
    TEST_NEAR(fir_gain_at(&r, 0.45 * 48000.0), 1.0, 0.05); /* upper pass   */
}

static void test_kaiser_auto_min_taps(void)
{
    /* a very relaxed attenuation spec used to collapse to N = 1 tap,
     * which made the window math divide by (N-1) = 0 and fail with a
     * numerical error; the estimate must clamp to a valid odd N >= 3 */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0; /* auto */
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 1.0; sp.transition_hz = 8000;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_taps >= 3u);
    TEST_ASSERT((r.num_taps & 1u) == 1u);
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
    test_antisymmetric_parity();
    test_even_tap_hp_bs();
    test_kaiser_auto_min_taps();
    printf("test_fir: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
