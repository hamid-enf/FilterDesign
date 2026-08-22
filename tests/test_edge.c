/* test_edge.c - boundary and edge cases. */
#include "filtercoeff.h"
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static uint8_t ws_mem[1 << 16];

static fce_status_t design(fce_spec_t* sp, fce_result_t* r)
{
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    return fce_generate(sp, r, &ws);
}

static void test_cutoff_near_zero(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 1.0; sp.num_taps = 101;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    /* narrow lowpass: DC gain still 1 */
    TEST_NEAR(r.dc_gain_db, 0.0, 1e-6);
}

static void test_cutoff_near_nyquist(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 23999.0; sp.num_taps = 101;
    sp.window = FCE_WIN_HAMMING;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_NEAR(r.dc_gain_db, 0.0, 1e-6);

    /* IIR with cutoff very close to Nyquist */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 23990.0; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.max_pole_radius < 1.0);
}

static void test_narrow_transition(void)
{
    /* Kaiser auto-taps for a very narrow transition: must clamp and flag.
     * The workspace must be sized for the clamped tap count. */
    uint8_t big[1 << 18];
    fce_workspace_t ws_big = { big, sizeof(big) };
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 1000; sp.num_taps = 0;
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 120; sp.transition_hz = 5.0;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws_big));
    TEST_ASSERT(r.num_taps == FCE_MAX_FIR_TAPS);
    TEST_ASSERT((r.flags & FCE_FLAG_ORDER_CLAMPED) != 0);
}

static void test_high_attenuation(void)
{
    fce_spec_t sp;
    fce_result_t r;
    /* elliptic with 120 dB stopband */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 8;
    sp.passband_ripple_db = 0.1; sp.stopband_atten_db = 120;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.max_pole_radius < 1.0);
    TEST_ASSERT(r.stopband_atten_measured_db > 100.0);
}

static void test_iir_near_instability(void)
{
    /* very narrow bandpass: poles close to the unit circle */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 9990; sp.fc2 = 10010; sp.order = 2;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.max_pole_radius < 1.0);
    /* very close to 1 */
    TEST_ASSERT(r.max_pole_radius > 0.99);
}

static void test_tiny_coefficients(void)
{
    /* very low fc FIR: coefficients span many orders of magnitude */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 10.0; sp.num_taps = 201;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 8.0;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.h_f64[100] > 1e-4);  /* center tap ~ 2fc/fs */
    TEST_ASSERT(r.h_f64[0] < 0.01);    /* edge taps decay slowly */
    TEST_ASSERT(isfinite(r.h_f64[0]));
}

static void test_float32_precision(void)
{
    /* float32 output: coefficients are the float-rounded doubles */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.precision = FCE_PRECISION_FLOAT32;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.h_f32 != NULL);
    {
        uint32_t i;
        for (i = 0; i < 101; i++)
        {
            /* float32 rounding error < 2^-23 relative */
            double d = fabs((double)r.h_f32[i] - r.h_f64[i]);
            TEST_ASSERT(d < 6e-8 * (1.0 + fabs(r.h_f64[i])));
        }
    }
}

static void test_odd_even_orders(void)
{
    uint16_t order;
    for (order = 1; order <= 12; order++)
    {
        fce_spec_t sp;
        fce_result_t r;
        fce_spec_defaults(&sp);
        sp.kind = FCE_KIND_IIR;
        sp.iir_family = FCE_IIR_ELLIPTIC;
        sp.iir_type = FCE_IIR_HIGHPASS;
        sp.fs = 48000; sp.fc1 = 3000; sp.order = order;
        sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
        sp.precision = FCE_PRECISION_FLOAT64;
        TEST_OK(design(&sp, &r));
        TEST_ASSERT(r.max_pole_radius < 1.0);
        TEST_ASSERT(r.num_sections == (uint16_t)((order + 1) / 2));
    }
}

static void test_null_and_abuse(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 33;
    /* null workspace */
    TEST_ASSERT(fce_generate(&sp, &r, NULL) == FCE_ERR_INVALID_ARGUMENT);
    /* workspace without data */
    {
        fce_workspace_t ws2 = { NULL, 100 };
        TEST_ASSERT(fce_generate(&sp, &r, &ws2) == FCE_ERR_INVALID_ARGUMENT);
    }
}

int main(void)
{
    test_cutoff_near_zero();
    test_cutoff_near_nyquist();
    test_narrow_transition();
    test_high_attenuation();
    test_iir_near_instability();
    test_tiny_coefficients();
    test_float32_precision();
    test_odd_even_orders();
    test_null_and_abuse();
    printf("test_edge: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
