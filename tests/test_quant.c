/* test_quant.c - Q15/Q31 conversion, scaling strategies, error metrics. */
#include "filtercoeff.h"
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static uint8_t ws_mem[1 << 16];

static void test_q15_basic(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    uint32_t i;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.q15 != NULL);
    /* all values within [-32767, 32767] */
    for (i = 0; i < 101; i++)
        TEST_ASSERT(r.q15[i] >= -32767 && r.q15[i] <= 32767);
    /* q = round(c * scale) */
    for (i = 0; i < 101; i++)
    {
        double expect = floor(r.h_f64[i] * r.scale + 0.5);
        TEST_ASSERT(abs((int)(expect - r.q15[i])) <= 1);
    }
    /* quantization error < 1/scale */
    TEST_ASSERT(r.q_max_abs_error < 1.001 / r.scale);
    TEST_ASSERT(r.q_max_abs_error > 0.0);
    /* reconstruction error in the response is small (measured above
     * -60 dB, where the response matters) */
    TEST_ASSERT(r.quant_response_max_error_db < 0.2);
    /* integer bits = 0 (all |c| <= 1) */
    TEST_ASSERT(r.q_int_bits == 0);
}

static void test_q31_fir(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    uint32_t i;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.qformat = FCE_QFORMAT_Q31;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    for (i = 0; i < 101; i++)
    {
        TEST_ASSERT(r.q31[i] >= -2147483647L && r.q31[i] <= 2147483647L);
        /* Q31 error much smaller than Q15 */
    }
    TEST_ASSERT(r.q_max_abs_error < 1e-8);
    TEST_ASSERT(r.quant_response_max_error_db < 1e-3);
}

static void test_q15_sos_symmetric(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.scale_strategy = FCE_SCALE_SYMMETRIC;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.q15 != NULL);
    TEST_ASSERT(r.scale > 0.0);
    /* quantized filter remains stable */
    TEST_ASSERT(r.quant_max_pole_radius < 1.0);
    TEST_ASSERT((r.flags & FCE_FLAG_QUANTIZATION_UNSTABLE) == 0);
    /* Q15 with |a1| up to ~1.3 needs 1 integer bit */
    TEST_ASSERT(r.q_int_bits == 1);
}

static void test_q15_sos_sectionwise(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.scale_strategy = FCE_SCALE_SECTION_WISE;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.section_scales != NULL);
    TEST_ASSERT(r.num_sections == 2);
    /* each section's coefficients fit Q15 with its own scale */
    {
        uint32_t s;
        for (s = 0; s < r.num_sections; s++)
        {
            uint32_t j;
            for (j = 0; j < 5; j++)
            {
                int32_t q = r.q15[5 * s + j];
                TEST_ASSERT(q >= -32767 && q <= 32767);
            }
        }
    }
    TEST_ASSERT(r.quant_stability_margin > 0.0);
    /* response error stays small */
    TEST_ASSERT(r.quant_response_max_error_db < 0.1);
}

static void test_q15_coefficientwise(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 51;
    sp.window = FCE_WIN_HAMMING;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.scale_strategy = FCE_SCALE_COEFFICIENT_WISE;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    TEST_ASSERT(r.coeff_scales != NULL);
    /* every nonzero coefficient saturates the format */
    {
        uint32_t i;
        for (i = 0; i < 51; i++)
        {
            if (r.h_f64[i] != 0.0)
                TEST_ASSERT(abs(r.q15[i]) == 32767);
        }
    }
}

static void test_quant_clamp(void)
{
    /* a coefficient > 1 in Q15 must not overflow: scale must shrink */
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_CHEBYSHEV1;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 2;
    sp.passband_ripple_db = 0.1;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(fce_generate(&sp, &r, &ws));
    /* |a1| for a narrow-ish cheby1 can exceed 1 */
    {
        uint32_t i;
        int32_t qmin = 32767, qmax = -32767;
        for (i = 0; i < 5u * r.num_sections; i++)
        {
            if (r.q15[i] < qmin) qmin = r.q15[i];
            if (r.q15[i] > qmax) qmax = r.q15[i];
        }
        TEST_ASSERT(qmin >= -32767 && qmax <= 32767);
    }
}

static void test_q15_unstable_detection(void)
{
    /* an extreme narrow-band filter pushed into Q15 may become unstable;
     * the library must flag it rather than silently deliver */
    fce_spec_t sp;
    fce_result_t r;
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    uint32_t s;
    int saw_flag = 0;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 9900; sp.fc2 = 10100; sp.order = 2;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    (void)fce_generate(&sp, &r, &ws);

    /* either the design is stable in Q15, or the flag/status is set */
    if ((r.flags & FCE_FLAG_QUANTIZATION_UNSTABLE) != 0)
        saw_flag = 1;
    if (r.status == FCE_ERR_QUANTIZATION)
        saw_flag = 1;
    if (!saw_flag)
    {
        /* stable case: verify the margin is reported */
        TEST_ASSERT(r.quant_stability_margin > -1e-9);
    }
    (void)s;
}

int main(void)
{
    test_q15_basic();
    test_q31_fir();
    test_q15_sos_symmetric();
    test_q15_sos_sectionwise();
    test_q15_coefficientwise();
    test_quant_clamp();
    test_q15_unstable_detection();
    printf("test_quant: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
