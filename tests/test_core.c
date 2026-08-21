/* test_core.c - workspace sizing, spec validation, error codes. */
#include "filtercoeff.h"
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static void test_status_str(void)
{
    TEST_ASSERT(strcmp(fce_status_str(FCE_OK), "OK") == 0);
    TEST_ASSERT(fce_status_str((fce_status_t)999) != NULL);
}

static void test_workspace_sizing(void)
{
    fce_spec_t sp;
    size_t need;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    need = fce_workspace_required(&sp);
    TEST_ASSERT(need > 0);
    TEST_ASSERT(need < 8192); /* 101 taps -> ~4.3 kB */

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 8;
    need = fce_workspace_required(&sp);
    TEST_ASSERT(need > 0);
    TEST_ASSERT(need < 16384);

    /* invalid spec -> 0 */
    fce_spec_defaults(&sp);
    sp.kind = (fce_kind_t)77;
    TEST_ASSERT(fce_workspace_required(&sp) == 0);

    TEST_ASSERT(fce_workspace_required(NULL) == 0);
    TEST_ASSERT(fce_workspace_required_max() > 0);
}

static void test_invalid_specs(void)
{
    fce_spec_t sp;
    fce_result_t r;
    uint8_t mem[16384];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_status_t st;

    /* fs = 0 */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fs = 0; sp.fc1 = 100; sp.num_taps = 32;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* cutoff above Nyquist */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fs = 1000; sp.fc1 = 600; sp.num_taps = 32;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* BP with fc2 <= fc1 */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fir_type = FCE_FIR_BANDPASS;
    sp.fs = 1000; sp.fc1 = 300; sp.fc2 = 200; sp.num_taps = 32;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* auto taps without Kaiser */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
    sp.window = FCE_WIN_HAMMING;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* auto taps with Kaiser but no transition */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
    sp.window = FCE_WIN_KAISER; sp.stopband_atten_db = 60;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* cheby1 without ripple */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_CHEBYSHEV1;
    sp.iir_type = FCE_IIR_LOWPASS; sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* ellip without attenuation */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS; sp.fs = 48000; sp.fc1 = 5000;
    sp.order = 4; sp.passband_ripple_db = 0.5;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_INVALID_SPEC);

    /* Bessel auto order not supported */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BESSEL;
    sp.iir_type = FCE_IIR_LOWPASS; sp.fs = 48000; sp.fc1 = 5000; sp.order = 0;
    sp.edge1_hz = 8000; sp.stopband_atten_db = 40;
    st = fce_generate(&sp, &r, &ws);
    TEST_ASSERT(st == FCE_ERR_UNSUPPORTED);

    /* too-small workspace */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    {
        uint8_t small[64];
        fce_workspace_t ws2 = { small, sizeof(small) };
        st = fce_generate(&sp, &r, &ws2);
        TEST_ASSERT(st == FCE_ERR_BUFFER_TOO_SMALL);
    }

    /* null args */
    TEST_ASSERT(fce_generate(NULL, &r, &ws) == FCE_ERR_INVALID_ARGUMENT);
    TEST_ASSERT(fce_generate(&sp, NULL, &ws) == FCE_ERR_INVALID_ARGUMENT);
}

int main(void)
{
    test_status_str();
    test_workspace_sizing();
    test_invalid_specs();
    printf("test_core: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
