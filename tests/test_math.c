/* test_math.c - elliptic functions, I0, Kaiser formulas, windows. */
#include "filtercoeff.h"
#include "fce_internal.h"
#include "test_harness.h"

int g_tests_run = 0;
int g_tests_failed = 0;

static void test_ellipk(void)
{
    TEST_NEAR(fce_ellipk(0.0), 1.5707963267948966, 1e-15);
    TEST_NEAR(fce_ellipk(0.5), 1.854074677301372, 1e-14);
    TEST_NEAR(fce_ellipk(0.99), 3.695637362989875, 1e-13);
    /* K(1-m) consistency */
    TEST_NEAR(fce_ellipkm1(0.5), fce_ellipk(0.5), 1e-15);
    TEST_NEAR(fce_ellipkm1(1e-10), 12.899219826387600, 1e-11);
    TEST_NEAR(fce_ellipkm1(0.99), fce_ellipk(0.01), 1e-12);
}

static void test_ellipj(void)
{
    double sn, cn, dn;
    /* m = 0: sn = sin */
    fce_ellipj(0.7, 0.0, &sn, &cn, &dn);
    TEST_NEAR(sn, sin(0.7), 1e-15);
    TEST_NEAR(cn, cos(0.7), 1e-15);
    TEST_NEAR(dn, 1.0, 1e-15);
    /* m = 1: sn = tanh */
    fce_ellipj(0.7, 1.0, &sn, &cn, &dn);
    TEST_NEAR(sn, tanh(0.7), 1e-15);
    /* known values (scipy.special.ellipj) */
    fce_ellipj(0.5, 0.5, &sn, &cn, &dn);
    TEST_NEAR(sn, 0.47075047365565736, 1e-14);
    TEST_NEAR(cn, 0.8822663948904402, 1e-14);
    TEST_NEAR(dn, 0.9429724257773857, 1e-14);
    fce_ellipj(1.0, 0.25, &sn, &cn, &dn);
    TEST_NEAR(sn, 0.8226355781298623, 1e-14);
    /* identity: sn^2 + cn^2 = 1, dn^2 = 1 - m sn^2 */
    fce_ellipj(1.234, 0.37, &sn, &cn, &dn);
    TEST_NEAR(sn * sn + cn * cn, 1.0, 1e-14);
    TEST_NEAR(dn * dn, 1.0 - 0.37 * sn * sn, 1e-14);
}

static void test_arc_jac_sn(void)
{
    /* asn(sn(u,m), m) = u for real u */
    {
        double u = 0.6, m = 0.3;
        double sn, cn, dn;
        fce_ellipj(u, m, &sn, &cn, &dn);
        {
            fce_cplx_t w = fce_arc_jac_sn(fce_cx(sn, 0.0), m);
            TEST_NEAR(w.re, u, 1e-12);
            TEST_NEAR(w.im, 0.0, 1e-12);
        }
    }
    /* scipy: asn(2j, 0.5) = 1.2190461596035318j */
    {
        fce_cplx_t w = fce_arc_jac_sn(fce_cx(0.0, 2.0), 0.5);
        TEST_NEAR(w.re, 0.0, 1e-12);
        TEST_NEAR(w.im, 1.2190461596035318, 1e-12);
    }
}

static void test_i0(void)
{
    TEST_NEAR(fce_i0(0.0), 1.0, 1e-15);
    TEST_NEAR(fce_i0(1.0), 1.2660658777520084, 1e-13);
    TEST_NEAR(fce_i0(5.0), 27.23987182360444, 1e-12);
    TEST_NEAR(fce_i0(10.0), 2815.716628466254, 1e-10);
}

static void test_kaiser(void)
{
    /* scipy.signal.kaiser_beta */
    TEST_NEAR(fce_kaiser_beta(65.0), 6.20426, 1e-4);
    TEST_NEAR(fce_kaiser_beta(80.0), 7.85726, 1e-4);
    TEST_NEAR(fce_kaiser_beta(30.0), 2.117, 1e-3);
    TEST_NEAR(fce_kaiser_beta(15.0), 0.0, 1e-15);
    /* scipy.signal.kaiserord: (167, 6.20426) for 65 dB, width 24/500 */
    {
        uint32_t taps = fce_kaiser_taps(65.0, 24.0, 1000.0);
        TEST_ASSERT(taps == 167u);
    }
    /* doc example: 80 dB, 1000 Hz transition at 48 kHz -> 243 taps */
    {
        uint32_t taps = fce_kaiser_taps(80.0, 1000.0, 48000.0);
        TEST_ASSERT(taps == 243u);
    }
}

static void test_windows(void)
{
    uint32_t n;
    /* window endpoints and symmetry */
    for (n = 2; n < 200; n += 7)
    {
        uint32_t i;
        for (i = 0; i < n; i++)
        {
            double a = fce_window_value(FCE_WIN_HANN, i, n, 0.0, 0.5);
            double b = fce_window_value(FCE_WIN_HANN, n - 1u - i, n, 0.0, 0.5);
            TEST_NEAR(a, b, 1e-14);
        }
        TEST_NEAR(fce_window_value(FCE_WIN_HANN, 0, n, 0.0, 0.5), 0.0, 1e-14);
        TEST_NEAR(fce_window_value(FCE_WIN_HAMMING, 0, n, 0.0, 0.5), 0.08,
                  1e-14);
        TEST_NEAR(fce_window_value(FCE_WIN_BLACKMAN, 0, n, 0.0, 0.5), 0.0,
                  1e-14);
        TEST_NEAR(fce_window_value(FCE_WIN_RECTANGULAR, 0, n, 0.0, 0.5), 1.0,
                  1e-15);
        TEST_NEAR(fce_window_value(FCE_WIN_BARTLETT, 0, n, 0.0, 0.5), 0.0,
                  1e-14);
        if (n & 1u)
            TEST_NEAR(fce_window_value(FCE_WIN_BARTLETT, (n - 1u) / 2u, n,
                                       0.0, 0.5), 1.0, 1e-14);
    }
    /* Kaiser window center = 1 (odd n only: exact center tap) */
    for (n = 3; n < 100; n += 6)
    {
        double w = fce_window_value(FCE_WIN_KAISER, (n - 1u) / 2u, n, 3.0,
                                    0.5);
        TEST_NEAR(w, 1.0, 1e-14);
    }
}

int main(void)
{
    test_ellipk();
    test_ellipj();
    test_arc_jac_sn();
    test_i0();
    test_kaiser();
    test_windows();
    printf("test_math: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
