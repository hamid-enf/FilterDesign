/* test_iir.c - IIR design: known filters, stability, auto-order. */
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

static void test_butter_lp2_closed_form(void)
{
    /* 2nd-order Butterworth LP at fc = fs/4 (warped):
     * scipy: b = [0.292893, 0.585786, 0.292893], a = [1, 0, 0.171573] */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 1000; sp.fc1 = 250; sp.order = 2;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_sections == 1);
    TEST_NEAR(r.sos_f64[0], 0.2928932188134525, 1e-12);
    TEST_NEAR(r.sos_f64[1], 0.5857864376269049, 1e-12);
    TEST_NEAR(r.sos_f64[2], 0.2928932188134525, 1e-12);
    TEST_NEAR(r.sos_f64[3], 0.0, 1e-14);
    TEST_NEAR(r.sos_f64[4], 0.1715728752538099, 1e-12);
    /* -3 dB at fc */
    TEST_NEAR(r.cutoff_measured_hz, 250.0, 1.0);
    TEST_ASSERT(r.max_pole_radius < 1.0);
}

static void test_all_families_stable(void)
{
    const fce_iir_family_t fams[5] = { FCE_IIR_BUTTERWORTH, FCE_IIR_CHEBYSHEV1,
                                       FCE_IIR_CHEBYSHEV2, FCE_IIR_ELLIPTIC,
                                       FCE_IIR_BESSEL };
    const fce_iir_type_t types[4] = { FCE_IIR_LOWPASS, FCE_IIR_HIGHPASS,
                                      FCE_IIR_BANDPASS, FCE_IIR_BANDSTOP };
    unsigned f, t;
    for (f = 0; f < 5; f++)
    {
        for (t = 0; t < 4; t++)
        {
            uint16_t order = (t < 2) ? (uint16_t)(3 + f) : (uint16_t)(2 + f);
            fce_spec_t sp;
            fce_result_t r;
            fce_spec_defaults(&sp);
            sp.kind = FCE_KIND_IIR;
            sp.iir_family = fams[f];
            sp.iir_type = types[t];
            sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 6000;
            sp.order = order;
            sp.passband_ripple_db = 0.5;
            sp.stopband_atten_db = 60;
            sp.precision = FCE_PRECISION_FLOAT64;
            TEST_OK(design(&sp, &r));
            TEST_ASSERT(r.max_pole_radius < 1.0);
            TEST_ASSERT(r.num_sections >= order / 2);
            TEST_ASSERT((r.flags & FCE_FLAG_UNSTABLE) == 0);
        }
    }
}

static void test_bessel_known_poles(void)
{
    /* 2nd-order Bessel (mag norm): H(s) = 3/(s^2 + 3s + 3) */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BESSEL;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 1000; sp.fc1 = 100; sp.order = 2;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_sections == 1);
    /* DC gain = 1 */
    {
        double h0 = (r.sos_f64[0] + r.sos_f64[1] + r.sos_f64[2]) /
                     (1.0 + r.sos_f64[3] + r.sos_f64[4]);
        TEST_NEAR(h0, 1.0, 1e-9);
    }
    /* -3 dB at fc */
    TEST_NEAR(r.cutoff_measured_hz, 100.0, 2.0);
}

static void test_auto_order_butter(void)
{
    /* the mission example: 1 kHz LP, Fs = 10 kHz, 60 dB at 2 kHz */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000; sp.order = 0;
    sp.edge1_hz = 2000;
    sp.passband_ripple_db = 3.0;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.order == 9); /* scipy buttord(1000,2000,3,60,fs=10k) */
    TEST_ASSERT(r.num_sections == 5);
    /* attenuation at 2 kHz must be >= 60 dB */
    {
        double w = 2.0 * FCE_PI * 2000.0 / 10000.0;
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
            double mb = sqrt(br * br + bi * bi);
            double ma = sqrt(ar * ar + ai * ai);
            g *= mb / ma;
        }
        TEST_ASSERT(-20.0 * log10(g) > 60.0);
    }
}

static void test_auto_order_ellip(void)
{
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000; sp.order = 0;
    sp.edge1_hz = 2000;
    sp.passband_ripple_db = 1.0;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.order >= 4 && r.order <= 6);
    TEST_ASSERT(r.max_pole_radius < 1.0);
}

static void test_iir_known_highpass(void)
{
    /* 1st-order Butterworth HP at fc = fs/4: the warped analog cutoff is
     * exactly 2*fs, so the pole lands at z = 0:
     * H(z) = 0.5 (1 - z^-1)  (scipy: b=[0.5,-0.5], a=[1,0]) */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_HIGHPASS;
    sp.fs = 1000; sp.fc1 = 250; sp.order = 1;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.num_sections == 1);
    TEST_NEAR(r.sos_f64[0], 0.5, 1e-12);
    TEST_NEAR(r.sos_f64[1], -0.5, 1e-12);
    TEST_NEAR(r.sos_f64[2], 0.0, 1e-14);
    TEST_NEAR(r.sos_f64[3], 0.0, 1e-14);
    TEST_NEAR(r.sos_f64[4], 0.0, 1e-14);
}

static void test_sos_sign_convention(void)
{
    /* verify the documented difference equation convention: simulate an
     * impulse and check y matches the direct-form recurrence */
    fce_spec_t sp;
    fce_result_t r;
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    {
        double x[16];
        double y[16];
        double state[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        uint32_t i;
        memset(x, 0, sizeof(x));
        x[0] = 1.0;
        fce_sim_sos(r.sos_f64, r.num_sections, x, y, 16, state);
        /* y[0] = b0 of section 0 (product with later sections = 1 at n=0) */
        TEST_NEAR(y[0], r.sos_f64[0], 1e-12);
        /* direct-form check for the first section alone */
        {
            double b0 = r.sos_f64[0], b1 = r.sos_f64[1], b2 = r.sos_f64[2];
            double a1 = r.sos_f64[3], a2 = r.sos_f64[4];
            double yd[16] = { 0 };
            for (i = 0; i < 16; i++)
            {
                double xi = (i == 0) ? 1.0 : 0.0;
                double yi = b0 * xi;
                if (i >= 1) yi += b1 * x[i - 1] - a1 * yd[i - 1];
                if (i >= 2) yi += b2 * x[i - 2] - a2 * yd[i - 2];
                yd[i] = yi;
            }
            /* feed the first section alone through sim */
            {
                double ys[16] = { 0 };
                double st2[2] = { 0, 0 };
                fce_sim_sos(r.sos_f64, 1, x, ys, 16, st2);
                for (i = 0; i < 16; i++)
                    TEST_NEAR(ys[i], yd[i], 1e-12);
            }
        }
    }
}

static void test_sos_ordering(void)
{
    fce_spec_t sp;
    fce_result_t r;
    double radius_asc[8], radius_desc[8];
    uint32_t s;

    /* ascending (default) */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 2000; sp.order = 8;
    sp.sos_order = FCE_SOS_ORDER_DEFAULT;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    for (s = 0; s < r.num_sections; s++)
    {
        double a2 = r.sos_f64[5 * s + 4];
        radius_asc[s] = sqrt(fabs(a2));
    }
    /* descending */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 2000; sp.order = 8;
    sp.sos_order = FCE_SOS_ORDER_POLE_RADIUS_DESC;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    for (s = 0; s < r.num_sections; s++)
    {
        double a2 = r.sos_f64[5 * s + 4];
        radius_desc[s] = sqrt(fabs(a2));
    }
    for (s = 0; s + 1 < r.num_sections; s++)
    {
        TEST_ASSERT(radius_asc[s] <= radius_asc[s + 1] + 1e-12);
        TEST_ASSERT(radius_desc[s] >= radius_desc[s + 1] - 1e-12);
    }
    /* the two orderings are permutations of each other */
    {
        int match = 1;
        for (s = 0; s < r.num_sections; s++)
        {
            uint32_t t;
            int found = 0;
            for (t = 0; t < r.num_sections; t++)
            {
                if (fabs(radius_desc[s] - radius_asc[t]) < 1e-12)
                    found = 1;
            }
            if (!found)
                match = 0;
        }
        TEST_ASSERT(match);
    }
}

static void test_high_order(void)
{
    fce_spec_t sp;
    fce_result_t r;
    /* order 16 butterworth LP (near the practical limit) */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 16;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.max_pole_radius < 1.0);
    TEST_ASSERT(r.num_sections == 8);

    /* high-order elliptic */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 10;
    sp.passband_ripple_db = 0.1; sp.stopband_atten_db = 80;
    sp.precision = FCE_PRECISION_FLOAT64;
    TEST_OK(design(&sp, &r));
    TEST_ASSERT(r.max_pole_radius < 1.0);
}

int main(void)
{
    test_butter_lp2_closed_form();
    test_all_families_stable();
    test_bessel_known_poles();
    test_auto_order_butter();
    test_auto_order_ellip();
    test_iir_known_highpass();
    test_sos_sign_convention();
    test_sos_ordering();
    test_high_order();
    printf("test_iir: %d run, %d failed\n", g_tests_run, g_tests_failed);
    return g_tests_failed ? 1 : 0;
}
