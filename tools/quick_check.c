/* quick numerical check vs scipy: prints JSON-ish output for a set of designs */
#include "filtercoeff.h"
#include <stdio.h>
#include <string.h>

static void print_fir(const char* label, fce_spec_t* sp)
{
    fce_result_t r;
    uint8_t mem[1 << 16];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_status_t st = fce_generate(sp, &r, &ws);
    uint32_t i;
    printf("FIR %s status=%d taps=%u flags=%x\n", label, (int)st,
           (unsigned)r.num_taps, (unsigned)r.flags);
    if (st != FCE_OK)
        return;
    for (i = 0; i < r.num_taps; i++)
        printf("%.17g%s", r.h_f64[i], (i + 1 < r.num_taps) ? "," : "");
    printf("\n");
}

static void print_iir(const char* label, fce_spec_t* sp)
{
    fce_result_t r;
    uint8_t mem[1 << 16];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_status_t st = fce_generate(sp, &r, &ws);
    uint32_t i, j;
    printf("IIR %s status=%d order=%u sec=%u flags=%x mr=%.12g\n", label,
           (int)st, (unsigned)r.order, (unsigned)r.num_sections,
           (unsigned)r.flags, r.max_pole_radius);
    if (st != FCE_OK)
        return;
    for (i = 0; i < r.num_sections; i++)
    {
        printf("  [");
        for (j = 0; j < 5; j++)
            printf("%.17g%s", r.sos_f64[5 * i + j], (j < 4) ? "," : "");
        printf("]\n");
    }
}

int main(void)
{
    fce_spec_t sp;

    /* FIR */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_fir("lp101-kaiser", &sp);

    /* FIR auto taps */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 80; sp.transition_hz = 1000;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_fir("lp-auto-kaiser", &sp);

    /* IIR butter 4th order LP */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("butter4-lp", &sp);

    /* IIR cheby1 5th order HP */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_CHEBYSHEV1;
    sp.iir_type = FCE_IIR_HIGHPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 5;
    sp.passband_ripple_db = 1.0;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("cheby1-5-hp", &sp);

    /* IIR cheby2 6th order BP */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_CHEBYSHEV2;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 6000; sp.order = 6;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("cheby2-6-bp", &sp);

    /* IIR elliptic 5th order BS */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_BANDSTOP;
    sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 6000; sp.order = 5;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("ellip-5-bs", &sp);

    /* IIR bessel 4th order LP */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BESSEL;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("bessel-4-lp", &sp);

    /* auto order butter */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000; sp.order = 0;
    sp.edge1_hz = 2000; sp.passband_ripple_db = 3.0;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("butter-auto-lp", &sp);

    /* auto order ellip */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000; sp.order = 0;
    sp.edge1_hz = 2000; sp.passband_ripple_db = 1.0;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    print_iir("ellip-auto-lp", &sp);

    return 0;
}
