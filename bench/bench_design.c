/*
 * bench_design.c - design-time benchmarks (host).
 * Measures the wall-clock time of fce_generate for representative specs.
 * Design-time cost matters when coefficients are computed on the target
 * MCU at boot; it is irrelevant for host-side generation.
 */
#include "filtercoeff.h"
#include <stdio.h>
#include <time.h>

static uint8_t ws_mem[1 << 18];

static double bench(const char* name, fce_spec_t* sp, int iters)
{
    fce_workspace_t ws = { ws_mem, sizeof(ws_mem) };
    fce_result_t r;
    fce_status_t st;
    clock_t t0, t1;
    int i;
    double ms;

    /* warm-up */
    fce_generate(sp, &r, &ws);

    t0 = clock();
    for (i = 0; i < iters; i++)
    {
        st = fce_generate(sp, &r, &ws);
        if (st != FCE_OK)
        {
            printf("%-28s FAILED (%s)\n", name, fce_status_str(st));
            return 0.0;
        }
    }
    t1 = clock();
    ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC / (double)iters;
    printf("%-28s %10.3f ms/design\n", name, ms);
    return ms;
}

int main(void)
{
    fce_spec_t sp;

    printf("FilterCoeff design-time benchmark (host, %d Hz clock)\n",
           (int)CLOCKS_PER_SEC);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("FIR LP 101 taps (Kaiser)", &sp, 200);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR; sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
    sp.window = FCE_WIN_KAISER;
    sp.stopband_atten_db = 80; sp.transition_hz = 1000;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("FIR LP auto (243 taps)", &sp, 200);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 8;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("IIR Butterworth LP8", &sp, 500);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 2000; sp.fc2 = 4000; sp.order = 8;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 80;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("IIR Elliptic BP8", &sp, 200);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_BESSEL;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 8;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("IIR Bessel LP8", &sp, 200);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR; sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 8;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 80;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    bench("IIR Elliptic LP8 + Q15", &sp, 200);

    return 0;
}
