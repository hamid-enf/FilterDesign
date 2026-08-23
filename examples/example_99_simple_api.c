/* example_99: the "simple API" - plain structs + one call, ready to hand
 * the coefficients to your own runtime (here: a fictional FilterLab). */
#include "common.h"

int main(void)
{
    FilterCoeffSpec spec;
    FilterCoeffResult result;
    static uint8_t mem[65536];
    FilterCoeffWorkspace ws = { mem, sizeof(mem) };

    fce_spec_defaults(&spec);
    spec.kind = FCE_KIND_FIR;
    spec.fir_type = FCE_FIR_LOWPASS;
    spec.fs = 48000;
    spec.fc1 = 5000;
    spec.num_taps = 101;
    spec.window = FCE_WIN_KAISER;
    spec.stopband_atten_db = 80;     /* beta is picked automatically */
    spec.precision = FCE_PRECISION_FLOAT32;

    printf("\n===== example_99: simple API =====\n");
    if (FilterCoeff_Generate(&spec, &result, &ws) == FCE_OK)
    {
        /* result.h_f32 -> FilterLab_FIR_SetCoefficients(...) */
        printf("OK: %u taps, beta %.3f, first tap %.8g\n",
               (unsigned)result.num_taps, result.kaiser_beta,
               (double)result.h_f32[0]);
        printf("Status: %s, flags 0x%02X\n",
               fce_status_str(result.status), (unsigned)result.flags);
    }
    else
    {
        printf("ERROR: %s\n", fce_status_str(result.status));
        return 1;
    }
    return 0;
}
