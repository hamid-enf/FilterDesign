/* example_09: IIR Bessel - maximally flat group delay (linear phase,
 * gentle magnitude). Shows the group-delay feedback point. */
#include "common.h"

static double g_dc_gd = 0.0;
static bool gd_probe(void* ctx, const fce_response_point_t* pt)
{
    (void)ctx;
    g_dc_gd = pt->group_delay;
    return false;
}

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;

    fce_spec_iir(&sp, FCE_IIR_BESSEL, FCE_IIR_LOWPASS,
                 48000, 2000, 0, 6, 0.0, 0.0, FCE_PRECISION_FLOAT64);
    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed\n");
        return 1;
    }
    ex_design_and_report("IIR Bessel LP - Fs=48k, Fc=2k, order 6", &sp);

    /* group delay at DC (samples) - the Bessel claim to fame */
    fce_response_sos(r.sos_f64, r.num_sections, r.fs, 1, 0.0, 0.0,
                     gd_probe, NULL);
    printf("Group delay at DC: %.2f samples (flat across the passband)\n",
           g_dc_gd);
    return 0;
}
