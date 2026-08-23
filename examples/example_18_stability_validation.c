/* example_18: stability validation - pole radii from the SOS, plus a
 * response scan through the streaming callback API. */
#include "common.h"

static uint32_t g_count = 0;
static bool print_every_64th(void* ctx, const fce_response_point_t* pt)
{
    (void)ctx;
    if ((g_count++ % 64u) == 0u)
        printf("  %8.0f Hz : %7.2f dB, GD %.1f samples\n",
               pt->f_hz, pt->mag_db, pt->group_delay);
    return true;
}

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;
    double mr, margin;

    fce_spec_iir(&sp, FCE_IIR_CHEBYSHEV1, FCE_IIR_HIGHPASS,
                 48000, 6000, 0, 8, 0.5, 0.0, FCE_PRECISION_FLOAT64);

    printf("\n===== example_18: stability validation =====\n");
    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed: %s\n", fce_status_str(r.status));
        return 1;
    }

    /* standalone stability check on the SOS (also done internally) */
    if (fce_stability_sos(r.sos_f64, r.num_sections, &mr, &margin) == FCE_OK)
        printf("STABLE: max pole radius %.6f, margin %.3g\n", mr, margin);
    else
        printf("UNSTABLE (this should not happen for a valid design)\n");

    /* transparency: the digital poles themselves are exposed */
    printf("digital poles/zeros exposed: %u / %u\n",
           (unsigned)r.iir_npoles, (unsigned)r.iir_nzeros);

    printf("response scan:\n");
    fce_response_sos(r.sos_f64, r.num_sections, r.fs, 512, 0.0, 24000.0,
                     print_every_64th, NULL);
    return 0;
}
