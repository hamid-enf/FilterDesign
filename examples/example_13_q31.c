/* example_13: Q31 fixed-point SOS with per-section scaling
 * (FCE_SCALE_SECTION_WISE): each biquad gets its own scale so the small
 * coefficients keep their precision. */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;
    uint32_t s;

    fce_spec_iir(&sp, FCE_IIR_ELLIPTIC, FCE_IIR_LOWPASS,
                 48000, 5000, 0, 6, 0.5, 80.0, FCE_PRECISION_FLOAT64);
    sp.qformat = FCE_QFORMAT_Q31;
    sp.scale_strategy = FCE_SCALE_SECTION_WISE;

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed\n");
        return 1;
    }

    printf("\n===== example_13: IIR + Q31 (section-wise) =====\n");
    for (s = 0; s < r.num_sections; s++)
    {
        const int32_t* q = r.q31 + 5u * s;
        printf("section %u: scale=%.0f  b=[%ld %ld %ld] a=[%ld %ld]\n",
               (unsigned)s, r.section_scales[s],
               (long)q[0], (long)q[1], (long)q[2], (long)q[3], (long)q[4]);
    }
    /* reconstructed section must be scaled back: c = q / section_scale */
    printf("response error: %.4f dB, quantized pole radius: %.9f\n",
           r.quant_response_max_error_db, r.quant_max_pole_radius);
    return 0;
}
