/* example_16: ESP32-C2 integration pattern (host-runnable demo).
 *
 * The C2 (RISC-V, small RAM, no FPU) wants Q15/Q31 coefficients; design at
 * init or pre-generate the arrays on the host and ship them in flash:
 *
 *   fce_generate(&sp, &r, &ws);          // Q15 SOS
 *   r.q15 / r.section_scales -> dsps_biquad_f32-style or your own Q15 loop
 */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;

    fce_spec_iir(&sp, FCE_IIR_BUTTERWORTH, FCE_IIR_LOWPASS,
                 16000, 2000, 0, 4, 0.0, 0.0, FCE_PRECISION_FLOAT32);
    sp.qformat = FCE_QFORMAT_Q15;

    printf("\n===== example_16: ESP32-C2 pattern =====\n");
    printf("workspace required: %lu bytes\n",
           (unsigned long)fce_workspace_required(&sp));

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(r.status));
        return 1;
    }
    printf("designed: order %u, %u sections, Q15 (scale %.1f)\n",
           (unsigned)r.order, (unsigned)r.num_sections, r.scale);
    printf("quantized poles stay at radius %.6f (< 1: stable on target)\n",
           r.quant_max_pole_radius);
    ex_print_sos(&sp);
    return 0;
}
