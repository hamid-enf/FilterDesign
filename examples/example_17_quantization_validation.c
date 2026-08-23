/* example_17: quantization validation - how much damage does Q15 do?
 * The library compares the quantized response against the float design
 * and re-checks the pole positions after quantization. */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;

    /* a sharp elliptic lowpass is the worst case for Q15: poles near |z|=1 */
    fce_spec_iir(&sp, FCE_IIR_ELLIPTIC, FCE_IIR_LOWPASS,
                 48000, 4000, 0, 8, 0.1, 100.0, FCE_PRECISION_FLOAT64);
    sp.qformat = FCE_QFORMAT_Q15;

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed: %s\n", fce_status_str(r.status));
        return 1;
    }

    printf("\n===== example_17: quantization validation =====\n");
    printf("float design      : pole radius %.9f (margin %.3g)\n",
           r.max_pole_radius, r.stability_margin);
    printf("quantized (Q15)   : pole radius %.9f (margin %.3g)\n",
           r.quant_max_pole_radius, r.quant_stability_margin);
    printf("coefficient errors: max %.3g, rms %.3g, max rel %.3g\n",
           r.q_max_abs_error, r.q_rms_error, r.q_max_rel_error);
    printf("response error    : %.4f dB (passband/transition)\n",
           r.quant_response_max_error_db);
    printf("flags             : 0x%02X%s\n", (unsigned)r.flags,
           (r.flags & FCE_FLAG_QUANTIZATION_WARNING) ? " QUANT_WARNING" : "");
    return 0;
}
