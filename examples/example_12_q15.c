/* example_12: Q15 fixed-point output for a FIR (symmetric scale), with
 * the full error report the library gives for free. */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;
    uint32_t i;

    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
                 FCE_WIN_KAISER, 0.0, FCE_PRECISION_FLOAT64);
    sp.stopband_atten_db = 60.0;
    sp.qformat = FCE_QFORMAT_Q15;

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed\n");
        return 1;
    }

    printf("\n===== example_12: FIR + Q15 =====\n");
    printf("scale             : %.3f  (q = round(c * scale))\n", r.scale);
    printf("integer bits      : %u\n", (unsigned)r.q_int_bits);
    printf("max abs error     : %.3g\n", r.q_max_abs_error);
    printf("rms error         : %.3g\n", r.q_rms_error);
    printf("max rel error     : %.3g\n", r.q_max_rel_error);
    printf("response error    : %.4f dB (float vs Q15)\n",
           r.quant_response_max_error_db);
    printf("first taps (Q15)  :");
    for (i = 0; i < 8u && i < r.num_taps; i++)
        printf(" %d", (int)r.q15[i]);
    printf(" ...\n");
    return 0;
}
