/* example_15: STM32H7 integration pattern (host-runnable demo).
 *
 * On the target you would design ONCE at init (no dynamic allocation,
 * workspace in fast RAM) and hand the arrays to CMSIS-DSP:
 *
 *   __attribute__((section(".ram_d1"))) static uint8_t ws_mem[4096];
 *
 *   arm_fir_instance_f32 S;
 *   static float32_t state[NUM_TAPS + BLOCK - 1];
 *   fce_generate(&sp, &r, &ws);
 *   arm_fir_init_f32(&S, r.num_taps, (float32_t*)r.h_f32, state, BLOCK);
 *   arm_fir_f32(&S, dma_in, dma_out, BLOCK);
 */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;

    /* ask exactly how much workspace the spec needs (could be static too) */
    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT32);
    printf("\n===== example_15: STM32H7 pattern =====\n");
    printf("workspace required: %lu bytes\n",
           (unsigned long)fce_workspace_required(&sp));

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(r.status));
        return 1;
    }
    printf("designed: %u taps, float32 -> feed r.h_f32 to arm_fir_init_f32\n",
           (unsigned)r.num_taps);
    printf("DC gain %.4f dB, stable output ready for DMA chunks\n",
           r.dc_gain_db);
    return 0;
}
