/* example_05: fully automatic Kaiser design: give the attenuation and the
 * transition width, get taps + beta back. Mission: 80 dB, 1 kHz transition. */
#include "common.h"

static double g_stop_db = 0.0;
static bool stopband_probe(void* ctx, const fce_response_point_t* pt)
{
    (void)ctx;
    g_stop_db = pt->mag_db;
    return false; /* first point only */
}

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;

    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 48000, 5000, 0, 0 /* auto taps */,
                 FCE_WIN_KAISER, 0.0 /* auto beta */, FCE_PRECISION_FLOAT64);
    sp.stopband_atten_db = 80.0;
    sp.transition_hz = 1000.0;

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed\n");
        return 1;
    }

    /* measured gain at the stopband edge (fc1 + transition = 6 kHz) */
    fce_response_fir(r.h_f64, r.num_taps, r.fs, 1, 6000.0, 6000.0,
                     stopband_probe, NULL);

    printf("\n===== example_05: auto Kaiser =====\n");
    printf("Taps      : %u\n", (unsigned)r.num_taps);
    printf("Kaiser b  : %.4f\n", r.kaiser_beta);
    printf("Stopband  : %.1f dB   (requested 80)\n", -g_stop_db);
    printf("Flags     : 0x%02X\n", (unsigned)r.flags);
    return 0;
}
