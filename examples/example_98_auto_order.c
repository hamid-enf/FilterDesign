/* example_98: automatic IIR order (the buttord/cheb1ord/ellipord
 * equivalent). Mission: LP 1 kHz, Fs = 10 kHz, >= 60 dB from 2 kHz. */
#include "common.h"

static double g_edge_db = 0.0;
static bool edge_probe(void* ctx, const fce_response_point_t* pt)
{
    (void)ctx;
    g_edge_db = pt->mag_db;
    return false;
}

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;
    fce_status_t st;

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000;      /* passband edge  */
    sp.edge1_hz = 2000;                /* stopband edge  */
    sp.stopband_atten_db = 60;         /* >= 60 dB there */
    sp.order = 0;                      /* automatic      */
    sp.precision = FCE_PRECISION_FLOAT64;

    st = fce_generate(&sp, &r, &ws);

    printf("\n===== example_98: auto order =====\n");
    printf("mission   : LP 1 kHz, Fs=10 kHz, >=60 dB from 2 kHz\n");
    if (st != FCE_OK)
    {
        printf("design failed: %s\n", fce_status_str(st));
        return 1;
    }
    printf("Order     : %u  (scipy.buttord gives the same)\n",
           (unsigned)r.order);
    printf("Sections  : %u\n", (unsigned)r.num_sections);

    /* verify the mission: gain at the 2 kHz stopband edge */
    fce_response_sos(r.sos_f64, r.num_sections, r.fs, 1, 2000.0, 2000.0,
                     edge_probe, NULL);
    printf("Gain@2kHz : %.1f dB (must be <= -60)\n", g_edge_db);
    printf("Validation: %s\n",
           (g_edge_db <= -60.0 && st == FCE_OK) ? "PASS" : "FAIL");
    return 0;
}
