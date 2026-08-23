/* example_14: exporting SOS to every format: C array, CSV, JSON, and the
 * Markdown report. Everything streams through fce_writer_t. */
#include "common.h"

int main(void)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_spec_t sp;
    static char buf[32768];
    fce_mem_writer_t mw;
    fce_writer_t w;

    fce_spec_iir(&sp, FCE_IIR_ELLIPTIC, FCE_IIR_BANDPASS,
                 48000, 2000, 4000, 6, 1.0, 50.0, FCE_PRECISION_FLOAT32);

    if (fce_generate(&sp, &r, &ws) != FCE_OK)
    {
        printf("design failed\n");
        return 1;
    }

    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_c_sos(&r, NULL, &w);
    printf("----- C array -----\n%s\n", buf);

    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_csv(&r, &w);
    printf("----- CSV -----\n%s\n", buf);

    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_json(&r, &w);
    printf("----- JSON -----\n%s\n", buf);
    return 0;
}
