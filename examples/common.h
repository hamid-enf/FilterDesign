/* common.h - shared helpers for the examples (host). */
#ifndef EXAMPLE_COMMON_H
#define EXAMPLE_COMMON_H

#include "filtercoeff.h"
#include <stdio.h>
#include <string.h>

/*
 * The examples call the design helpers sequentially, so they only need one
 * workspace.  Keep the size configurable for embedded targets, for example:
 *
 *     -DEX_WS_SIZE=8192
 *
 * For a 101-tap FIR, fce_workspace_required() reports about 4.3 KiB.  The
 * default below covers fce_workspace_required_max() with the repository's
 * default limits (86,016 bytes) without reserving an unnecessarily large
 * quarter megabyte on an MCU.
 */
#ifndef EX_WS_SIZE
#define EX_WS_SIZE (96u * 1024u)
#endif

/* The report/export helpers are also used by embedded examples. Keep their
 * temporary output out of the target stack; override this if a larger export
 * is required. */
#ifndef EX_OUTPUT_SIZE
#define EX_OUTPUT_SIZE (16u * 1024u)
#endif

#if defined(__GNUC__) || defined(__clang__)
#define EX_UNUSED __attribute__((unused))
#define EX_ALIGN8 __attribute__((aligned(8)))
#else
#define EX_UNUSED
#define EX_ALIGN8
#endif

/* One shared buffer avoids a separate 256 KiB .bss allocation per helper. */
static uint8_t ex_workspace_mem[EX_WS_SIZE] EX_ALIGN8;
static char ex_output_buf[EX_OUTPUT_SIZE];

/* design + print the markdown report to stdout */
static EX_UNUSED void ex_design_and_report(const char* title, fce_spec_t* sp)
{
    fce_workspace_t ws = { ex_workspace_mem, sizeof(ex_workspace_mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    fce_mem_writer_t mw;
    fce_writer_t w;

    printf("\n===== %s =====\n", title);
    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, ex_output_buf, sizeof(ex_output_buf));
    fce_export_report(&r, &w);
    fwrite(ex_output_buf, 1, mw.pos, stdout);
}

/* print the float coefficients of a FIR design */
static EX_UNUSED void ex_print_fir(fce_spec_t* sp)
{
    fce_workspace_t ws = { ex_workspace_mem, sizeof(ex_workspace_mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    uint32_t i;
    fce_mem_writer_t mw;
    fce_writer_t w;

    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, ex_output_buf, sizeof(ex_output_buf));
    fce_export_c_fir(&r, NULL, &w);
    printf("%s", ex_output_buf);
    printf("(first 8 of %u taps: ", (unsigned)r.num_taps);
    for (i = 0; i < 8 && i < r.num_taps; i++)
        printf("%.6f ", r.h_f64[i]);
    printf("...)\n");
}

/* print the SOS of an IIR design */
static EX_UNUSED void ex_print_sos(fce_spec_t* sp)
{
    fce_workspace_t ws = { ex_workspace_mem, sizeof(ex_workspace_mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    fce_mem_writer_t mw;
    fce_writer_t w;

    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, ex_output_buf, sizeof(ex_output_buf));
    fce_export_c_sos(&r, NULL, &w);
    printf("%s", ex_output_buf);
}

#endif /* EXAMPLE_COMMON_H */
