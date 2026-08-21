/* common.h - shared helpers for the examples (host). */
#ifndef EXAMPLE_COMMON_H
#define EXAMPLE_COMMON_H

#include "filtercoeff.h"
#include <stdio.h>
#include <string.h>

#define EX_WS_SIZE (1u << 18)

/* design + print the markdown report to stdout */
static void ex_design_and_report(const char* title, fce_spec_t* sp)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;

    printf("\n===== %s =====\n", title);
    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_report(&r, &w);
    fwrite(buf, 1, mw.pos, stdout);
}

/* print the float coefficients of a FIR design */
static void ex_print_fir(fce_spec_t* sp)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    uint32_t i;
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;

    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_c_fir(&r, NULL, &w);
    printf("%s", buf);
    printf("(first 8 of %u taps: ", (unsigned)r.num_taps);
    for (i = 0; i < 8 && i < r.num_taps; i++)
        printf("%.6f ", r.h_f64[i]);
    printf("...)\n");
}

/* print the SOS of an IIR design */
static void ex_print_sos(fce_spec_t* sp)
{
    static uint8_t mem[EX_WS_SIZE];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_result_t r;
    fce_status_t st = fce_generate(sp, &r, &ws);
    char buf[16384];
    fce_mem_writer_t mw;
    fce_writer_t w;

    if (st != FCE_OK)
    {
        printf("ERROR: %s\n", fce_status_str(st));
        return;
    }
    fce_writer_mem_init(&w, &mw, buf, sizeof(buf));
    fce_export_c_sos(&r, NULL, &w);
    printf("%s", buf);
}

#endif /* EXAMPLE_COMMON_H */
