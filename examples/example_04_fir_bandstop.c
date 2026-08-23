/* example_04: FIR bandstop (notch band 1 kHz .. 3 kHz). */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_BANDSTOP, 48000, 1000, 3000, 121,
                 FCE_WIN_BLACKMAN, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("FIR Bandstop - Fs=48k, stop 1k..3k, 121 taps, Blackman", &sp);
    ex_print_fir(&sp);
    return 0;
}
