/* example_03: FIR bandpass. fc1/fc2 are the two band edges. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_BANDPASS, 48000, 4000, 8000, 151,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("FIR Bandpass - Fs=48k, passband 4k..8k, 151 taps", &sp);
    ex_print_fir(&sp);
    return 0;
}
