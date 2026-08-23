/* example_01: FIR lowpass (Hamming window) - the hello-world design. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("FIR Lowpass - Fs=48k, Fc=5k, 101 taps, Hamming", &sp);
    ex_print_fir(&sp);
    return 0;
}
