/* example_02: FIR highpass - note the odd tap count (Type I): an even-tap
 * highpass has a forced null at Nyquist and is flagged by the library. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_HIGHPASS, 48000, 2000, 0, 101,
                 FCE_WIN_HANN, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("FIR Highpass - Fs=48k, Fc=2k, 101 taps, Hann", &sp);
    ex_print_fir(&sp);
    return 0;
}
