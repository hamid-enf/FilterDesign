/* example_06: IIR Butterworth lowpass (order 6), SOS output. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_iir(&sp, FCE_IIR_BUTTERWORTH, FCE_IIR_LOWPASS,
                 48000, 5000, 0, 6, 0.0, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("IIR Butterworth LP - Fs=48k, Fc=5k, order 6", &sp);
    ex_print_sos(&sp);
    return 0;
}
