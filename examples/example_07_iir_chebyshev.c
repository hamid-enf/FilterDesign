/* example_07: Chebyshev I (passband ripple) and Chebyshev II (stopband
 * ripple) side by side - both 5th order lowpass at 5 kHz. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;

    fce_spec_iir(&sp, FCE_IIR_CHEBYSHEV1, FCE_IIR_LOWPASS,
                 48000, 5000, 0, 5, 1.0 /* rp = 1 dB */, 0.0,
                 FCE_PRECISION_FLOAT64);
    ex_design_and_report("IIR Chebyshev I LP - rp=1 dB, order 5", &sp);
    ex_print_sos(&sp);

    fce_spec_iir(&sp, FCE_IIR_CHEBYSHEV2, FCE_IIR_LOWPASS,
                 48000, 5000, 0, 5, 0.0, 60.0 /* rs = 60 dB */,
                 FCE_PRECISION_FLOAT64);
    ex_design_and_report("IIR Chebyshev II LP - rs=60 dB, order 5", &sp);
    ex_print_sos(&sp);
    return 0;
}
