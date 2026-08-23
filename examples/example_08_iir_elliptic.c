/* example_08: IIR Elliptic (Cauer) - sharpest transition per order. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;

    fce_spec_iir(&sp, FCE_IIR_ELLIPTIC, FCE_IIR_LOWPASS,
                 48000, 5000, 0, 6, 0.5 /* rp */, 80.0 /* rs */,
                 FCE_PRECISION_FLOAT64);
    ex_design_and_report("IIR Elliptic LP - rp=0.5 dB, rs=80 dB, order 6", &sp);
    ex_print_sos(&sp);

    /* elliptic bandpass: 6th-order prototype -> 6 SOS sections */
    fce_spec_iir(&sp, FCE_IIR_ELLIPTIC, FCE_IIR_BANDPASS,
                 48000, 4000, 8000, 6, 0.5, 60.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("IIR Elliptic BP - 4k..8k, order 6", &sp);
    ex_print_sos(&sp);
    return 0;
}
