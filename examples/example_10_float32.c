/* example_10: float32 output - internal math stays float64; float32 is
 * only the final rounding, ready for a float32 FIR runtime. */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
                 FCE_WIN_HAMMING, 0.0, FCE_PRECISION_FLOAT32);
    ex_design_and_report("FIR LP / float32 output", &sp);
    ex_print_fir(&sp); /* exports a static const float32_t array */
    return 0;
}
