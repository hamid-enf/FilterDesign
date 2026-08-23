/* example_11: float64 output - full 17-digit coefficients for offline /
 * host use (also the upper bound for accuracy comparisons). */
#include "common.h"

int main(void)
{
    fce_spec_t sp;
    fce_spec_fir(&sp, FCE_FIR_LOWPASS, 16000, 1000, 0, 63,
                 FCE_WIN_BLACKMAN_HARRIS, 0.0, FCE_PRECISION_FLOAT64);
    ex_design_and_report("FIR LP / float64 output (Blackman-Harris)", &sp);
    ex_print_fir(&sp); /* exports a static const float64_t array, %.17g */
    return 0;
}
