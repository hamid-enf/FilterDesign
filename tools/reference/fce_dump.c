/*
 * fce_dump.c - reference-test dump tool.
 *
 * Designs a fixed matrix of filters and prints one JSON object per line:
 *   { "id": ..., "kind": "FIR"|"IIR", "spec": {...}, "coeffs": [...],
 *     "sos": [[...]], "order": n, "taps": n, "beta": x, "norm_factor": x,
 *     "status": n, "flags": n, "design_fc1": x, "design_fc2": x }
 *
 * The Python reference tool (compare_reference.py) reads this and compares
 * against SciPy. This tool is host-only (uses stdio).
 */
#include "filtercoeff.h"
#include <stdio.h>
#include <string.h>

static void print_json_str(const char* s)
{
    putchar('"');
    for (; *s; s++)
    {
        if (*s == '"' || *s == '\\')
            putchar('\\');
        putchar(*s);
    }
    putchar('"');
}

static void dump(const char* id, fce_spec_t* sp)
{
    fce_result_t r;
    uint8_t mem[1 << 17];
    fce_workspace_t ws = { mem, sizeof(mem) };
    fce_status_t st = fce_generate(sp, &r, &ws);
    uint32_t i, j;

    printf("{\"id\":");
    print_json_str(id);
    printf(",\"kind\":\"%s\"", r.kind == FCE_KIND_FIR ? "FIR" : "IIR");
    if (r.kind == FCE_KIND_FIR)
    {
        printf(",\"type\":\"%s\"",
               r.fir_type == FCE_FIR_LOWPASS ? "Lowpass" :
               r.fir_type == FCE_FIR_HIGHPASS ? "Highpass" :
               r.fir_type == FCE_FIR_BANDPASS ? "Bandpass" :
               r.fir_type == FCE_FIR_BANDSTOP ? "Bandstop" :
               r.fir_type == FCE_FIR_HILBERT ? "Hilbert" : "Differentiator");
        printf(",\"window\":\"%s\"",
               r.window == FCE_WIN_RECTANGULAR ? "Rectangular" :
               r.window == FCE_WIN_HANN ? "Hann" :
               r.window == FCE_WIN_HAMMING ? "Hamming" :
               r.window == FCE_WIN_BLACKMAN ? "Blackman" :
               r.window == FCE_WIN_KAISER ? "Kaiser" :
               r.window == FCE_WIN_BLACKMAN_HARRIS ? "Blackman-Harris" :
               r.window == FCE_WIN_BARTLETT ? "Bartlett" : "Tukey");
        printf(",\"norm\":\"%s\"",
               r.normalization == FCE_NORM_DC ? "DC" :
               r.normalization == FCE_NORM_NYQUIST ? "Nyquist" :
               r.normalization == FCE_NORM_PASSBAND_PEAK ? "Peak" : "None");
    }
    else
    {
        printf(",\"type\":\"%s\"",
               r.iir_type == FCE_IIR_LOWPASS ? "Lowpass" :
               r.iir_type == FCE_IIR_HIGHPASS ? "Highpass" :
               r.iir_type == FCE_IIR_BANDPASS ? "Bandpass" : "Bandstop");
        printf(",\"family\":\"%s\"",
               r.iir_family == FCE_IIR_BUTTERWORTH ? "Butterworth" :
               r.iir_family == FCE_IIR_CHEBYSHEV1 ? "Chebyshev I" :
               r.iir_family == FCE_IIR_CHEBYSHEV2 ? "Chebyshev II" :
               r.iir_family == FCE_IIR_ELLIPTIC ? "Elliptic" : "Bessel");
        printf(",\"rp\":%.17g,\"rs\":%.17g", r.passband_ripple_db,
               r.stopband_atten_db);
    }
    printf(",\"fs\":%.17g,\"fc1\":%.17g,\"fc2\":%.17g", r.fs, r.fc1, r.fc2);
    printf(",\"design_fc1\":%.17g,\"design_fc2\":%.17g",
           r.design_fc1, r.design_fc2);
    printf(",\"status\":%d,\"flags\":%u", (int)st, (unsigned)r.flags);
    printf(",\"precision\":\"%s\"",
           r.precision == FCE_PRECISION_FLOAT32 ? "float32" : "float64");
    printf(",\"norm_factor\":%.17g", r.norm_factor);
    printf(",\"order\":%u,\"taps\":%u,\"sections\":%u",
           (unsigned)r.order, (unsigned)r.num_taps, (unsigned)r.num_sections);
    printf(",\"beta\":%.17g", r.kaiser_beta);
    printf(",\"max_pole_radius\":%.17g", r.max_pole_radius);
    printf(",\"stopband_atten_measured_db\":%.17g",
           r.stopband_atten_measured_db);

    if (r.kind == FCE_KIND_FIR && r.h_f64 != NULL)
    {
        printf(",\"coeffs\":[");
        for (i = 0; i < r.num_taps; i++)
            printf("%.17g%s", r.h_f64[i], (i + 1 < r.num_taps) ? "," : "");
        printf("]");
    }
    if (r.kind == FCE_KIND_IIR && r.sos_f64 != NULL)
    {
        printf(",\"sos\":[");
        for (i = 0; i < r.num_sections; i++)
        {
            printf("[");
            for (j = 0; j < 5; j++)
                printf("%.17g%s", r.sos_f64[5 * i + j],
                       (j < 4) ? "," : "");
            printf("]%s", (i + 1 < r.num_sections) ? "," : "");
        }
        printf("]");
    }
    if (r.q15 != NULL)
    {
        uint32_t n = (r.kind == FCE_KIND_FIR) ? r.num_taps
                                              : 5u * r.num_sections;
        printf(",\"q15\":[");
        for (i = 0; i < n; i++)
            printf("%d%s", (int)r.q15[i], (i + 1 < n) ? "," : "");
        printf("]");
    }
    if (r.q31 != NULL)
    {
        uint32_t n = (r.kind == FCE_KIND_FIR) ? r.num_taps
                                              : 5u * r.num_sections;
        printf(",\"q31\":[");
        for (i = 0; i < n; i++)
            printf("%ld%s", (long)r.q31[i], (i + 1 < n) ? "," : "");
        printf("]");
    }
    printf(",\"scale\":%.17g,\"q_int_bits\":%u", r.scale,
           (unsigned)r.q_int_bits);
    printf(",\"q_max_abs_error\":%.17g,\"q_rms_error\":%.17g",
           r.q_max_abs_error, r.q_rms_error);
    printf(",\"q_resp_err_db\":%.17g", r.quant_response_max_error_db);
    printf("}\n");
}

int main(void)
{
    fce_spec_t sp;

    /* ---------------- FIR ---------------- */
#define FIR_CASE(id_, type_, fs_, fc1_, fc2_, taps_, win_, beta_, norm_) \
    do { \
        fce_spec_defaults(&sp); \
        sp.kind = FCE_KIND_FIR; \
        sp.fir_type = (type_); sp.fs = (fs_); sp.fc1 = (fc1_); \
        sp.fc2 = (fc2_); sp.num_taps = (taps_); sp.window = (win_); \
        sp.kaiser_beta = (beta_); sp.normalization = (norm_); \
        sp.precision = FCE_PRECISION_FLOAT64; \
        dump(id_, &sp); \
    } while (0)

    FIR_CASE("fir_lp_hann_odd", FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
             FCE_WIN_HANN, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_hamming_even", FCE_FIR_LOWPASS, 8000, 1000, 0, 64,
             FCE_WIN_HAMMING, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_blackman", FCE_FIR_LOWPASS, 44100, 8000, 0, 151,
             FCE_WIN_BLACKMAN, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_kaiser", FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
             FCE_WIN_KAISER, 7.86, FCE_NORM_DC);
    FIR_CASE("fir_lp_bh", FCE_FIR_LOWPASS, 48000, 12000, 0, 129,
             FCE_WIN_BLACKMAN_HARRIS, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_bartlett", FCE_FIR_LOWPASS, 16000, 2000, 0, 81,
             FCE_WIN_BARTLETT, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_tukey", FCE_FIR_LOWPASS, 48000, 5000, 0, 101,
             FCE_WIN_TUKEY, 0, FCE_NORM_DC);
    FIR_CASE("fir_lp_rect", FCE_FIR_LOWPASS, 48000, 5000, 0, 51,
             FCE_WIN_RECTANGULAR, 0, FCE_NORM_DC);
    FIR_CASE("fir_hp_hamming", FCE_FIR_HIGHPASS, 48000, 5000, 0, 101,
             FCE_WIN_HAMMING, 0, FCE_NORM_NYQUIST);
    FIR_CASE("fir_hp_kaiser", FCE_FIR_HIGHPASS, 48000, 1000, 0, 201,
             FCE_WIN_KAISER, 5.0, FCE_NORM_NYQUIST);
    FIR_CASE("fir_bp_hamming", FCE_FIR_BANDPASS, 48000, 3000, 6000, 121,
             FCE_WIN_HAMMING, 0, FCE_NORM_PASSBAND_PEAK);
    FIR_CASE("fir_bs_kaiser", FCE_FIR_BANDSTOP, 48000, 3000, 6000, 161,
             FCE_WIN_KAISER, 6.2, FCE_NORM_DC);
    FIR_CASE("fir_hilbert", FCE_FIR_HILBERT, 48000, 0, 0, 65,
             FCE_WIN_HAMMING, 0, FCE_NORM_PASSBAND_PEAK);
    FIR_CASE("fir_diff", FCE_FIR_DIFFERENTIATOR, 48000, 0, 0, 64,
             FCE_WIN_HAMMING, 0, FCE_NORM_NYQUIST);
    {
        fce_spec_defaults(&sp);
        sp.kind = FCE_KIND_FIR;
        sp.fir_type = FCE_FIR_LOWPASS;
        sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 0;
        sp.window = FCE_WIN_KAISER;
        sp.stopband_atten_db = 80; sp.transition_hz = 1000;
        sp.precision = FCE_PRECISION_FLOAT64;
        dump("fir_lp_kaiser_auto", &sp);
    }

    /* ---------------- IIR ---------------- */
#define IIR_CASE(id_, fam_, type_, fs_, fc1_, fc2_, order_, rp_, rs_, norm_) \
    do { \
        fce_spec_defaults(&sp); \
        sp.kind = FCE_KIND_IIR; \
        sp.iir_family = (fam_); sp.iir_type = (type_); \
        sp.fs = (fs_); sp.fc1 = (fc1_); sp.fc2 = (fc2_); \
        sp.order = (order_); \
        sp.passband_ripple_db = (rp_); sp.stopband_atten_db = (rs_); \
        sp.precision = FCE_PRECISION_FLOAT64; \
        dump(id_, &sp); \
    } while (0)

    IIR_CASE("iir_butter_lp2", FCE_IIR_BUTTERWORTH, FCE_IIR_LOWPASS,
             48000, 5000, 0, 2, 0, 0, 0);
    IIR_CASE("iir_butter_lp4", FCE_IIR_BUTTERWORTH, FCE_IIR_LOWPASS,
             48000, 5000, 0, 4, 0, 0, 0);
    IIR_CASE("iir_butter_lp7", FCE_IIR_BUTTERWORTH, FCE_IIR_LOWPASS,
             10000, 1500, 0, 7, 0, 0, 0);
    IIR_CASE("iir_butter_hp3", FCE_IIR_BUTTERWORTH, FCE_IIR_HIGHPASS,
             48000, 1000, 0, 3, 0, 0, 0);
    IIR_CASE("iir_butter_hp8", FCE_IIR_BUTTERWORTH, FCE_IIR_HIGHPASS,
             44100, 5000, 0, 8, 0, 0, 0);
    IIR_CASE("iir_butter_bp4", FCE_IIR_BUTTERWORTH, FCE_IIR_BANDPASS,
             48000, 2000, 4000, 4, 0, 0, 0);
    IIR_CASE("iir_butter_bs6", FCE_IIR_BUTTERWORTH, FCE_IIR_BANDSTOP,
             48000, 3000, 6000, 6, 0, 0, 0);
    IIR_CASE("iir_cheby1_lp5", FCE_IIR_CHEBYSHEV1, FCE_IIR_LOWPASS,
             48000, 5000, 0, 5, 1.0, 0, 0);
    IIR_CASE("iir_cheby1_hp6", FCE_IIR_CHEBYSHEV1, FCE_IIR_HIGHPASS,
             48000, 2000, 0, 6, 0.5, 0, 0);
    IIR_CASE("iir_cheby1_bp3", FCE_IIR_CHEBYSHEV1, FCE_IIR_BANDPASS,
             48000, 1500, 3500, 3, 0.1, 0, 0);
    IIR_CASE("iir_cheby1_bs4", FCE_IIR_CHEBYSHEV1, FCE_IIR_BANDSTOP,
             48000, 2500, 4500, 4, 1.0, 0, 0);
    IIR_CASE("iir_cheby2_lp4", FCE_IIR_CHEBYSHEV2, FCE_IIR_LOWPASS,
             48000, 5000, 0, 4, 0, 60, 0);
    IIR_CASE("iir_cheby2_hp7", FCE_IIR_CHEBYSHEV2, FCE_IIR_HIGHPASS,
             48000, 1000, 0, 7, 0, 50, 0);
    IIR_CASE("iir_cheby2_bp6", FCE_IIR_CHEBYSHEV2, FCE_IIR_BANDPASS,
             48000, 3000, 6000, 6, 0, 60, 0);
    IIR_CASE("iir_cheby2_bs5", FCE_IIR_CHEBYSHEV2, FCE_IIR_BANDSTOP,
             48000, 3000, 6000, 5, 0, 40, 0);
    IIR_CASE("iir_ellip_lp4", FCE_IIR_ELLIPTIC, FCE_IIR_LOWPASS,
             48000, 5000, 0, 4, 0.5, 60, 0);
    IIR_CASE("iir_ellip_hp5", FCE_IIR_ELLIPTIC, FCE_IIR_HIGHPASS,
             48000, 3000, 0, 5, 0.1, 80, 0);
    IIR_CASE("iir_ellip_bp6", FCE_IIR_ELLIPTIC, FCE_IIR_BANDPASS,
             48000, 2000, 4000, 6, 1.0, 50, 0);
    IIR_CASE("iir_ellip_bs5", FCE_IIR_ELLIPTIC, FCE_IIR_BANDSTOP,
             48000, 3000, 6000, 5, 0.5, 60, 0);
    IIR_CASE("iir_bessel_lp3", FCE_IIR_BESSEL, FCE_IIR_LOWPASS,
             48000, 5000, 0, 3, 0, 0, 0);
    IIR_CASE("iir_bessel_lp6", FCE_IIR_BESSEL, FCE_IIR_LOWPASS,
             48000, 2000, 0, 6, 0, 0, 0);
    IIR_CASE("iir_bessel_hp4", FCE_IIR_BESSEL, FCE_IIR_HIGHPASS,
             48000, 3000, 0, 4, 0, 0, 0);
    IIR_CASE("iir_bessel_bp4", FCE_IIR_BESSEL, FCE_IIR_BANDPASS,
             48000, 1500, 3500, 4, 0, 0, 0);
    IIR_CASE("iir_bessel_bs4", FCE_IIR_BESSEL, FCE_IIR_BANDSTOP,
             48000, 2500, 4500, 4, 0, 0, 0);

    /* auto-order cases */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 10000; sp.fc1 = 1000; sp.order = 0;
    sp.edge1_hz = 2000; sp.passband_ripple_db = 3.0;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_butter_auto_lp", &sp);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_CHEBYSHEV1;
    sp.iir_type = FCE_IIR_HIGHPASS;
    sp.fs = 10000; sp.fc1 = 2000; sp.order = 0;
    sp.edge1_hz = 1000; sp.passband_ripple_db = 1.0;
    sp.stopband_atten_db = 50;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_cheby1_auto_hp", &sp);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_CHEBYSHEV2;
    sp.iir_type = FCE_IIR_BANDPASS;
    sp.fs = 48000; sp.fc1 = 2000; sp.fc2 = 4000; sp.order = 0;
    sp.edge1_hz = 1500; sp.edge2_hz = 5000;
    sp.passband_ripple_db = 1.0; sp.stopband_atten_db = 40;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_cheby2_auto_bp", &sp);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_BANDSTOP;
    sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 5000; sp.order = 0;
    sp.edge1_hz = 2000; sp.edge2_hz = 7000;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_ellip_auto_bs", &sp);

    /* fixed point cases */
    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_FIR;
    sp.fir_type = FCE_FIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
    sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("fir_lp_q15", &sp);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_BUTTERWORTH;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.qformat = FCE_QFORMAT_Q31;
    sp.scale_strategy = FCE_SCALE_SECTION_WISE;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_butter_lp4_q31_secwise", &sp);

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = FCE_IIR_ELLIPTIC;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.passband_ripple_db = 0.5; sp.stopband_atten_db = 60;
    sp.qformat = FCE_QFORMAT_Q15;
    sp.precision = FCE_PRECISION_FLOAT64;
    dump("iir_ellip_lp4_q15", &sp);

    return 0;
}
