# فصل ۱۷ — مثال‌ها

## مفهوم

۲۰ مثال کامل در `examples/` وجود دارد. همه با `make examples` ساخته می‌شوند و خروجی‌شان گزارش + ضرایب است.

## فهرست مثال‌ها

| # | فایل | موضوع |
|---|---|---|
| ۱ | `example_01_fir_lowpass.c` | FIR پایین‌گذر (Hamming) |
| ۲ | `example_02_fir_highpass.c` | FIR بالاگذر |
| ۳ | `example_03_fir_bandpass.c` | FIR باندگذر |
| ۴ | `example_04_fir_bandstop.c` | FIR باندقطع |
| ۵ | `example_05_fir_kaiser.c` | Kaiser خودکار (taps + beta) |
| ۶ | `example_06_iir_butterworth.c` | Butterworth LP |
| ۷ | `example_07_iir_chebyshev.c` | Chebyshev I و II |
| ۸ | `example_08_iir_elliptic.c` | Elliptic |
| ۹ | `example_09_iir_bessel.c` | Bessel |
| ۱۰ | `example_10_float32.c` | خروجی float32 |
| ۱۱ | `example_11_float64.c` | خروجی float64 |
| ۱۲ | `example_12_q15.c` | تبدیل Q15 (FIR) |
| ۱۳ | `example_13_q31.c` | تبدیل Q31 (IIR, section-wise) |
| ۱۴ | `example_14_sos_export.c` | خروجی SOS برای کتابخانهٔ خارجی |
| ۱۵ | `example_15_stm32h7.c` | الگوی یکپارچه‌سازی STM32H7 |
| ۱۶ | `example_16_esp32c2.c` | الگوی یکپارچه‌سازی ESP32-C2 |
| ۱۷ | `example_17_quantization_validation.c` | اعتبارسنجی کمی‌سازی |
| ۱۸ | `example_18_stability_validation.c` | اعتبارسنجی پایداری |
| ۹۸ | `example_98_auto_order.c` | مرتبهٔ خودکار (مثال mission) |
| ۹۹ | `example_99_simple_api.c` | Simple API یک‌خطی |

## اجرا

```sh
make examples
./examples/example_05_fir_kaiser
./examples/example_98_auto_order
```

## مثال ۵ — Kaiser خودکار (خروجی)

```
Taps      : 243
Kaiser β  : 7.8573
Stopband  : 80.1 dB   (مطلوب 80 ✓)
```

## مثال ۹۸ — مرتبهٔ خودکار (خروجی)

مشخصات mission: LP 1kHz، fs=10k، تضعیف ≥60dB در 2kHz

```
Order     : 9  (مطابق scipy.buttord)
Sections  : 5
Validation: PASS
```

## مثال ۹۹ — Simple API

```c
FilterCoeffSpec spec;
FilterCoeffResult result;
static uint8_t mem[65536];
FilterCoeffWorkspace ws = { mem, sizeof(mem) };

fce_spec_defaults(&spec);
spec.kind = FCE_KIND_FIR;
spec.fir_type = FCE_FIR_LOWPASS;
spec.fs = 48000;
spec.fc1 = 5000;
spec.num_taps = 101;
spec.window = FCE_WIN_KAISER;
spec.stopband_atten_db = 80;
spec.precision = FCE_PRECISION_FLOAT32;

FilterCoeff_Generate(&spec, &result, &ws);
/* result.h_f32 → FilterLab_FIR_SetCoefficients(...) */
```

## استفادهٔ عملی

- مثال‌ها را به‌عنوان الگو کپی کنید؛ همه از `common.h` استفاده می‌کنند که فقط یک helper ساده است.
- برای پروژهٔ خودتان، `fce_workspace_required()` را صدا بزنید و فقط همان اندازه حافظه بدهید.

[فصل بعد: عیب‌یابی](18_troubleshooting.md)
