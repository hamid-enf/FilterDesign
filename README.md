# FilterCoeff — تولیدکننده ضرایب فیلتر دیجیتال

> **«مشخصات فیلتر را به من بده؛ ضرایب دقیق و آماده استفاده را تحویل بگیر.»**

یک کتابخانهٔ حرفه‌ای، دقیق و Embedded-Friendly به زبان **C99** برای **محاسبه و تولید ضرایب فیلترهای دیجیتال** (FIR و IIR). این کتابخانه فقط *تولیدکنندهٔ ضریب* است — فیلتر اجرا نمی‌کند. خروجی مستقیماً قابل استفاده در کتابخانه‌های FIR/IIR موجود شماست.

[![CI](https://github.com/hamid-enf/FilterDesign/actions/workflows/ci.yml/badge.svg)](https://github.com/hamid-enf/FilterDesign/actions/workflows/ci.yml)

---

## چرا FilterCoeff؟

| نیاز شما | راه‌حل FilterCoeff |
|---|---|
| ضرایب FIR دقیق (Windowed-Sinc + Kaiser) | `fce_generate()` + مشخصات ساده |
| ضرایب IIR (Butterworth/Chebyshev I-II/Elliptic/Bessel) | خروجی SOS استاندارد |
| خروجی float32 / float64 | محاسبات داخلی float64، گرد کردن نهایی |
| خروجی Q15 / Q31 | مقیاس‌گذاری + تحلیل خطا + بررسی پایداری |
| اعتبارسنجی ضرایب | پاسخ فرکانسی، پایداری، خطای کمی‌سازی |
| خروجی آمادهٔ کدنویسی | C / CSV / JSON / گزارش Markdown |
| اجرا روی STM32H7 و ESP32-C2 | بدون وابستگی به HAL/IDF/FreeRTOS |

## شروع سریع

```c
#include "filtercoeff.h"

FilterCoeffSpec spec;            /* یا fce_spec_t */
FilterCoeffResult result;
static uint8_t mem[8192];
FilterCoeffWorkspace ws = { mem, sizeof(mem) };

fce_spec_defaults(&spec);
spec.kind    = FCE_KIND_FIR;
spec.fir_type = FCE_FIR_LOWPASS;
spec.fs      = 48000;            /* نرخ نمونه‌برداری */
spec.fc1     = 5000;             /* فرکانس قطع */
spec.num_taps = 101;             /* تعداد ضریب */
spec.window  = FCE_WIN_KAISER;
spec.stopband_atten_db = 80;     /* برای محاسبهٔ خودکار beta */

if (fce_generate(&spec, &result, &ws) == FCE_OK)
{
    /* result.h_f32 / result.h_f64 -> مستقیم به کتابخانه FIR شما */
    FilterLab_FIR_SetCoefficients(result.h_f32, result.num_taps);
}
```

## نصب و ساخت

```sh
make            # ساخت libfiltercoeff.a
make test       # اجرای ۵۹۰۰+ تست واحد
make ref        # مقایسه با SciPy (نیازمند python3 + scipy)
make examples   # ساخت ۲۰ مثال
make bench      # بنچمارک زمان طراحی
```

فقط فایل‌های `include/` و `src/` را به پروژهٔ خود اضافه کنید؛ هیچ وابستگی خارجی جز `<stdint.h> <stddef.h> <stdbool.h> <math.h>` نیست.

## ساختار مخزن

```text
FilterCoeff/
├── include/          ← فایل‌های سرآیند عمومی (filtercoeff.h)
├── src/              ← پیاده‌سازی (C99، بدون تخصیص پویا)
├── examples/         ← ۲۰ مثال کامل
├── tests/            ← تست‌های واحد (بیش از ۵۹۰۰ بررسی)
├── bench/            ← بنچمارک زمان طراحی
├── tools/reference/  ← ابزار مقایسه با SciPy
├── docs/             ← مستندات کامل فارسی (۱۸ فصل)
└── .github/workflows/← CI
```

## مستندات

| فصل | موضوع |
|---|---|
| [فصل ۱](docs/01_introduction.md) | معرفی چارچوب |
| [فصل ۲](docs/02_specification.md) | تعریف Specification |
| [فصل ۳](docs/03_fir_generation.md) | تولید ضرایب FIR |
| [فصل ۴](docs/04_fir_windows.md) | پنجره‌های FIR |
| [فصل ۵](docs/05_iir_generation.md) | تولید ضرایب IIR |
| [فصل ۶](docs/06_iir_families.md) | خانواده‌های IIR |
| [فصل ۷](docs/07_sos_biquad.md) | ضرایب SOS/Biquad |
| [فصل ۸](docs/08_normalization.md) | نرمال‌سازی |
| [فصل ۹](docs/09_float_coefficients.md) | ضرایب اعشاری |
| [فصل ۱۰](docs/10_fixed_point.md) | Fixed-Point / Q15 / Q31 |
| [فصل ۱۱](docs/11_quantization_error.md) | خطای کمی‌سازی |
| [فصل ۱۲](docs/12_numerical_stability.md) | پایداری عددی |
| [فصل ۱۳](docs/13_validation.md) | اعتبارسنجی |
| [فصل ۱۴](docs/14_export.md) | خروجی (Export) |
| [فصل ۱۵](docs/15_stm32h7.md) | STM32H7 |
| [فصل ۱۶](docs/16_esp32c2.md) | ESP32-C2 |
| [فصل ۱۷](docs/17_examples.md) | مثال‌ها |
| [فصل ۱۸](docs/18_troubleshooting.md) | عیب‌یابی |
| [پژوهش فنی](docs/ALGORITHM_REFERENCES.md) | منابع و فرمول‌های الگوریتم‌ها |
| [اعتبارسنجی مرجع](docs/reference_validation.md) | روش مقایسه با SciPy |
| CI | فایل `ci/ci.yml` — برای فعال‌سازی به `.github/workflows/` کپی کنید |
| Adapter | `include/filtercoeff_adapter.h` — تبدیل به فرمت scipy/MATLAB/CMSIS |

## ویژگی‌ها

- **FIR:** Lowpass / Highpass / Bandpass / Bandstop / Hilbert / Differentiator
  با پنجره‌های Rectangular, Hann, Hamming, Blackman, Kaiser, Blackman-Harris, Bartlett, Tukey
- **IIR:** Butterworth / Chebyshev I / Chebyshev II / Elliptic (Cauer) / Bessel
  در هر چهار نوع پاسخ (LP/HP/BP/BS) با خروجی SOS
- **محاسبهٔ خودکار:** تعداد taps برای Kaiser (فرمول Kaiser) و مرتبهٔ IIR (معادل `buttord`/`cheb1ord`/`cheb2ord`/`ellipord`)
- **Fixed-point:** Q15/Q31 با سه استراتژی مقیاس‌گذاری (سراسری / بخش‌به‌بخش / ضریب‌به‌ضریب)، گزارش scale، خطا و پایداری پس از کمی‌سازی
- **اعتبارسنجی:** پایداری (شعاع قطب‌ها)، پاسخ فرکانسی، ripple، تضعیف، cutoff، خطای ضرایب
- **بدون تخصیص پویا:** حافظهٔ کاری را کاربر تأمین می‌کند (`fce_workspace_t`)
- **بدون black box:** ضرایب میانی (پنجره، پاسخ ایده‌آل، قطب‌ها/صفرها، مقیاس‌ها) در دسترس‌اند

## امضای قرارداد ضرایب

IIR (SOS) — هر سکشن به صورت `{ b0, b1, b2, a1, a2 }` با `a0 = 1` ضمنی:

```
y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
```

## اعتبارسنجی در برابر مرجع

پس از هر تغییر، تمام ضرایب با SciPy مقایسه می‌شوند:

```sh
make ref
# FilterCoeff vs SciPy: 105 checks, 0 failures
```

همچنین ۶۰ طراحی تصادفی (مرتبهٔ ۲ تا ۱۰، هر ۵ خانواده، هر ۴ نوع پاسخ) با حداکثر خطای ضرایب ~1e-13 بررسی شده است.

## مجوز

آزاد برای استفاده در پروژه‌های شخصی و تجاری. (هیچ‌گونه وابستگی به کتابخانه‌های DSP خارجی ندارد.)
