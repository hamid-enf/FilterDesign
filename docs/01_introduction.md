# فصل ۱ — معرفی چارچوب

## مفهوم (Concept)

**FilterCoeff** یک کتابخانهٔ تولیدکنندهٔ ضرایب فیلتر دیجیتال است. ورودی آن یک *مشخصات فیلتر* (Filter Specification) و خروجی آن مجموعه‌ای از ضرایب *معتبر، دقیق و آمادهٔ استفاده* است — در قالب float32، float64، Q15 یا Q31 و در قالب آرایهٔ FIR یا SOS برای IIR.

این کتابخانه **اجراکنندهٔ فیلتر نیست**. پردازش realtime، DMA، ADC و streaming در حیطهٔ این پروژه نیستند. هر ابزار شبیه‌سازی یا تحلیل پاسخ که در این مخزن می‌بینید، صرفاً برای **اثبات صحت ضرایب** وجود دارد.

## شهود (Intuition)

کاربر می‌گوید:

> «یک فیلتر Low-pass با نرخ نمونه‌برداری ۴۸ کیلوهرتز، فرکانس قطع ۵ کیلوهرتز و ۸۰ دسی‌بل تضعیف می‌خواهم.»

و کتابخانه برمی‌گرداند:

```c
static const float32_t coeffs[243] = {
    -4.2799849e-06, ... /* ۲۴۳ ضریب دقیق */
};
```

بین این دو فقط یک فراخوانی API فاصله است:

```c
FilterCoeff_Generate(&spec, &result);
```

## معماری

```text
User Specification
        ↓
fce_generate()
   ├── اعتبارسنجی مشخصات
   ├── محاسبهٔ خودکار (taps / order)
   ├── طراحی FIR / IIR  ← محاسبات float64
   ├── کمی‌سازی Q15/Q31 (اختیاری)
   ├── اعتبارسنجی (پایداری، پاسخ، خطا)
   └── نتیجه + ضرایب
        ↓
float32 / float64 / Q15 / Q31  ←  آماده برای کتابخانهٔ FIR/IIR شما
```

## جعبه‌ابزار (ماژول‌ها)

| ماژول | فایل | نقش |
|---|---|---|
| Core | `fce_generate.c` | هماهنگی کل pipeline |
| Math | `fce_math.c` | توابع بیضوی، I0، ریشه‌یابی، پنجره‌ها |
| FIR | `fce_fir.c` | طراحی پنجره‌ای-سینک |
| IIR | `fce_iir.c` | پروتوتایپ آنالوگ → تبدیل → SOS |
| Quantization | `fce_quant.c` | Q15/Q31 با مقیاس‌گذاری |
| Validation | `fce_validate.c` | پاسخ فرکانسی، پایداری، خطا |
| Simulation | `fce_sim.c` | سیگنال‌های تست (فقط اعتبارسنجی) |
| Export | `fce_export.c` | C / CSV / JSON / Report |

## API

ساده‌ترین استفاده:

```c
#include "filtercoeff.h"

FilterCoeffSpec    spec;
FilterCoeffResult  result;
static uint8_t     mem[65536];
FilterCoeffWorkspace ws = { mem, sizeof(mem) };

fce_spec_defaults(&spec);
spec.kind = FCE_KIND_FIR;
spec.fir_type = FCE_FIR_LOWPASS;
spec.fs = 48000;
spec.fc1 = 5000;
spec.num_taps = 101;
spec.window = FCE_WIN_HAMMING;

if (FilterCoeff_Generate(&spec, &result, &ws) == FCE_OK)
{
    /* ضرایب آماده‌اند */
}
```

## قراردادهای مهم

1. **حافظه:** هیچ تخصیص پویایی وجود ندارد؛ آرایه‌های خروجی داخل `ws.data` قرار می‌گیرند و تا وقتی به ضرایب نیاز دارید باید زنده بماند.
2. **دقت:** تمام محاسبات داخلی `double` (float64) است؛ float32 فقط گرد کردن نهایی است.
3. **قرارداد علامت (SOS):** `y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] − a1 y[n-1] − a2 y[n-2]` — آرایهٔ هر سکشن `{ b0, b1, b2, a1, a2 }`.
4. **کد خطا:** `result.status` همیشه پر می‌شود؛ هشدارهای غیرمرگبار در `result.flags` (بیت‌ماسک) می‌آیند.

## مثال عددی کوچک

فیلتر میانگین‌گیر ۳ ضریبی در `fc = fs/4` (پنجرهٔ مستطیلی):

```
h = [0.28005, 0.43990, 0.28005]   (نرمال‌شده به بهرهٔ DC = 1)
```

## استفادهٔ عملی

- ضرایب را یک بار در زمان راه‌اندازی (boot) محاسبه کنید و به کتابخانهٔ runtime خود بدهید؛
- یا در Host (PC/Linux) محاسبه و به‌صورت `static const` در کد هدف قرار دهید (فصل ۱۴ — Export).

## مطالعهٔ بعدی

[فصل ۲: تعریف Specification](02_specification.md)
