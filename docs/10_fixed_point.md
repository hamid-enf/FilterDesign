# فصل ۱۰ — Fixed-Point / Q15 / Q31

## مفهوم

تبدیل fixed-point در این کتابخانه **هرگز یک cast ساده نیست**:

```c
int16_t q = (int16_t)float_value;   /* ❌ هرگز */
```

بلکه یک pipeline کامل است:

```text
ضرایب float
    ↓ تعیین مقیاس (Scale)
    ↓ تحلیل محدوده (Range Analysis)
    ↓ کمی‌سازی (round-to-nearest)
    ↓ بررسی اشباع / سرریز (Saturation / Overflow)
    ↓ محاسبهٔ خطای کمی‌سازی
    ↓ اعتبارسنجی (پاسخ، پایداری)
```

## شهود

فرمت Qm.F یعنی F بیت کسری. برای ضریب `c`:

```
q = round(c · 2^F)          (با محدودیت به ±(2^F − 1))
c̃ = q / 2^F                 (بازسازی)
خطای هر ضریب ≤ 0.5 / 2^F
```

اگر `|c| > 1`، به بیت‌های صحیح (m) نیاز داریم: `q = round(c · 2^F)` در Q15 دیگر جا نمی‌شود. راه‌حل: مقیاس سراسری کوچک‌تر — یعنی `scale = (2^F − 1) / max|c|` و `q = round(c · scale)`. این مقیاس در `result.scale` گزارش می‌شود.

## فرمت‌ها و پارامترهای گزارش‌شده

| فیلد | معنی |
|---|---|
| `result.scale` | مقیاس (ضریب) — `q = round(c · scale)` |
| `result.q_int_bits` | تعداد بیت‌های صحیح مؤثر (Qm.F) |
| `result.q_max_abs_error` | بیشینهٔ خطای مطلق ضریب |
| `result.q_rms_error` | خطای RMS |
| `result.q_max_rel_error` | بیشینهٔ خطای نسبی |
| `result.section_scales` | مقیاس هر سکشن (در SECTION_WISE) |
| `result.coeff_scales` | مقیاس هر ضریب (در COEFFICIENT_WISE) |
| `result.quant_response_max_error_db` | بیشینهٔ اختلاف پاسخ (dB) |
| `result.quant_max_pole_radius` | شعاع قطب پس از کمی‌سازی |

## استراتژی‌های مقیاس‌گذاری

### ۱. سراسری (Symmetric) — پیش‌فرض

```c
scale = (2^F − 1) / max|all coefficients|
```

یک مقیاس برای همه. ساده و قابل حمل؛ اما اگر یک ضریب بزرگ باشد (مثل a1 در فیلتر باریک)، بقیه ضریب‌ها دقت کمتری می‌گیرند.

### ۲. بخش‌به‌بخش (Section-wise) — توصیه‌شده برای IIR

```c
scale_s = (2^F − 1) / max|coefficients of section s|
```

هر سکشن SOS مقیاس خودش را دارد. دقت هر سکشن بیشینه است؛ runtime باید مقیاس‌ها را مدیریت کند (معمولاً با شیفت در فرم Direct-Form Transposed).

### ۳. ضریب‌به‌ضریب (Coefficient-wise)

هر ضریب مقیاس خودش را دارد و تقریباً همیشه به بیشینهٔ فرمت می‌رسد (storage-only؛ در عمل به‌ندرت قابل استفاده در runtime).

## مثال عددی — FIR Kaiser Q15

`fs=48k, fc=5k, N=101`:

```
max|c| = 0.20834
scale  = 32767 / 0.20834 = 157276.24
q[50]  = round(0.20834 × 157276.24) = 32767
q[0]   = round(1.639e-5 × 157276.24) = 3
خطای بیشینه = 3.09e-6   (< 1/scale = 6.36e-6)
```

## کد

```c
sp.qformat = FCE_QFORMAT_Q15;                  /* یا Q31 */
sp.scale_strategy = FCE_SCALE_SECTION_WISE;    /* برای IIR */

fce_generate(&sp, &result, &ws);

/* ضرایب آماده */
const int16_t* q15 = result.q15;               /* FIR: N ضریب، IIR: 5×S ضریب */
```

## خطاها و پرچم‌ها

| شرط | رفتار |
|---|---|
| سرریز واقعی (بعد از اشباع) | `FCE_FLAG_COEFFICIENT_OVERFLOW` + `FCE_ERR_OVERFLOW` |
| فیلتر کمی‌شده ناپایدار | `FCE_FLAG_QUANTIZATION_UNSTABLE` + `FCE_ERR_QUANTIZATION` |
| خطای کمی‌سازی بالا | `FCE_FLAG_QUANTIZATION_WARNING` |
| حاشیهٔ پایداری کم | `FCE_FLAG_SPEC_MARGINAL` |

## استفادهٔ عملی — CMSIS-DSP

```c
/* Q15 FIR */
arm_fir_q15(&S, result.q15, num_taps, ...);    /* ضرایب با scale = 2^15 */

/* Q15 IIR: با مقیاس بخش‌به‌بخش، هر سکشن را جداگانه با شیفت اعمال کنید */
```

[فصل بعد: خطای کمی‌سازی](11_quantization_error.md)
