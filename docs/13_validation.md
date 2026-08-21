# فصل ۱۳ — اعتبارسنجی

## مفهوم

هیچ ضرایبی بدون اثبات تحویل داده نمی‌شود. پس از تولید، کتابخانه به‌صورت خودکار (و قابل تنظیم) این موارد را بررسی می‌کند:

```text
پایداری        → شعاع قطب‌ها
پاسخ فرکانسی   → بهرهٔ DC/Nyquist، ripple باند عبور، تضعیف باند قطع، cutoff
انطباق با spec → تضعیف در لبهٔ باند قطع (حالت خودکار)
خطای کمی‌سازی  → float در برابر quantized (ضرایب + پاسخ + پایداری)
```

## شهود

اعتبارسنجی یعنی «آیا این ضرایب همان کاری را می‌کنند که کاربر خواسته؟». پاسخ فرکانسی روی یک شبکهٔ فرکانسی ارزیابی می‌شود و اعداد اندازه‌گیری‌شده در `fce_result_t` و در گزارش (فصل ۱۴) می‌آیند.

## فرمول‌ها

ارزیابی یک سکشن در فرکانس نرمال‌شدهٔ `ω`:

```
H(e^jω) = (b0 + b1 e^−jω + b2 e^−j2ω) / (1 + a1 e^−jω + a2 e^−j2ω)
|H|کل = Π|H_section|
20 log10|H|  ← برای dB
```

## API — اسکن پاسخ

```c
/* callback برای دریافت نقاط */
static bool my_cb(void* ctx, const fce_response_point_t* pt)
{
    printf("%.1f Hz  %.3f dB  %.1f deg  GD=%.3f\n",
           pt->f_hz, pt->mag_db, pt->phase_deg, pt->group_delay);
    return true;   /* false = توقف اسکن */
}

fce_response_sos(result.sos_f64, result.num_sections, result.fs,
                 512, 0.0, 24000.0, my_cb, NULL);

/* برای FIR: */
fce_response_fir(result.h_f64, result.num_taps, result.fs,
                 512, 0.0, 24000.0, my_cb, NULL);
```

## خروجی‌های اندازه‌گیری‌شده

| فیلد | معنی |
|---|---|
| `dc_gain_db` | بهره در DC |
| `nyquist_gain_db` | بهره در Nyquist |
| `passband_ripple_measured_db` | پهنای نوسان باند عبور |
| `stopband_atten_measured_db` | تضعیف باند قطع (dB مثبت) |
| `cutoff_measured_hz` | فرکانس قطع اندازه‌گیری‌شده (−3dB یا −gpass) |
| `quant_response_max_error_db` | اختلاف پاسخ float/quantized |

## خطای ضرایب (ابزار عمومی)

```c
double max_abs, rms, max_rel;
fce_coeff_error(ref_coeffs, test_coeffs, n, &max_abs, &rms, &max_rel);
```

## کنترل subset اعتبارسنجی

```c
sp.validate = FCE_VALIDATE_STABILITY | FCE_VALIDATE_SPEC;
/* یا */
sp.validate = FCE_VALIDATE_NONE;   /* فقط طراحی، بدون بررسی */
```

## مثال — گزارش یک طراحی

```
Passband ripple      : 0.0002 dB
Stopband attenuation : 80.4 dB   (مطلوب: 80 dB ✓)
Measured cutoff      : 5000.3 Hz
```

## استفادهٔ عملی

- در حالت پیش‌فرض همه‌چیز فعال است؛ هزینهٔ آن یک اسکن ۵۱۲ نقطه‌ای (~چند میلی‌ثانیه روی H7) است.
- در تولید انبوه روی MCU، می‌توانید `FCE_VALIDATE_NONE` بگذارید و فقط یک بار در تست‌ها اعتبارسنجی کامل انجام دهید.

[فصل بعد: Export](14_export.md)
