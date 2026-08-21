# مستندات FilterCoeff

| # | فصل | موضوع |
|---|---|---|
| ۱ | [معرفی چارچوب](01_introduction.md) | هدف، معماری، API ساده |
| ۲ | [تعریف Specification](02_specification.md) | ساختار `fce_spec_t` و قواعد اعتبارسنجی |
| ۳ | [تولید ضرایب FIR](03_fir_generation.md) | Windowed-Sinc، انواع فیلتر، taps خودکار |
| ۴ | [پنجره‌های FIR](04_fir_windows.md) | ۸ پنجره، فرمول‌ها و مقایسه |
| ۵ | [تولید ضرایب IIR](05_iir_generation.md) | پروتوتایپ → تبدیل → دوخطی → SOS |
| ۶ | [خانواده‌های IIR](06_iir_families.md) | Butterworth/Chebyshev/Elliptic/Bessel |
| ۷ | [SOS/Biquad](07_sos_biquad.md) | قرارداد علامت، جفت‌کردن، مرتب‌سازی |
| ۸ | [نرمال‌سازی](08_normalization.md) | استراتژی‌های بهره |
| ۹ | [ضرایب اعشاری](09_float_coefficients.md) | float32/float64 |
| ۱۰ | [Fixed-Point](10_fixed_point.md) | Q15/Q31، سه استراتژی مقیاس |
| ۱۱ | [خطای کمی‌سازی](11_quantization_error.md) | معیارهای خطا و مثال‌ها |
| ۱۲ | [پایداری عددی](12_numerical_stability.md) | قطب‌ها، حاشیهٔ پایداری |
| ۱۳ | [اعتبارسنجی](13_validation.md) | پاسخ فرکانسی، ripple، تضعیف |
| ۱۴ | [خروجی (Export)](14_export.md) | C/CSV/JSON/Report |
| ۱۵ | [STM32H7](15_stm32h7.md) | راه‌اندازی و نکات |
| ۱۶ | [ESP32-C2](16_esp32c2.md) | راه‌اندازی و نکات |
| ۱۷ | [مثال‌ها](17_examples.md) | ۲۰ مثال کامل |
| ۱۸ | [عیب‌یابی](18_troubleshooting.md) | خطاها و راه‌حل‌ها |
| — | [پژوهش فنی](ALGORITHM_REFERENCES.md) | منابع و فرمول‌ها |
| — | [اعتبارسنجی مرجع](reference_validation.md) | روش مقایسه با SciPy |

## ساختار هر فصل

هر فصل از الگوی زیر پیروی می‌کند:

```text
مفهوم (Concept)
↓
شهود (Intuition)
↓
فرمول (Formula)
↓
مثال عددی (Numerical Example)
↓
API
↓
کد (Code)
↓
ضرایب تولیدشده (Generated Coefficients)
↓
استفادهٔ عملی (Practical Usage)
```
