# فصل ۹ — ضرایب اعشاری (float32 / float64)

## مفهوم

کتابخانه دو سطح دقت خروجی دارد:

```text
مشخصات ورودی
    ↓
محاسبات داخلی float64 (همیشه)
    ↓
گرد کردن نهایی به float32 (فقط اگر خواسته شود)
```

## شهود

اگر از ابتدا با float32 محاسبه کنید، خطاها در مراحل میانی (مثل توابع بیضوی یا تبدیل دوخطی) انباشته می‌شوند. قانون طلایی: **با بیشترین دقت ممکن محاسبه کن، در انتها گرد کن.**

## جزئیات

| مورد | float32 | float64 |
|---|---|---|
| محاسبات داخلی | float64 | float64 |
| ضرایب خروجی | `float` (نزدیک‌ترین گرد شدن) | `double` |
| خطای هر ضریب | ≤ 2^-24 نسبی (~6e-8) | ~1e-16 |
| حافظهٔ ضرایب | ۴ بایت × N | ۸ بایت × N |
| فیلد result | `h_f32` / `sos_f32` | `h_f64` / `sos_f64` |

توجه: `h_f64`/`sos_f64` **همیشه** موجود است — حتی وقتی `precision = FCE_PRECISION_FLOAT32` — چون محاسبات داخلی float64 است و نسخهٔ double در workspace نگهداری می‌شود.

## مثال عددی

FIR Kaiser، `fs=48k, fc=5k, N=101`:

```c
sp.precision = FCE_PRECISION_FLOAT32;
```

ضریب مرکزی:
- float64: `0.2083404271070948`
- float32: `0.20834042`

خطای گرد کردن: `~1.5e-9` (کاملاً در حد انتظار float32).

## API

```c
sp.precision = FCE_PRECISION_FLOAT32;   /* یا FCE_PRECISION_FLOAT64 */

const double* h  = result.h_f64;        /* همیشه */
const float*  hf = result.h_f32;        /* فقط با FLOAT32 */
```

## کد — تولید مستقیم آرایهٔ C

```c
char buf[32768];
fce_mem_writer_t mw;
fce_writer_t w;
fce_writer_mem_init(&w, &mw, buf, sizeof(buf));

fce_export_c_fir(&result, NULL, &w);    /* float32_t coeffs[] = {...}; */
```

خروجی:

```c
static const float32_t coeffs[101] = {
    1.6393232e-05f,
    1.8045962e-05f,
    ...
};
```

## استفادهٔ عملی

- **H7 / PC:** float64 رایگان است (FPU double روی Cortex-M7).
- **ESP32-C2:** بدون FPU سخت‌افزاری؛ float64 نرم‌افزاری است. اگر ضرایب را در Host تولید و به‌صورت ثابت embed می‌کنید، float32 کافی است (خطای 6e-8 در ضرایب معمولاً بی‌اثر است).
- قانون کلی: **اگر ضرایب روی target محاسبه می‌شوند، float64؛ اگر از Host می‌آیند، float32.**

[فصل بعد: Fixed-Point](10_fixed_point.md)
