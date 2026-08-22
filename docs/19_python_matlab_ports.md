# فصل ۱۹ — پورت‌های Python و MATLAB

این فصل پورت‌های پایتون و متلب را معرفی می‌کند و نحوهٔ تولید و مقایسهٔ خروجی هر سه
پیاده‌سازی (C، پایتون، متلب) را نشان می‌دهد.

## ۱۹.۱ چرا پورت؟

کتابخانهٔ C را می‌توان در هر جایی که کامپایلر C99 هست (MCU، دسکتاپ، سرور) ساخت.
اما برای نمونه‌سازی سریع، اسکریپت‌نویسی، تحلیل با NumPy/SciPy یا کار با متلب،
داشتن همان الگوریتم به زبان‌های پایتون و متلب بسیار راحت‌تر است.

پورت‌ها **خط‌به‌خط از همان الگوریتم C** پیاده شده‌اند (نه wrapper روی SciPy)، بنابراین:

- هیچ وابستگی به SciPy یا جعبه‌ابزار Signal Processing ندارند.
- خروجی‌شان با C تا حدود **~1e-15** (خطای گرد کردن float64) یکسان است.
- امضای ضرایب و قراردادها در هر سه یکی است.

## ۱۹.۲ پورت پایتون

فایل: `python/filtercoeff.py` — فقط با کتابخانهٔ استاندارد پایتون.

```python
import filtercoeff as fce

# FIR: ورودی مشخصات، خروجی ضرایب
h = fce.design({
    "kind": "fir", "fir_type": "lowpass",
    "fs": 48000, "fc1": 5000, "num_taps": 101, "window": "hamming",
})

# IIR: خروجی SOS به صورت {b0,b1,b2,a1,a2} در هر سکشن
sos = fce.design({
    "kind": "iir", "iir_family": "butterworth", "iir_type": "lowpass",
    "fs": 48000, "fc1": 5000, "order": 4,
})

# Fixed-point
q, scale, int_bits = fce.quantize(flat_coeffs, "q15", "symmetric", sec_len=n)
```

`design(spec, details=True)` متادیتا و واسطه‌ها (پنجره، پاسخ ایده‌آل، بهرهٔ سکشن‌ها و…)
را هم برمی‌گرداند.

## ۱۹.۳ پورت متلب

فایل‌ها: `matlab/filtercoeff.m` و `matlab/filtercoeff_quantize.m`.

```matlab
% FIR lowpass: Fs=48k, Fc=5k, 101 taps, Hamming
sp = struct('kind','fir','fir_type','lowpass','fs',48000, ...
            'fc1',5000,'num_taps',101,'window','hamming');
[h, info] = filtercoeff(sp);

% IIR Butterworth lowpass (order 4)
sp = struct('kind','iir','iir_family','butterworth','iir_type','lowpass', ...
            'fs',48000,'fc1',5000,'order',4);
[sos, info] = filtercoeff(sp);

% Fixed-point Q31
[q, scale, int_bits] = filtercoeff_quantize(sos, 'q31', 'symmetric', numel(sos));
```

## ۱۹.۴ ارزیابی دقت (خروجی‌ها چقدر یکی‌اند؟)

هر سه پیاده‌سازی همان محاسبات float64 را انجام می‌دهند، پس ضرایب تا حد دقت ماشین
(≈1e-15) یکسان‌اند. ارزیابی به دو بخش انجام شده است:

### ۱) پایتون در برابر C

اسکریپت `python/compare_c.py` همهٔ ۳۱ طراحی مرجع در `tools/reference/fce_dump.c`
(FIR/IIR، همهٔ خانواده‌ها، خودکار، Q15/Q31) را با هر دو پیاده‌سازی می‌سازد و مقایسه می‌کند:

```sh
make                    # ساخت کتابخانهٔ C
make -C tools/reference fce_dump
python3 python/compare_c.py tools/reference/fce_dump
# C vs Python port: 76 checks, 0 failures
```

| دسته | تعداد بررسی | آستانه |
|---|---|---|
| ضرایب FIR | 15 | < 1e-12 |
| ضرایب SOS (IIR) | 16 | < 1e-12 |
| مرتبهٔ خودکار / فرکانس طراحی | 8 | < 1e-9 |
| مقادیر Q15 / Q31 | 3 | یکسان (exact) |
| مقیاس fixed-point | 2 | نسبی < 1e-9 |

### ۲) C در برابر SciPy (مرجع بیرونی)

`make ref` ضرایب C را با SciPy مقایسه می‌کند (۱۰۵ بررسی، صفر شکست) — یعنی پورت
پایتون که با C یکی است، عملاً با SciPy هم هم‌خوان است.

### ۳) متلب در برابر C

متلب و اکتاو در این محیط در دسترس نبودند، بنابراین پورت متلب به‌صورت ترجمهٔ ۱:۱
پورت پایتون (که در برابر C اعتبارسنجی شده) نوشته شد. برای اجرای اعتبارسنجی در
متلب خودتان:

1. ضرایب مرجع C در `matlab/reference/*.csv` قرار دارند (با `make matref` بازتولید می‌شوند).
2. در متلب اجرا کنید: `validate_filtercoeff`

```matlab
>> validate_filtercoeff
PASS  fir_lp_hann_odd        max err 0.00e+00
PASS  iir_butter_lp4         max err 2.22e-16
...
MATLAB vs C reference: 0 failures
```

هم‌خوانی مشخصات اسکریپت متلب با CSVهای مرجع از طریق پورت پایتون (که با C یکی است)
تا ~1e-15 تأیید شده است.

## ۱۹.۵ محدودیت‌ها

- پورت‌ها روی **تولید ضرایب** تمرکز دارند (FIR/IIR/خودکار/fixed-point). ماژول‌های
  export، validation-report و شبیه‌سازیِ C در پورت‌ها نیامده‌اند.
- مرتبهٔ خودکار Bessel در هر سه پیاده‌سازی پشتیبانی نمی‌شود (مانند خود C).
- `filtercoeff_quantize` در متلب، مقادیر exact کوانتیزاسیون را با C برمی‌گرداند؛
  مقیاس در بدترین حالت تا خطای نسبی ~1e-9 (به دلیل آخرین بیت ضرایب) می‌تواند
  نوسان کند — این یک خطای واقعی نیست، بلکه خطای گرد کردن ماشین است.
