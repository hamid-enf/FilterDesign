# فصل ۳ — تولید ضرایب FIR

## مفهوم

تولید FIR در این کتابخانه با روش **Windowed-Sinc** انجام می‌شود: پاسخ ضربهٔ ایده‌آل (سینک) محاسبه، با پنجره محدود، و سپس نرمال می‌شود.

```text
پاسخ ضربهٔ ایده‌آل (ideal h[n])
        ↓ ضرب در پنجره w[n]
پاسخ ضربهٔ پنجره‌شده
        ↓ نرمال‌سازی بهره
ضرایب نهایی
```

## شهود

فیلتر پایین‌گذر ایده‌آل پاسخ ضربه‌ای نامتناهی و غیرسببی دارد (سینک). با «بریدن» آن به N ضریب (ضرب در پنجره) یک FIR عملی می‌سازیم؛ هزینهٔ این بریدن، ripple در باند عبور و تضعیف محدود در باند قطع است. انتخاب پنجره، تعادل بین این دو را کنترل می‌کند.

## فرمول‌ها

فیلتر پایین‌گذر با فرکانس قطع `fc` و `M = (N−1)/2`:

```
h_ideal[n] = (2 fc / fs) · sinc(2 fc (n − M) / fs)

سایر انواع:
  High-pass : h = δ[n−M] − h_LP
  Band-pass : h = h_LP(fc2) − h_LP(fc1)
  Band-stop : h = δ[n−M] − h_BP
  Hilbert   : h[n] = (cos(w1 m) − cos(w2 m)) / (π m)        m = n−M
  Diff.     : h[n] = (w2 cos(w2 m) − w1 cos(w1 m) − (sin(w2 m) − sin(w1 m))/m) / (π m)
```

سپس:

```
h[n] = h_ideal[n] · w[n]        (w = پنجره)
h[n] = h[n] / G                 (G = بهرهٔ مرجع نرمال‌سازی)
```

## مثال عددی

`fs = 48000`، `fc = 5000`، `N = 101`، پنجرهٔ Kaiser با `beta = 7.86`:

```c
fce_spec_t sp;
fce_spec_defaults(&sp);
sp.kind = FCE_KIND_FIR;
sp.fir_type = FCE_FIR_LOWPASS;
sp.fs = 48000; sp.fc1 = 5000; sp.num_taps = 101;
sp.window = FCE_WIN_KAISER; sp.kaiser_beta = 7.86;
sp.precision = FCE_PRECISION_FLOAT64;
```

ضرایب تولیدشده (۸ ضریب اول):

```c
static const float64_t coeffs[101] = {
    1.6393231856371196e-05, 1.8045962035311434e-05,
   -5.695443794375429e-20, -4.1619956366239744e-05,
   -9.2849397644657658e-05, -0.00012081213308045559,
   -8.6689229419262297e-05, 2.9387675147424378e-05, ...
};
```

این ضرایب **با `scipy.signal.firwin` یکسان‌اند** (خطای ~1e-17).

## API

```c
fce_status_t fce_generate(const fce_spec_t* spec,
                          fce_result_t* result,
                          fce_workspace_t* ws);
```

فیلدهای مربوط به FIR در result:

| فیلد | معنی |
|---|---|
| `result.h_f64` | ضرایب float64 (همیشه موجود) |
| `result.h_f32` | ضرایب float32 (اگر precision=float32) |
| `result.num_taps` | تعداد ضریب مؤثر |
| `result.symmetry` | نوع تقارن (I/II/III/IV) |
| `result.fir_ideal` | پاسخ ایده‌آل پیش از پنجره (شفافیت) |
| `result.fir_window` | مقادیر پنجره |
| `result.norm_factor` | ضریب نرمال‌سازی اعمال‌شده |

## محاسبهٔ خودکار تعداد taps (Kaiser)

اگر `num_taps = 0` باشد:

```
A = تضعیف باند قطع [dB]
Δf = پهنای باند انتقال [Hz]

N ≈ (A − 7.95) / (2.285 · 2π Δf/fs) + 1     (گرد به بالا + فرد شدن)
β = 0.1102 (A − 8.7)                 اگر A > 50
  = 0.5842 (A − 21)^0.4 + 0.07886 (A − 21)   اگر 21 < A ≤ 50
  = 0                              اگر A ≤ 21
```

مثال: `fs = 48000`، `Δf = 1000 Hz`، `A = 80 dB` ← **N = 243**، **β = 7.857** (مطابق `scipy.signal.kaiserord`).

## انواع تقارن (Symmetry)

| نوع | تعداد ضریب | تقارن | موارد استفاده |
|---|---|---|---|
| I | فرد | متقارن | LP/HP/BP/BS |
| II | زوج | متقارن | LP/BP (دارای صفر در Nyquist) |
| III | فرد | پادمتقارن | Hilbert |
| IV | زوج | پادمتقارن | Differentiator |

کتابخانه به‌صورت خودکار نوع مناسب را انتخاب می‌کند؛ اگر HP/BS با تعداد ضریب زوج بخواهید، پرچم `FCE_FLAG_SYMMETRY_WARNING` ست می‌شود (چون Type II در Nyquist صفر دارد).

## نرمال‌سازی

| استراتژی | بهرهٔ مرجع |
|---|---|
| `FCE_NORM_DC` | `|H(0)| = 1` |
| `FCE_NORM_NYQUIST` | `|H(π fs)| = 1` |
| `FCE_NORM_PASSBAND_PEAK` | بیشینهٔ `|H|` در باند عبور = ۱ (با جستجوی دقیق golden-section) |
| `FCE_NORM_AUTO` | پیش‌فرض بر اساس نوع فیلتر |

## استفادهٔ عملی

```c
/* خروجی مستقیم برای کتابخانهٔ FIR شما */
FilterLab_FIR_SetCoefficients(result.h_f32, result.num_taps);
```

[فصل بعد: پنجره‌ها](04_fir_windows.md)
