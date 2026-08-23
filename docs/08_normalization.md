# فصل ۸ — نرمال‌سازی

## مفهوم

«نرمال‌سازی» یعنی تعیین اینکه بهرهٔ واحد در چه فرکانسی تعریف شود. بعد از طراحی، ضرایب طوری مقیاس می‌شوند که `|H(f_ref)| = 1`.

## شهود

یک فیلتر پایین‌گذر پنجره‌شده به‌صورت طبیعی بهرهٔ DC حدود ۱ دارد اما نه دقیقاً. اگر بخواهید سیگنال بعد از فیلتر هم‌دامنه بماند، باید DC را دقیقاً ۱ کنید. برای فیلتر بالاگذر، مرجع Nyquist است؛ برای باندگذر، پیک باند عبور.

## استراتژی‌ها

| مقدار | مرجع | کاربرد |
|---|---|---|
| `FCE_NORM_AUTO` | بر اساس نوع فیلتر | پیش‌فرض |
| `FCE_NORM_DC` | `|H(0)| = 1` | LP، BS |
| `FCE_NORM_NYQUIST` | `|H(π·fs)| = 1` | HP، Differentiator |
| `FCE_NORM_PASSBAND_PEAK` | بیشینهٔ `|H|` در باند عبور = ۱ | BP، Hilbert |
| `FCE_NORM_NONE` | بدون مقیاس | نیاز به بهرهٔ خام |

Auto برای هر نوع: LP→DC، HP→Nyquist، BP→PassbandPeak، BS→DC، Hilbert→PassbandPeak، Differentiator→Nyquist.

## فرمول‌ها

### FIR

```
G = |Σ h[n] e^(−jω_ref n)|          (ω_ref = 0 یا π یا فرکانس پیک)
h'[n] = h[n] / G
norm_factor = 1/G  →  در result.norm_factor گزارش می‌شود
```

پیک باند عبور با اسکن روی شبکه + بهینه‌سازی golden-section پیدا می‌شود (دقت ~1e-12). عرض شبکه متناسب با ریپل فیلتر است — `max(256, 4·N)` نقطه (سقف 4096) — چون دورهٔ تناوب ریپل حدود `fs/(N−1)` است و شبکهٔ ثابتِ خیلی درشت ممکن است پیک واقعی را برای فیلترهای بلند (N≳100) جا بیندازد.

## مرجع روی صفر فیلتر (مورد خاص)

اگر مرجع نرمال‌سازی دقیقاً روی یک صفرِ پاسخ بیفتد، بهرهٔ مرجع صفر است و تقسیم بر آن بی‌معنی:

- **HP/BS با تعداد ضریب زوج**: فیلتر Type-II حتماً در Nyquist صفر دارد.
- **Differentiator با تعداد ضریب فرد**: Type-III در Nyquist صفر دارد.
- هر مرجع صریحی که در عمق باند قطع بیفتد.

رفتار کتابخانه:

| وضعیت | رفتار |
|---|---|
| `FCE_NORM_AUTO` | سقوط خودکار به `FCE_NORM_PASSBAND_PEAK` + ست شدن `FCE_FLAG_SYMMETRY_WARNING` (طراحی معتبر می‌ماند) |
| مرجع صریح روی صفر دقیق | `FCE_ERR_NUMERICAL` |
| مرجع صریح در عمق باند قطع (تقویت > ~1e8×) | انجام می‌شود ولی `FCE_FLAG_NUMERICAL_WARNING` ست می‌شود |

`result.normalization` همیشه مرجعی را نشان می‌دهد که واقعاً استفاده شده (بعد از سقوط AUTO ممکن است با `sp.normalization` فرق کند).

### IIR

بهرهٔ کل داخل ضرایب SOS است (`result.norm_factor = 1`). نرمال‌سازی به صورت ضمنی در محاسبهٔ بهرهٔ `k` در زنجیرهٔ zpk انجام می‌شود: پروتوتایپ‌ها طوری تعریف می‌شوند که DC (یا فرکانس مرجع خانواده) بهرهٔ ۱ داشته باشد.

## مثال عددی — FIR LP با Hamming

`fs=48k, fc=5k, N=101`، بدون نرمال‌سازی: `Σh = 0.99997` → با `FCE_NORM_DC`:

```
norm_factor = 1.00003
Σh' = 1.0000000000
```

## API

```c
sp.normalization = FCE_NORM_DC;         /* یا سایر مقادیر */
...
double f = result.norm_factor;          /* ضریب اعمال‌شده */
```

## کد — مقایسهٔ استراتژی‌ها

```c
fce_spec_t sp;
fce_result_t r;
uint8_t mem[16384];
fce_workspace_t ws = { mem, sizeof(mem) };

fce_spec_defaults(&sp);
sp.kind = FCE_KIND_FIR;
sp.fir_type = FCE_FIR_BANDPASS;
sp.fs = 48000; sp.fc1 = 3000; sp.fc2 = 6000; sp.num_taps = 121;
sp.window = FCE_WIN_HAMMING;
sp.normalization = FCE_NORM_PASSBAND_PEAK;   /* پیک باند عبور = ۱ */
sp.precision = FCE_PRECISION_FLOAT64;
fce_generate(&sp, &r, &ws);
printf("norm factor = %.6f\n", r.norm_factor);   /* ~0.9973 */
```

## استفادهٔ عملی

- برای فیلترهایی که بعداً بهرهٔ دیجیتال/آنالوگ دارند (مثل اکولایزر)، `FCE_NORM_NONE` بگذارید و بهره را خودتان مدیریت کنید.
- برای fixed-point، نرمال‌سازی روی بهرهٔ DC/پیک اثر مستقیم روی استفاده از محدودهٔ Q دارد: پیک ۱ یعنی استفادهٔ کامل از 0..32767.

[فصل بعد: ضرایب اعشاری](09_float_coefficients.md)
