# فصل ۶ — خانواده‌های IIR

## مفهوم

هر خانواده یک «پروتوتایپ آنالوگ» متفاوت دارد — یعنی مجموعه قطب‌های متفاوت در صفحهٔ s که رفتار متفاوتی (ripple، تاخیر گروهی، ...) می‌سازند.

## شهود

| خانواده | ایده | بهترین کاربرد |
|---|---|---|
| Butterworth | صاف‌ترین پاسخ ممکن (بدون ripple) | عمومی |
| Chebyshev I | ripple فقط در باند عبور → انتقال تیزتر با همان مرتبه | نیاز به تیز بودن |
| Chebyshev II | ripple فقط در باند قطع | مانند Chebyshev I با رفتار بهتر در باند عبور |
| Elliptic (Cauer) | ripple در هر دو باند → تیزترین انتقال با کمترین مرتبه | کمترین هزینهٔ محاسباتی |
| Bessel | صاف‌ترین تاخیر گروهی | فیلترهای حساس به فاز/پالس |

## فرمول‌ها

### Butterworth — قطب‌های روی دایرهٔ واحد سمت چپ

```
p_k = −exp(j π (2k−1−n) / 2n)      k = 1..n
```

### Chebyshev I — بیضی

```
ε = √(10^(0.1 rp) − 1) ،  μ = asinh(1/ε) / n
p_k = −sinh(μ)·sin(θ_k) − j cosh(μ)·cos(θ_k)      θ_k = π(2k−1)/(2n)
```

### Chebyshev II — معکوس Chebyshev I

```
صفرها:  z_k = j / sin(θ_k)
قطب‌ها: p_k = −1 / sinh(μ + jθ_k)      μ = asinh(1/ε')/n ، ε' = 1/√(10^(0.1 rs) − 1)
```

### Elliptic — توابع بیضوی ژاکوبی (Orfanidis / SciPy)

```
ε² = 10^(0.1 rp) − 1 ،  ε1² = ε² / (10^(0.1 rs) − 1)
m = حل معادلهٔ درجه: n·K(m)/K(1−m) = K(ε1²)/K'(ε1²)      (روش nome)
صفرها:  z_k = j / (√m · sn(j·k·K(m)/n, m))
قطب‌ها:  p_k = −(c·d·sv·cv + j·s·dv) / (1 − (d·sv)²)
بهره:  k = Re(Π(−p)/Π(−z)) ،  و برای n زوج: k /= √(1+ε²)
```

### Bessel — چندجمله‌ای‌های Bessel معکوس

```
قطب‌ها = ریشه‌های θ_n(s) = Σ_k (2n−k)!/(2^(n−k) k! (n−k)!) · s^k
نرمال‌سازی 'mag': مقیاس طوری که |H(jω)| = −3dB در ω = 1
```

ریشه‌یابی با روش **Aberth–Ehrlich** و نقاط شروع Campos–Calderon انجام می‌شود (دقت ~1e-14 برای مرتبه‌های ≤ 8).

## مثال عددی — مقایسهٔ مرتبه برای مشخصات یکسان

مشخصات: `fs=10k، fc=1k، تضعیف ≥60dB در 2k، ripple=1dB`

| خانواده | مرتبهٔ خودکار |
|---|---|
| Butterworth | 10 |
| Chebyshev I | 6 |
| Chebyshev II | 6 |
| Elliptic | 5 |

(مطابق `buttord`/`cheb1ord`/`cheb2ord`/`ellipord`)

## API

```c
sp.iir_family = FCE_IIR_BUTTERWORTH;   /* یا CHEBYSHEV1/CHEBYSHEV2/ELLIPTIC/BESSEL */
sp.passband_ripple_db = 0.5;           /* Chebyshev I و Elliptic */
sp.stopband_atten_db  = 60;            /* Chebyshev II و Elliptic */
```

## کد — طراحی هر ۵ خانواده

```c
for (int fam = FCE_IIR_BUTTERWORTH; fam <= FCE_IIR_BESSEL; fam++)
{
    fce_spec_t sp;
    fce_result_t r;
    uint8_t mem[16384];
    fce_workspace_t ws = { mem, sizeof(mem) };

    fce_spec_defaults(&sp);
    sp.kind = FCE_KIND_IIR;
    sp.iir_family = (fce_iir_family_t)fam;
    sp.iir_type = FCE_IIR_LOWPASS;
    sp.fs = 48000; sp.fc1 = 5000; sp.order = 4;
    sp.passband_ripple_db = 0.5;
    sp.stopband_atten_db = 60;
    sp.precision = FCE_PRECISION_FLOAT64;

    if (fce_generate(&sp, &r, &ws) == FCE_OK)
        printf("family %d: max pole radius = %.9f\n", fam, r.max_pole_radius);
}
```

## ضرایب تولیدشده — نمونه (مرتبهٔ ۴، LP، fs=48k، fc=5k)

| خانواده | سکشن ۱ (b0,b1,b2,a1,a2) |
|---|---|
| Butterworth | `0.00554, 0.01108, 0.00554, -1.01554, 0.28006` |
| Chebyshev I (0.5dB) | `0.00274, 0.00548, 0.00274, -1.13504, 0.31273` |
| Chebyshev II (60dB) | `0.09334, 0.00107, 0.09334, -1.16645, 0.31251` |
| Elliptic (0.5dB/60dB) | `0.00751, 0.01027, 0.00751, -1.42817, 0.56034` |
| Bessel | `0.02899, 0.05798, 0.02899, -1.28065, 0.48290` |

## استفادهٔ عملی

- **دقت فاز:** Bessel؛ **حداقل هزینه:** Elliptic؛ **بدون ripple باند عبور:** Butterworth یا Chebyshev II.
- همیشه `result.max_pole_radius < 1` را چک کنید (کتابخانه این کار را می‌کند).

[فصل بعد: SOS](07_sos_biquad.md)
