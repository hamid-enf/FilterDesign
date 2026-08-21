# پژوهش فنی — منابع و فرمول‌های الگوریتم‌ها

این سند حاصل «پژوهش پیش از پیاده‌سازی» است: برای هر الگوریتم، منبع، فرمول، روش عددی، دقت، محدودیت‌ها و امکان‌سنجی embedded مشخص شده است.

## ۱. Elliptic (Cauer) — پروتوتایپ آنالوگ

| مورد | توضیح |
|---|---|
| منبع | Orfanidis, *Lecture Notes on Elliptic Filter Design*, Rutgers (ece521); Lutovac–Tosic–Evans, *Filter Design for Signal Processing*, ch. 5/12 |
| روش | توابع بیضوی ژاکوبی + معادلهٔ درجه با روش nome |
| فرمول‌ها | `ε²=10^(rp/10)−1`؛ `m` از `n·K(m)/K(1−m)=K(m1)/K(1−m1)` با `q1=exp(−πK'/K)`، `q=q1^(1/n)`، `m=16q(Σq^(k(k+1))/(1+2Σq^(k²)))⁴`؛ صفرها `j/(√m·sn(u,m))`؛ قطب‌ها با `arc_jac_sc1` |
| پیاده‌سازی مرجع | `scipy.signal.ellipap` (بررسی مستقیم کد منبع) |
| دقت | ~1e-15 در پروتوتایپ (مقایسه با SciPy) |
| محدودیت‌ها | `ck1² ≥ 1` نامعتبر (rp ≥ rs)؛ مرتبه‌های خیلی بالا (> 20) حساس |
| Embedded | مناسب؛ فقط چند توابع بیضوی (AGM + Landen) لازم است |

## ۲. Chebyshev I / II

| مورد | توضیح |
|---|---|
| منبع | Parks–Burrus, *Digital Filter Design*; استاندارد |
| فرمول | قطب‌های بیضوی: `−sinh(μ)sin(θ) ∓ j cosh(μ)cos(θ)`؛ Chebyshev II معکوس با صفرهای `j/sin(θ)` |
| دقت | بسته (closed-form)؛ ~1e-15 |
| محدودیت | Chebyshev II مرتبهٔ فرد یک صفر در بی‌نهایت دارد (در SOS با سکشن مرتبهٔ اول نمایان می‌شود) |

## ۳. Butterworth

| مورد | توضیح |
|---|---|
| منبع | استاندارد |
| فرمول | قطب‌ها روی دایرهٔ واحد: `−exp(jπ(2k−1−n)/2n)` |
| دقت | بسته؛ ~1e-16 |

## ۴. Bessel (Thomson)

| مورد | توضیح |
|---|---|
| منبع | Thomson 1949; Bond, *Bessel Filter Constants*; Campos–Calderon 2011 (arXiv:1105.0957) |
| روش | ریشه‌یابی چندجمله‌ای Bessel معکوس `θ_n(s)` با روش **Aberth–Ehrlich** و نقاط شروع Campos؛ نرمال‌سازی 'mag' با دوبخشی (bisection) روی |H(jω)|=−3dB |
| دقت | ~1e-14 (مرتبه ≤ 8)؛ برای مرتبه‌های بالاتر دقت کاهش می‌یابد (~1e-6 در مرتبهٔ ۲۰) |
| محدودیت‌ها | مرتبهٔ خودکار تعریف‌نشده (مانند SciPy)؛ ریشه‌های دور از مبدأ نیازمند شروع خوب |
| Embedded | ریشه‌یابی Aberth چند ده تکرار است — مناسب طراحی‌های گاه‌به‌گاه |

## ۵. تبدیل‌های فرکانس آنالوگ (lp2lp/hp/bp/bs)

| مورد | توضیح |
|---|---|
| منبع | استاندارد (زیرجایگشت s)؛ مطابق `scipy.signal.lp2*_zpk` |
| فرمول | `s→s/ω0`، `s→ω0/s`، `s→(s²+ω0²)/(s·BW)`، `s→s·BW/(s²+ω0²)` |
| نکتهٔ عددی | بهره با مقادیر **پروتوتایپ اصلی** محاسبه می‌شود (قبل از تبدیل) — اشتباه رایجی که در اولین نسخه‌ها رخ داد و با تست مرجع پیدا شد |

## ۶. تبدیل دوخطی (Bilinear/Tustin)

| مورد | توضیح |
|---|---|
| منبع | استاندارد؛ مطابق `scipy.signal.bilinear_zpk` |
| فرمول | `s = 2fs(z−1)/(z+1)`؛ `k' = k·Re(Π(2fs−z)/Π(2fs−p))` — بهره **قبل از** نگاشت |
| دقت | ~1e-15 |
| محدودیت | warp فرکانسی — جبران با prewarp در طراحی‌های بر پایهٔ فرکانس |

## ۷. zpk → SOS (جفت‌کردن nearest)

| مورد | توضیح |
|---|---|
| منبع | `scipy.signal.zpk2sos` (pairing='nearest') — بازتولید مستقل الگوریتم |
| روش | هر بار «بدترین» قطب (نزدیک‌ترین به دایرهٔ واحد) + نزدیک‌ترین صفر؛ موارد خاص: آخرین قطب حقیقی، حفظ صفر حقیقی برای مرتبهٔ فرد |
| نکتهٔ عددی | مزدوج‌ها باید «مصرف» شوند؛ مقادیر حقیقی فقط یک نسخه مصرف می‌کنند (باگ‌های اولیه با تست مرجع پیدا شدند) |
| خروجی | سکشن‌های مرتبهٔ دوم + احتمالاً یک سکشن مرتبهٔ اول |

## ۸. مرتبهٔ خودکار (Order Selection)

| مورد | توضیح |
|---|---|
| منبع | فرمول‌های کلاسیک؛ مطابق `buttord/cheb1ord/cheb2ord/ellipord` |
| فرمول‌ها | Butterworth: `ceil(log10((Gst−1)/(Gp−1))/(2 log10 Ω))`؛ Chebyshev: `ceil(acosh(√((Gst−1)/(Gp−1)))/acosh(Ω))`؛ Elliptic: نسبت انتگرال‌های بیضوی کامل |
| نکته | فرکانس‌های prewarp شده؛ برای bandpass/bandstop از `nat = min(گزینه‌های دو لبه)` |
| دقت | تطابق کامل با SciPy در تست‌های مرجع |

## ۹. Kaiser (تخمین taps و β)

| مورد | توضیح |
|---|---|
| منبع | Oppenheim & Schafer, *Discrete-Time Signal Processing*, pp. 475-476; `scipy.signal.kaiserord/kaiser_beta` |
| فرمول | `N = (A−7.95)/(2.285·Δω) + 1`؛ β سه‌تکه‌ای |
| دقت | تطابق کامل (`kaiserord(80, ...)` → 243) |

## ۱۰. توابع بیضوی (K, sn/cn/dn, asn)

| مورد | توضیح |
|---|---|
| منبع | Abramowitz & Stegun 17.6 (AGM/Landen); Orfanidis eq. 56 |
| روش | K با AGM؛ sn/cn/dn با Landen نزولی؛ asn با Landen صعودی + asin مختلط |
| نکتهٔ عددی | `K(1−m)` نباید با `1−m` محاسبه شود (خطای 4.7e-7 برای m کوچک) — از AGM با `√m` استفاده می‌شود |
| دقت | ~1e-15 در مقابل `scipy.special` |

## ۱۱. ریشه‌یابی چندجمله‌ای (Bessel)

| مورد | توضیح |
|---|---|
| منبع | Aberth 1973; Ehrlich 1967; Campos–Calderon 2011 |
| روش | Aberth–Ehrlich هم‌زمان: `r_i ← r_i − w/(1 − w·Σ1/(r_i−r_j))` |
| چرا نه Durand-Kerner | DK روی چندجمله‌ای‌های Bessel دچار رکود (stagnation) می‌شود (با آزمایش تجربی نشان داده شد) |
| دقت | همگرایی در ۳–۸ تکرار؛ ~1e-14 |

## ۱۲. Fixed-Point

| مورد | توضیح |
|---|---|
| روش | round-to-nearest + اشباع؛ سه استراتژی مقیاس (سراسری/بخش‌به‌بخش/ضریب‌به‌ضریب) |
| محدودیت ذاتی | فیلترهای با `1−|p| < 2^−F` در Q15 ریسک ناپایداری دارند — کتابخانه تشخیص و گزارش می‌دهد |
| اعتبارسنجی | بازسازی `q/scale` و مقایسه با float؛ پایداری قطب‌های کمی‌شده |

## منابع اصلی

1. Orfanidis, *Lecture Notes on Elliptic Filter Design*, https://www.ece.rutgers.edu/~orfanidi/ece521/notes.pdf
2. Oppenheim & Schafer, *Discrete-Time Signal Processing*, 3rd ed.
3. Parks & Burrus, *Digital Filter Design*
4. Lutovac, Tosic, Evans, *Filter Design for Signal Processing*
5. Thomson, *Delay Networks having Maximally Flat Frequency Characteristics*, Proc. IEE 1949
6. Bond, *Bessel Filter Constants*, crbond.com/papers/bsf.pdf
7. Campos & Calderon, *Approximate closed-form formulas for the zeros of the Bessel Polynomials*, arXiv:1105.0957
8. Aberth, *Iteration Methods for Finding all Zeros of a Polynomial Simultaneously*, Math. Comp. 1973
9. SciPy source (`scipy/signal/_filter_design.py`) — مرجع پیاده‌سازی و تست
