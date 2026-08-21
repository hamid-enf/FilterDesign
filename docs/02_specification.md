# فصل ۲ — نحوه تعریف Specification

## مفهوم

همه‌چیز با یک ساختار `fce_spec_t` شروع می‌شود. این ساختار تمام خواسته‌های کاربر را توصیف می‌کند؛ سپس `fce_generate()` آن را به ضرایب تبدیل می‌کند.

```c
typedef struct fce_spec
{
    fce_kind_t kind;              /* FIR یا IIR */
    double fs;                    /* نرخ نمونه‌برداری [Hz] */
    double fc1, fc2;              /* فرکانس‌های قطع/لبه */
    fce_precision_t precision;    /* float32 یا float64 */
    fce_norm_t normalization;     /* استراتژی نرمال‌سازی */
    fce_qformat_t qformat;        /* Q15 / Q31 / none */
    fce_scale_strategy_t scale_strategy;
    uint32_t validate;            /* کدام بررسی‌ها اجرا شوند */

    fce_fir_type_t fir_type;      /* FIR */
    uint16_t num_taps;            /* 0 = خودکار (Kaiser) */
    fce_window_t window;
    double kaiser_beta;           /* 0 = خودکار از روی تضعیف */
    double stopband_atten_db;     /* برای Kaiser خودکار */
    double transition_hz;         /* برای تعداد taps خودکار */

    fce_iir_family_t iir_family;  /* IIR */
    fce_iir_type_t iir_type;
    uint16_t order;               /* 0 = خودکار */
    double passband_ripple_db;    /* gpass */
    double stopband_atten_db;     /* gstop */
    double edge1_hz, edge2_hz;    /* لبهٔ باند مقابل (مرتبهٔ خودکار) */
    fce_sos_order_t sos_order;
} fce_spec_t;
```

## شهود

سه راه برای پر کردن spec وجود دارد:

1. **کاملاً دستی** — همهٔ فیلدها را خودتان ست کنید.
2. **کمک‌کننده‌ها** — `fce_spec_fir()` و `fce_spec_iir()` برای پر کردن سریع.
3. **خودکار** — `num_taps = 0` یا `order = 0` بگذارید تا کتابخانه محاسبه کند.

## فرمول‌ها و قواعد

| فیلد | قانون اعتبارسنجی |
|---|---|
| `fs` | باید `> 0` باشد |
| `fc1` | برای FIR/IIR عادی: `0 < fc1 < fs/2` |
| `fc2` | فقط برای BP/BS: `fc1 < fc2 < fs/2` |
| `num_taps` | اگر ۰ باشد باید `window == FCE_WIN_KAISER` و `stopband_atten_db > 0` و `transition_hz > 0` باشد |
| `order` | اگر ۰ باشد باید `edge1_hz` (و برای BP/BS: `edge2_hz`) معتبر باشد |
| `passband_ripple_db` | برای Chebyshev I و Elliptic الزامی است |
| `stopband_atten_db` | برای Chebyshev II و Elliptic الزامی است |

## مثال‌های کد

### FIR با کمک‌کننده

```c
fce_spec_t sp;
fce_spec_fir(&sp,
             FCE_FIR_LOWPASS,   /* نوع */
             48000.0,           /* fs  */
             5000.0,            /* fc1 */
             0.0,               /* fc2 (بی‌استفاده) */
             101,               /* taps */
             FCE_WIN_KAISER,    /* پنجره */
             7.86,              /* beta */
             FCE_PRECISION_FLOAT32);
```

### IIR با مرتبهٔ خودکار

```c
fce_spec_t sp;
fce_spec_defaults(&sp);
sp.kind = FCE_KIND_IIR;
sp.iir_family = FCE_IIR_BUTTERWORTH;
sp.iir_type  = FCE_IIR_LOWPASS;
sp.fs        = 10000;
sp.fc1       = 1000;      /* لبهٔ passband */
sp.order     = 0;         /* خودکار */
sp.edge1_hz  = 2000;      /* لبهٔ stopband */
sp.passband_ripple_db = 3.0;
sp.stopband_atten_db  = 60;
```

نتیجه: مرتبهٔ ۹ (مطابق `scipy.signal.buttord`).

### خروجی Q15

```c
sp.qformat = FCE_QFORMAT_Q15;
sp.scale_strategy = FCE_SCALE_SYMMETRIC;  /* یا SECTION_WISE */
```

## اعتبارسنجی ورودی

`fce_generate()` ابتدا spec را بررسی می‌کند. اگر نامعتبر باشد:

| کد | معنی |
|---|---|
| `FCE_ERR_INVALID_ARGUMENT` | اشاره‌گر NULL |
| `FCE_ERR_INVALID_SPEC` | مقادیر نامعتبر (مثل `fs <= 0`) |
| `FCE_ERR_UNSUPPORTED` | ویژگی‌ای که در build غیرفعال است |
| `FCE_ERR_BUFFER_TOO_SMALL` | workspace کافی نیست |

## نکتهٔ مهم: تعیین اندازهٔ workspace

```c
size_t need = fce_workspace_required(&spec);  /* دقیق برای این spec */
size_t need = fce_workspace_required_max();   /* برای بدترین حالت */
```

## استفادهٔ عملی

- همیشه اول `fce_spec_defaults(&sp)` را صدا بزنید تا مقادیر پیش‌فرض (precision=float32، normalization=Auto، اعتبارسنجی کامل) ست شوند.
- برای پروژه‌های embedded کوچک، `FCE_MAX_FIR_TAPS` و `FCE_MAX_IIR_ORDER` را در `filtercoeff_config.h` کم کنید تا workspace کوچک‌تر شود.

[فصل بعد: FIR](03_fir_generation.md)
