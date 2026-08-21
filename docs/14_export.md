# فصل ۱۴ — خروجی (Export)

## مفهوم

ضرایب باید «آمادهٔ کدنویسی» باشند. چهار خروجی وجود دارد:

| خروجی | کاربرد |
|---|---|
| **C arrays** | قرار دادن مستقیم در کد هدف |
| **CSV** | Excel / ابزارهای خارجی |
| **JSON** | ابزارها و اسکریپت‌ها |
| **Markdown report** | مستندسازی و بازبینی انسانی |

## شهود

همهٔ خروجی‌ها از یک **نویسندهٔ متنی** (`fce_writer_t`) استفاده می‌کنند؛ یعنی روی هر پلتفرمی کار می‌کنند: بافر RAM روی MCU، `FILE*` روی Host، UART و...

## API

### حافظه (embedded)

```c
char buf[16384];
fce_mem_writer_t mw;
fce_writer_t w;
fce_writer_mem_init(&w, &mw, buf, sizeof(buf));

fce_export_c_fir(&result, NULL, &w);
/* buf آمادهٔ استفاده است (NUL-terminated) */
```

### فایل (host)

```c
FILE* f = fopen("coeffs.h", "w");
fce_writer_t w;
fce_writer_file_init(&w, f);
fce_export_c_sos(&result, NULL, &w);
fclose(f);
```

## خروجی C — FIR

```c
static const float32_t coeffs[101] = {
    1.6393232e-05f,
    1.8045962e-05f,
    ...
};

static const int16_t coeffs_q15[101] = {   /* اگر Q15 خواسته باشید */
    3, 3, 0, -7, ...
};
```

## خروجی C — SOS

```c
/* y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2] */
static const float32_t sos[][5] = {
    { 0.00554177f, 0.01108354f, 0.00554177f, -1.01554283f, 0.28006372f },
    { 1.00000000f, 2.00000000f, 1.00000000f, -1.28690544f, 0.62210881f }
};
```

## خروجی CSV

```
section,b0,b1,b2,a1,a2
0,0.0055417715,0.0110835431,0.0055417715,-1.0155428256,0.2800637169
1,1,2,1,-1.2869054403,0.6221088069
```

## خروجی JSON

```json
{
  "library": "FilterCoeff",
  "kind": "IIR",
  "family": "Elliptic",
  "order": 4,
  "num_sections": 2,
  "sos_f64": [ [...], [...] ],
  "qformat": "Q15",
  ...
}
```

## خروجی Markdown (گزارش)

```markdown
# FILTER COEFFICIENT REPORT

**Type:** IIR Lowpass
**Family:** Elliptic
**Sample rate:** 48000.0000 Hz
**Order / Sections:** 4 / 2
**Stopband attenuation:** 60.3 dB
...
**Validation: PASS**
```

## گزینه‌ها

```c
fce_export_opts_t o;
memset(&o, 0, sizeof(o));
o.name = "my_coeffs";              /* نام آرایه */
o.precision = 9;                   /* رقم‌های اعشار (0 = خودکار) */
o.include_quantized = 1;           /* آرایه‌های Q15/Q31 هم تولید شود */
```

## استفادهٔ عملی — اتصال به کتابخانهٔ خارجی

```c
/* ۱) تولید خروجی */
fce_export_c_fir(&result, NULL, &w);

/* ۲) کپی خروجی در کد هدف */
static const float32_t coeffs[101] = { ... };

/* ۳) تحویل به runtime شما */
FilterLab_FIR_SetCoefficients(coeffs, 101);
```

[فصل بعد: STM32H7](15_stm32h7.md)
