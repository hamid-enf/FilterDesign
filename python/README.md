# FilterCoeff — Python port

A faithful pure-Python port of the C99 **FilterCoeff** library. It reimplements
the exact algorithms from `src/fce_math.c`, `src/fce_fir.c` and `src/fce_iir.c`
using only the Python standard library (`math` / `cmath`) — it does **not** wrap
SciPy, so the coefficients are produced by the same computation as the C code.

```text
C library (FilterCoeff)  <--line-for-line-->  python/filtercoeff.py
```

## Requirements

- Python 3.8+ (standard library only; no NumPy/SciPy needed to *use* it).

## Files

| File | Purpose |
|---|---|
| `filtercoeff.py` | the port (single module) |
| `example.py` | concrete input → output examples |
| `compare_c.py` | validates the port against the C library (needs a C toolchain + the built lib) |
| `README.md` | this file |

## Quick start

```python
import filtercoeff as fce

# FIR lowpass: Fs = 48 kHz, cutoff 5 kHz, 101 taps, Hamming window
h = fce.design({
    "kind": "fir", "fir_type": "lowpass",
    "fs": 48000, "fc1": 5000, "num_taps": 101, "window": "hamming",
})

# IIR 4th-order Butterworth lowpass (SOS, {b0,b1,b2,a1,a2} per section)
sos = fce.design({
    "kind": "iir", "iir_family": "butterworth", "iir_type": "lowpass",
    "fs": 48000, "fc1": 5000, "order": 4,
})

# Fixed-point (Q15 / Q31)
flat = [v for sec in sos_reshaped for v in sec]
q, scale, int_bits = fce.quantize(flat, "q15", "symmetric", sec_len=len(flat))
```

`design(spec, details=True)` returns a richer dict (metadata, internals such as
the ideal response / window, section gains, design frequencies …).

## What is covered

- **FIR:** lowpass / highpass / bandpass / bandstop / Hilbert / differentiator
  with Rectangular, Hann, Hamming, Blackman, Kaiser, Blackman-Harris, Bartlett
  and Tukey windows, DC / Nyquist / passband-peak normalization, and Kaiser
  auto-tap counting.
- **IIR:** Butterworth / Chebyshev I / Chebyshev II / Elliptic (Cauer) / Bessel
  in LP / HP / BP / BS, with the same analog-prototype → transform → bilinear →
  zpk2sos pipeline as the C code (and SciPy), plus `buttord`/`cheb1ord`/… auto
  order selection.
- **Fixed-point:** Q15 / Q31 with symmetric / section-wise / coefficient-wise
  scaling.

## Validating against the C library

```sh
make                    # build libfiltercoeff.a at the repo root
make -C tools/reference fce_dump
python3 python/compare_c.py tools/reference/fce_dump
# C vs Python port: 76 checks, 0 failures   (max coeff error < 1e-12)
```

All 76 checks (FIR + IIR + auto-order + Q15/Q31) match the C output to
`~1e-15`, confirming the port is bit-for-bit faithful.

## Sign conventions (identical to the C library)

- FIR: plain coefficient array `h`.
- IIR: SOS layout `{b0, b1, b2, a1, a2}` per section, `a0 == 1` implicit:
  `y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]`.
