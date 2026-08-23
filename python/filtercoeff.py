"""
filtercoeff.py - faithful pure-Python port of the FilterCoeff C99 library.

This module reproduces the coefficient-generation algorithms in
``src/fce_math.c``, ``src/fce_fir.c`` and ``src/fce_iir.c`` using nothing
but the Python standard library (``math`` / ``cmath``). It is a *port* of
the C code, not a wrapper around SciPy, so every design is computed by the
same algorithm and produces the same coefficients (to ~1e-15).

Design flow (mirrors ``fce_generate``):

    spec -> auto taps/order -> FIR/IIR design -> coefficients

Sign conventions (identical to the C library):

    FIR:      plain coefficient array ``h`` (windowed-sinc).
    IIR:      SOS with per-section layout ``{b0, b1, b2, a1, a2}``
              (``a0 == 1`` implicit):
              y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
                     - a1*y[n-1] - a2*y[n-2]

API
---

    coeffs = design(spec)   # spec is a dict; see below

    result = design(spec, details=True)
    # -> {"kind":..., "h":..., "sos":..., "num_taps":..., "order":...,
    #      "num_sections":..., "kaiser_beta":..., "norm_factor":...,
    #      "symmetry":..., "design_fc1":..., "design_fc2":...,
    #      "order_clamped":... / "taps_clamped":...}  (auto-design limits hit)

    # Low-level helpers are exposed too: windows, prototypes, transforms,
    # so you can build custom pipelines exactly like the C library.

Spec fields (a plain dict; sensible defaults for omitted keys)::

    kind         : "fir" | "iir"
    fs, fc1, fc2 : sample rate and cutoff / band edges in Hz
    precision    : "float32" | "float64"  (float64 default; float32 only
                   rounds the final coefficients)
    normalization: "auto" | "dc" | "nyquist" | "passband_peak" | "none"
    # FIR
    fir_type     : "lowpass" | "highpass" | "bandpass" | "bandstop"
                   | "hilbert" | "differentiator"
    num_taps     : int (0 = auto via Kaiser)
    window       : "rectangular"|"hann"|"hamming"|"blackman"|"kaiser"
                   |"blackman_harris"|"bartlett"|"tukey"
    kaiser_beta  : float (0 = auto from stopband_atten_db)
    transition_hz: Hz for Kaiser auto-taps
    stopband_atten_db : dB for Kaiser beta/auto-taps
    # IIR
    iir_family   : "butterworth"|"chebyshev1"|"chebyshev2"|"elliptic"|"bessel"
    iir_type     : "lowpass"|"highpass"|"bandpass"|"bandstop"
    order        : int (0 = auto)
    passband_ripple_db : dB (cheby1/ellip)
    stopband_atten_db  : dB (cheby2/ellip + auto-order)
    edge1_hz, edge2_hz : opposite-band edges (auto-order only)
    sos_order    : "default"|"pole_radius_asc"|"pole_radius_desc"|"internal_gain"
"""

import math
import cmath

# ----------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------
PI = 3.14159265358979323846264338327950288
CPLX_TOL = 1e-9

MAX_FIR_TAPS = 2048
MAX_IIR_ORDER = 32
MAX_SECTIONS = (MAX_IIR_ORDER + 1) // 2
MAX_AUTO_ORDER = 24

_NAN = complex(float("nan"), float("nan"))


def _isfinite_c(z):
    return math.isfinite(z.real) and math.isfinite(z.imag)


# ======================================================================
# math core (port of fce_math.c)
# ======================================================================

def sinc(x):
    """sinc(x) = sin(pi*x)/(pi*x),  sinc(0) = 1  (numpy convention)."""
    if x == 0.0:
        return 1.0
    px = PI * x
    return math.sin(px) / px


def i0(x):
    """Modified Bessel I0 (series), matches fce_i0."""
    s = 1.0
    term = 1.0
    x2 = 0.25 * x * x
    k = 1
    while k < 200:
        term *= x2 / (float(k) * float(k))
        s += term
        if term < 1e-18 * s:
            break
        k += 1
    return s


def ellipk(m):
    """Complete elliptic integral K(m), 0 <= m < 1 (AGM)."""
    if m <= 0.0:
        return PI / 2.0
    if m >= 1.0:
        return 1.0e300
    a = 1.0
    b = math.sqrt(1.0 - m)
    while True:
        t = 0.5 * (a + b)
        b = math.sqrt(a * b)
        a = t
        if a - b < 1e-15 * a:
            break
    return PI / (2.0 * a)


def ellipkm1(m):
    """K(1-m) computed stably for small m (AGM(1, sqrt(m)))."""
    if m <= 0.0:
        return 1.0e300
    if m >= 1.0:
        return PI / 2.0
    a = 1.0
    b = math.sqrt(m)
    while True:
        t = 0.5 * (a + b)
        b = math.sqrt(a * b)
        a = t
        if a - b < 1e-15 * a:
            break
    return PI / (2.0 * a)


def ellipj(u, m):
    """Jacobi elliptic functions (sn, cn, dn) for real u, 0 <= m <= 1."""
    if m <= 0.0:
        return (math.sin(u), math.cos(u), 1.0)
    if m >= 1.0:
        th = math.tanh(u)
        c = 1.0 / math.cosh(u)
        return (th, c, c)

    a = [1.0]
    b = math.sqrt(1.0 - m)
    c = [math.sqrt(m)]
    n = 0
    while n < 63:
        an = 0.5 * (a[n] + b)
        cn = 0.5 * (a[n] - b)
        b = math.sqrt(a[n] * b)
        a.append(an)
        c.append(cn)
        n += 1
        if c[n] < 1e-16 * a[n]:
            break

    phi = math.ldexp(a[n], n) * u
    for i in range(n, 0, -1):
        t = (c[i] / a[i]) * math.sin(phi)
        if t > 1.0:
            t = 1.0
        if t < -1.0:
            t = -1.0
        phi = 0.5 * (math.asin(t) + phi)

    sn = math.sin(phi)
    cn = math.cos(phi)
    dn = math.sqrt(1.0 - m * sn * sn)
    return (sn, cn, dn)


def _cx_complement(kx):
    """sqrt(1 - kx^2), stable for small kx."""
    return math.sqrt((1.0 - kx) * (1.0 + kx))


def arc_jac_sn(w, m):
    """Inverse Jacobi sn for complex w (ascending Landen), as in SciPy."""
    k = math.sqrt(m)
    if k > 1.0:
        return _NAN
    if k == 1.0:
        return cmath.atanh(w)  # sn(z,1)=tanh(z) -> z=atanh(w)

    ks = [k]
    niter = 0
    while ks[niter] != 0.0:
        kp = _cx_complement(ks[niter])
        ks.append((1.0 - kp) / (1.0 + kp))
        niter += 1
        if niter >= 62:
            return _NAN

    K = 1.0
    for i in range(1, niter + 1):
        K *= (1.0 + ks[i])
    K *= PI / 2.0

    wns = [w]
    for i in range(niter):
        num = 2.0 * wns[i]
        arg = wns[i] * ks[i]
        cplx_c = cmath.sqrt(1.0 - arg * arg)
        den = (1.0 + ks[i + 1]) * (1.0 + cplx_c)
        wns.append(num / den)
        if not _isfinite_c(wns[i + 1]):
            return _NAN

    u = (2.0 / PI) * cmath.asin(wns[niter])
    return u * K


def ellipdeg(n, m1):
    """Solve the elliptic degree equation n*K(m)/K(1-m) = K(m1)/K(1-m1)."""
    if n == 0:
        return float("nan")
    K1 = ellipk(m1)
    K1p = ellipkm1(m1)
    q1 = math.exp(-PI * K1p / K1)
    q = q1 ** (1.0 / float(n))
    num = 0.0
    for k in range(0, 9):
        num += q ** (float(k * (k + 1)))
    den = 1.0
    for k in range(1, 10):
        den += 2.0 * q ** (float(k * k))
    return 16.0 * q * (num / den) ** 4.0


def _polyval(x, c):
    """c[0] + c[1] x + ... + c[n-1] x^(n-1) (Horner)."""
    s = 0.0
    for ci in reversed(c):
        s = s * x + ci
    return s


def roots_durand_kerner(coeff, deg, start=None):
    """Durand-Kerner / Aberth roots of a real-coefficient polynomial.
    coeff[k] = coefficient of x^k. Returns a list of `deg` complex roots.
    Raises ValueError on non-convergence."""
    if deg == 0 or deg > MAX_IIR_ORDER:
        raise ValueError("bad degree")

    lead = coeff[deg]
    if lead == 0.0:
        raise ValueError("leading coefficient is zero")
    a = [coeff[i] / lead for i in range(deg + 1)]

    if start is not None:
        r = list(start)
    else:
        scale = sum(abs(a[i]) for i in range(deg + 1)) / float(deg + 1)
        if not scale > 0.0:
            raise ValueError("zero polynomial")
        r = []
        for i in range(deg):
            ang = 2.0 * PI * float(i) / float(deg)
            rad = 0.7 + 0.3 * float((i * 7) % 5) / 4.0
            r.append(complex(rad * math.cos(ang) * scale,
                             rad * math.sin(ang) * scale))

    converged = False
    for _ in range(200):
        converged = True
        for i in range(deg):
            pv = 1.0 + 0.0j
            dp = 0.0 + 0.0j
            for j in range(deg - 1, -1, -1):
                dp = dp * r[i] + pv
                pv = pv * r[i] + a[j]
            ssum = 0.0 + 0.0j
            for j in range(deg):
                if j != i:
                    ssum += 1.0 / (r[i] - r[j])
            w = pv / dp
            den = 1.0 - w * ssum
            if abs(den) < 1e-300:
                raise ValueError("denominator too small")
            step = w / den
            r[i] = r[i] - step
            if abs(step) > 1e-12 * (1.0 + abs(r[i])):
                converged = False
        if converged:
            break

    if not converged:
        raise ValueError("root finder did not converge")
    return r


def campos_zeros(n):
    """Approximate zeros of the Bessel polynomial (Campos & Calderon)."""
    if n == 1:
        return [complex(-1.0, 0.0)]
    s = _polyval(float(n), [0.0, 0.0, 2.0, 0.0, -3.0, 1.0])
    b3 = _polyval(float(n), [16.0, -8.0]) / s
    b2 = _polyval(float(n), [-24.0, -12.0, 12.0]) / s
    b1 = _polyval(float(n), [8.0, 24.0, -12.0, -2.0]) / s
    b0 = _polyval(float(n), [0.0, -6.0, 0.0, 5.0, -1.0]) / s
    r = _polyval(float(n), [0.0, 0.0, 2.0, 1.0])
    a1 = _polyval(float(n), [-6.0, -6.0]) / r
    a2 = 6.0 / r
    z = []
    for k in range(1, n + 1):
        kk = float(k)
        z.append(complex(a2 * kk * kk + a1 * kk,
                         ((b3 * kk + b2) * kk + b1) * kk + b0))
    return z


def bessel_poles(n):
    """Poles of the reverse Bessel polynomial (delay normalization)."""
    if n == 0 or n > MAX_IIR_ORDER:
        raise ValueError("bad n")
    a = [0.0] * (n + 1)
    for k in range(n + 1):
        num = 1.0
        for i in range(2 * n - k):
            num *= float(2 * n - k - i)
        den = math.ldexp(1.0, n - k)
        for i in range(2, k + 1):
            den *= float(i)
        for i in range(2, n - k + 1):
            den *= float(i)
        a[k] = num / den
    yz = campos_zeros(n)
    start = [1.0 / y for y in yz]
    poles = roots_durand_kerner(a, n, start)
    ssum = sum(p.real for p in poles)
    if abs(ssum + a[n - 1] / a[n]) > 1e-8 * (1.0 + abs(ssum)):
        raise ValueError("bessel poles sanity failed")
    return poles


# ----------------------------------------------------------------------
# biquad helpers
# ----------------------------------------------------------------------

def eval_biquad(c, w):
    """H(z) = (b0+b1 z^-1+b2 z^-2)/(1+a1 z^-1+a2 z^-2) at z^-1 = e^-jw.
    c = [b0, b1, b2, a1, a2]."""
    cw = math.cos(w)
    sw = math.sin(w)
    c2w = math.cos(2.0 * w)
    s2w = math.sin(2.0 * w)
    num = complex(c[0] + c[1] * cw + c[2] * c2w,
                  -(c[1] * sw + c[2] * s2w))
    den = complex(1.0 + c[3] * cw + c[4] * c2w,
                  -(c[3] * sw + c[4] * s2w))
    return num / den


def biquad_poles(a1, a2):
    """Roots of z^2 + a1 z + a2 (our sign convention). Returns (p1, p2)."""
    disc = a1 * a1 - 4.0 * a2
    if disc >= 0.0:
        s = math.sqrt(disc)
        return (complex(0.5 * (-a1 + s), 0.0),
                complex(0.5 * (-a1 - s), 0.0))
    s = math.sqrt(-disc)
    return (complex(-0.5 * a1, 0.5 * s),
            complex(-0.5 * a1, -0.5 * s))


def biquad_peak_gain(c, grid=256):
    """Peak |H| of a biquad over [0, pi] on a `grid`-point scan."""
    best = 0.0
    for i in range(grid + 1):
        w = PI * float(i) / float(grid)
        g = abs(eval_biquad(c, w))
        if g > best:
            best = g
    return best


def biquad_group_delay(c, w):
    """Finite-difference-free group delay [samples] of a biquad at w."""
    cw = math.cos(w)
    sw = math.sin(w)
    c2w = math.cos(2.0 * w)
    s2w = math.sin(2.0 * w)
    br = c[0] + c[1] * cw + c[2] * c2w
    bi = -(c[1] * sw + c[2] * s2w)
    bdr = -(c[1] * sw + 2.0 * c[2] * s2w)
    bdi = -(c[1] * cw + 2.0 * c[2] * c2w)
    ar = 1.0 + c[3] * cw + c[4] * c2w
    ai = -(c[3] * sw + c[4] * s2w)
    adr = -(c[3] * sw + 2.0 * c[4] * s2w)
    adi = -(c[3] * cw + 2.0 * c[4] * c2w)
    bm2 = br * br + bi * bi
    am2 = ar * ar + ai * ai
    gd_b = (br * bdi - bi * bdr) / bm2 if bm2 > 0.0 else 0.0
    gd_a = (ar * adi - ai * adr) / am2 if am2 > 0.0 else 0.0
    return gd_b - gd_a


# ----------------------------------------------------------------------
# Kaiser parameters
# ----------------------------------------------------------------------

def kaiser_beta(atten_db):
    a = abs(atten_db)
    if a > 50.0:
        return 0.1102 * (a - 8.7)
    if a > 21.0:
        return 0.5842 * (a - 21.0) ** 0.4 + 0.07886 * (a - 21.0)
    return 0.0


def kaiser_taps(atten_db, transition_hz, fs):
    a = abs(atten_db)
    dw = 2.0 * PI * transition_hz / fs
    if dw <= 0.0 or fs <= 0.0:
        return 0
    n = (a - 7.95) / (2.285 * dw) + 1.0
    n = max(3.0, n)  # N >= 2 required (windows divide by N-1); stay odd
    taps = int(math.ceil(n))
    if taps & 1 == 0:
        taps += 1
    return taps


# ======================================================================
# FIR (port of fce_fir.c)
# ======================================================================

def window_value(win, n, N, kaiser_beta, tukey_alpha=0.5):
    """Value of the window `win` at index n of an N-point window."""
    if win == "rectangular":
        return 1.0
    if win == "hann":
        return 0.5 * (1.0 - math.cos(2.0 * PI * float(n) / float(N - 1)))
    if win == "hamming":
        return 0.54 - 0.46 * math.cos(2.0 * PI * float(n) / float(N - 1))
    if win == "blackman":
        return (0.42 - 0.5 * math.cos(2.0 * PI * float(n) / float(N - 1))
                + 0.08 * math.cos(4.0 * PI * float(n) / float(N - 1)))
    if win == "kaiser":
        t = 2.0 * float(n) / float(N - 1) - 1.0
        return i0(kaiser_beta * math.sqrt(1.0 - t * t)) / i0(kaiser_beta)
    if win == "blackman_harris":
        return (0.35875
                - 0.48829 * math.cos(2.0 * PI * float(n) / float(N - 1))
                + 0.14128 * math.cos(4.0 * PI * float(n) / float(N - 1))
                - 0.01168 * math.cos(6.0 * PI * float(n) / float(N - 1)))
    if win == "bartlett":
        x = 2.0 * float(n) / float(N - 1) - 1.0
        return 1.0 - abs(x)
    if win == "tukey":
        alpha = tukey_alpha
        nn = float(n)
        nmax = float(N - 1)
        if alpha <= 0.0:
            return 1.0
        if alpha >= 1.0:
            return 0.5 * (1.0 - math.cos(2.0 * PI * nn / nmax))
        if nn < alpha * nmax * 0.5:
            return 0.5 * (1.0 + math.cos(PI * (2.0 * nn / (alpha * nmax) - 1.0)))
        if nn >= nmax * (1.0 - alpha * 0.5):
            return 0.5 * (1.0 + math.cos(PI * (2.0 * nn / (alpha * nmax)
                                               - 2.0 / alpha + 1.0)))
        return 1.0
    raise ValueError("unknown window: " + str(win))


def _fir_ideal(fir_type, fs, fc1, fc2, N):
    """Ideal (pre-window) impulse response."""
    M = 0.5 * float(N - 1)
    h = [0.0] * N
    for n in range(N):
        m = float(n) - M
        v = 0.0
        if fir_type == "lowpass":
            v = (2.0 * fc1 / fs) * sinc(2.0 * fc1 * m / fs)
        elif fir_type == "highpass":
            # delta - LP; the band-limited delta is sinc(m), which is
            # nonzero at the half-integer m of even tap counts
            v = sinc(m) - (2.0 * fc1 / fs) * sinc(2.0 * fc1 * m / fs)
        elif fir_type == "bandpass":
            v = ((2.0 * fc2 / fs) * sinc(2.0 * fc2 * m / fs)
                 - (2.0 * fc1 / fs) * sinc(2.0 * fc1 * m / fs))
        elif fir_type == "bandstop":
            v = (sinc(m)
                 - (2.0 * fc2 / fs) * sinc(2.0 * fc2 * m / fs)
                 + (2.0 * fc1 / fs) * sinc(2.0 * fc1 * m / fs))
        elif fir_type == "hilbert":
            w1 = (2.0 * PI * fc1 / fs) if (fc1 > 0.0 and fc2 > fc1) else 0.0
            w2 = (2.0 * PI * fc2 / fs) if (fc1 > 0.0 and fc2 > fc1) else PI
            if m == 0.0:
                v = 0.0
            else:
                v = (math.cos(w1 * m) - math.cos(w2 * m)) / (PI * m)
        elif fir_type == "differentiator":
            w1 = (2.0 * PI * fc1 / fs) if (fc1 > 0.0 and fc2 > fc1) else 0.0
            w2 = (2.0 * PI * fc2 / fs) if (fc1 > 0.0 and fc2 > fc1) else PI
            if m == 0.0:
                v = 0.0
            else:
                v = (w2 * math.cos(w2 * m) - w1 * math.cos(w1 * m)
                     - (math.sin(w2 * m) - math.sin(w1 * m)) / m) / (PI * m)
        else:
            raise ValueError("unknown fir_type: " + str(fir_type))
        h[n] = v
    return h


def _fir_gain_at(h, w):
    """|H(e^jw)| = |sum h[n] e^{-jwn}|."""
    re = 0.0
    im = 0.0
    for n, hn in enumerate(h):
        a = w * float(n)
        re += hn * math.cos(a)
        im -= hn * math.sin(a)
    return math.sqrt(re * re + im * im)


def _fir_peak_in_band(h, fs, f_lo, f_hi):
    # grid scales with N so one bracket stays below a ripple period
    # (~fs/(N-1) Hz); golden-section needs a unimodal bracket
    grid = max(256, 4 * len(h))
    grid = min(grid, 4096)
    f_lo_c = 0.0 if f_lo <= 0.0 else f_lo
    f_hi_c = (0.5 * fs) if f_hi >= 0.5 * fs else f_hi
    if f_hi_c <= f_lo_c:
        f_hi_c = 0.5 * fs

    f_best = 0.0
    best = 0.0
    for i in range(grid + 1):
        f = f_lo_c + (f_hi_c - f_lo_c) * float(i) / float(grid)
        g = _fir_gain_at(h, 2.0 * PI * f / fs)
        if g > best:
            best = g
            f_best = f

    gr = 0.6180339887498948482
    a = f_best - (f_hi_c - f_lo_c) / float(grid)
    b = f_best + (f_hi_c - f_lo_c) / float(grid)
    if a < f_lo_c:
        a = f_lo_c
    if b > f_hi_c:
        b = f_hi_c
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = _fir_gain_at(h, 2.0 * PI * c / fs)
    fd = _fir_gain_at(h, 2.0 * PI * d / fs)
    for _ in range(60):
        if fc > fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = _fir_gain_at(h, 2.0 * PI * c / fs)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = _fir_gain_at(h, 2.0 * PI * d / fs)
        if b - a < 1e-12 * (1.0 + abs(b)):
            break
    return 0.5 * (fc + fd)


_FIR_TYPE_NORM = {
    "lowpass": "dc",
    "highpass": "nyquist",
    "bandpass": "passband_peak",
    "bandstop": "dc",
    "hilbert": "passband_peak",
    "differentiator": "nyquist",
}


def _fir_passband(spec):
    """Passband interval [f_lo, f_hi] for passband-peak normalization."""
    fs = spec["fs"]
    fc1 = spec["fc1"]
    fc2 = spec.get("fc2", 0.0)
    fir_type = spec["fir_type"]
    if fir_type == "bandpass":
        return fc1, fc2
    if fir_type == "highpass":
        return fc1, 0.5 * fs
    if fir_type == "hilbert":
        if fc1 > 0.0 and fc2 > fc1:
            return fc1, fc2
        return 0.0, 0.5 * fs
    # LP / BS / differentiator: full band (BS passbands are outside
    # [fc1, fc2]; scanning the stopband there was a bug)
    return 0.0, 0.5 * fs


def _norm_gain_bad(norm_gain, hmax):
    return (not (norm_gain > 0.0)) or (not (norm_gain < 1e300)) or \
           (not (norm_gain > 1e-12 * hmax))


def fir_design(spec):
    """Port of fce_fir_design. Returns dict with h, internals, metadata."""
    fs = spec["fs"]
    fir_type = spec["fir_type"]
    window = spec.get("window", "hamming")
    fc1 = spec["fc1"]
    fc2 = spec.get("fc2", 0.0)
    N = spec["num_taps"]
    beta = 0.0

    if window == "kaiser":
        kb = spec.get("kaiser_beta", 0.0)
        if kb > 0.0:
            beta = kb
        elif spec.get("stopband_atten_db", 0.0) > 0.0:
            beta = kaiser_beta(spec["stopband_atten_db"])
        else:
            beta = 0.0

    ideal = _fir_ideal(fir_type, fs, fc1, fc2, N)
    win = [window_value(window, n, N, beta, spec.get("tukey_alpha", 0.5))
           for n in range(N)]
    h = [ideal[n] * win[n] for n in range(N)]

    norm = spec.get("normalization", "auto")
    if norm == "auto":
        norm = _FIR_TYPE_NORM.get(fir_type, "dc")

    if norm == "dc":
        norm_gain = _fir_gain_at(h, 0.0)
    elif norm == "nyquist":
        norm_gain = _fir_gain_at(h, PI)
    elif norm == "passband_peak":
        f_lo, f_hi = _fir_passband(spec)
        norm_gain = _fir_peak_in_band(h, fs, f_lo, f_hi)
    else:  # "none"
        norm_gain = 1.0

    # degenerate reference (a symmetry null on the AUTO reference edge):
    # AUTO falls back to the passband peak; an explicit degenerate
    # request is an error. (mirrors fce_fir_design)
    hmax = max(abs(hn) for hn in h)
    if (_norm_gain_bad(norm_gain, hmax)
            and spec.get("normalization", "auto") == "auto"
            and norm != "passband_peak"):
        norm = "passband_peak"
        f_lo, f_hi = _fir_passband(spec)
        norm_gain = _fir_peak_in_band(h, fs, f_lo, f_hi)
    if _norm_gain_bad(norm_gain, hmax):
        raise ValueError("numerical error in normalization")
    h = [hn / norm_gain for hn in h]

    precision = spec.get("precision", "float64")
    if precision == "float32":
        h = [float(x) for x in h]

    # anti-symmetric responses follow tap parity exactly like C:
    # odd taps -> Type III, even taps -> Type IV
    if N & 1:
        sym = "I" if fir_type in ("lowpass", "highpass", "bandpass",
                                  "bandstop") else "III"
    else:
        sym = "II" if fir_type in ("lowpass", "highpass", "bandpass",
                                   "bandstop") else "IV"

    return {"kind": "fir",
            "h": h,
            "num_taps": N,
            "kaiser_beta": beta,
            "norm_factor": 1.0 / norm_gain,
            "normalization": norm,
            "symmetry": sym,
            "fir_ideal": ideal,
            "fir_window": win}


# ======================================================================
# IIR prototypes (port of fce_iir.c)
# ======================================================================

def proto_butter(n):
    p = []
    for i in range(n):
        m = 2 * i - n + 1
        th = PI * float(m) / (2.0 * float(n))
        p.append(complex(-math.cos(th), -math.sin(th)))
    return p, 1.0


def proto_cheb1(n, rp):
    eps = math.sqrt(10.0 ** (0.1 * rp) - 1.0)
    mu = math.asinh(1.0 / eps) / float(n)
    prod = 1.0
    p = []
    for i in range(n):
        m = 2 * i - n + 1
        th = PI * float(m) / (2.0 * float(n))
        sh = math.sinh(mu)
        ch = math.cosh(mu)
        p.append(complex(-sh * math.cos(th), -ch * math.sin(th)))
        prod *= (sh * math.cos(th)) ** 2 + (ch * math.sin(th)) ** 2
    k = math.sqrt(prod)
    if n % 2 == 0:
        k /= math.sqrt(1.0 + eps * eps)
    return p, k


def proto_cheb2(n, rs):
    de = 1.0 / math.sqrt(10.0 ** (0.1 * rs) - 1.0)
    mu = math.asinh(1.0 / de) / float(n)
    z = []
    p = []
    nz_total = (n - 1) if (n & 1) else n
    for i in range(nz_total):
        if n & 1:
            half = (n - 1) // 2
            if i < half:
                m = -(2 * (half - i))
            else:
                m = 2 * (i - half + 1)
        else:
            m = 2 * i - n + 1
        th = PI * float(m) / (2.0 * float(n))
        z.append(complex(0.0, 1.0 / math.sin(th)))
    for i in range(n):
        m = 2 * i - n + 1
        th = PI * float(m) / (2.0 * float(n))
        sh = math.sinh(mu)
        ch = math.cosh(mu)
        s = complex(sh * math.cos(th), ch * math.sin(th))
        p.append(-1.0 / s)
    prod_p = 1.0
    prod_z = 1.0
    for pi in p:
        prod_p *= (-pi)
    for zi in z:
        prod_z *= (-zi)
    k = (prod_p / prod_z).real
    return p, z, k


def proto_ellip(n, rp, rs):
    if n == 1:
        p0 = -math.sqrt(1.0 / (10.0 ** (0.1 * rp) - 1.0))
        return [complex(p0, 0.0)], [], -p0

    eps_sq = 10.0 ** (0.1 * rp) - 1.0
    eps = math.sqrt(eps_sq)
    ck1_sq = eps_sq / (10.0 ** (0.1 * rs) - 1.0)
    if not (ck1_sq > 0.0 and ck1_sq < 1.0):
        raise ValueError("invalid elliptic spec (ck1_sq)")

    val0 = ellipk(ck1_sq)
    val1 = ellipkm1(ck1_sq)
    m = ellipdeg(n, ck1_sq)
    if not (m > 0.0 and m < 1.0):
        raise ValueError("elliptic degree equation failed")
    capk = ellipk(m)

    # zeros
    z = []
    for j in range((n + 1) // 2):
        jv = (2 * j) if (n & 1) else (2 * j + 1)
        u = float(jv) * capk / float(n)
        s, c, d = ellipj(u, m)
        if abs(s) > 1e-14:
            z.append(complex(0.0, 1.0 / (math.sqrt(m) * s)))
    cnt = len(z)
    for j in range(cnt):
        z.append(z[j].conjugate())

    # r = Im(asn(j/eps, ck1_sq))
    w = arc_jac_sn(complex(0.0, 1.0 / eps), ck1_sq)
    if not _isfinite_c(w) or abs(w.real) > 1e-9:
        raise ValueError("arc_jac failed")
    r = w.imag
    v0 = capk * r / (float(n) * val0)

    # poles
    sv, cv, dv = ellipj(v0, 1.0 - m)
    p = []
    for j in range((n + 1) // 2):
        jv = (2 * j) if (n & 1) else (2 * j + 1)
        u = float(jv) * capk / float(n)
        s, c, d = ellipj(u, m)
        num = complex(-(c * d * sv * cv), -(s * dv))
        den = 1.0 - (d * sv) * (d * sv)
        if abs(den) < 1e-300:
            raise ValueError("elliptic pole denominator")
        p.append(num / den)

    if n & 1:
        norm2 = sum(abs(pi) ** 2 for pi in p)
        keep = [pi for pi in p if abs(pi.imag) > 1e-14 * math.sqrt(norm2)]
        for pi in keep:
            p.append(pi.conjugate())
    else:
        base = len(p)
        for j in range(base):
            p.append(p[j].conjugate())

    prod_p = 1.0
    prod_z = 1.0
    for pi in p:
        prod_p *= (-pi)
    for zi in z:
        prod_z *= (-zi)
    k = (prod_p / prod_z).real
    if n % 2 == 0:
        k /= math.sqrt(1.0 + eps_sq)
    return p, z, k


def proto_bessel(n):
    p = bessel_poles(n)
    a_last = 1.0
    for i in range(n):
        a_last *= float(2 * n - i)
    for i in range(n):
        a_last /= 2.0

    target = 1.0 / math.sqrt(2.0)
    lo = 0.1
    hi = 1.5
    for _ in range(200):
        prod = 1.0
        for pi in p:
            prod *= (1j * hi - pi)
        if a_last / abs(prod) < target:
            break
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        prod_lo = 1.0
        for pi in p:
            prod_lo *= (1j * lo - pi)
        g_lo = a_last / abs(prod_lo)
        prod_mid = 1.0
        for pi in p:
            prod_mid *= (1j * mid - pi)
        g_mid = a_last / abs(prod_mid)
        if (g_lo - target) * (g_mid - target) <= 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-14 * hi:
            break
    norm_factor = 0.5 * (lo + hi)

    p = [pi / norm_factor for pi in p]
    k = a_last * norm_factor ** (-float(n))
    return p, k


# ----------------------------------------------------------------------
# analog frequency transformations (port)
# ----------------------------------------------------------------------

def tr_lp2lp(z, p, k, wo):
    nz, np_ = len(z), len(p)
    z = [zi * wo for zi in z]
    p = [pi * wo for pi in p]
    k *= wo ** float(np_ - nz)
    return z, p, k


def tr_lp2hp(z, p, k, wo):
    nz, np_ = len(z), len(p)
    deg = np_ - nz
    pz = 1.0
    pp = 1.0
    for zi in z:
        pz *= (-zi)
    for pi in p:
        pp *= (-pi)
    k *= (pz / pp).real
    z = [wo / zi for zi in z]
    p = [wo / pi for pi in p]
    z += [0.0 + 0.0j] * deg
    return z, p, k


def tr_lp2bp(z, p, k, wo, bw):
    nz, np_ = len(z), len(p)
    deg = np_ - nz
    zo = []
    po = []
    for zi in z:
        a = zi * (0.5 * bw)
        s = cmath.sqrt(a * a - wo * wo)
        zo.append(a + s)
    for zi in z:
        a = zi * (0.5 * bw)
        s = cmath.sqrt(a * a - wo * wo)
        zo.append(a - s)
    for pi in p:
        a = pi * (0.5 * bw)
        s = cmath.sqrt(a * a - wo * wo)
        po.append(a + s)
    for pi in p:
        a = pi * (0.5 * bw)
        s = cmath.sqrt(a * a - wo * wo)
        po.append(a - s)
    nzo = 2 * nz + deg
    npo = 2 * np_
    while len(zo) < nzo:
        zo.append(0.0 + 0.0j)
    ko = k * bw ** float(deg)
    return zo, po, ko


def tr_lp2bs(z, p, k, wo, bw):
    nz, np_ = len(z), len(p)
    deg = np_ - nz
    zo = []
    po = []
    for zi in z:
        a = (0.5 * bw) / zi
        s = cmath.sqrt(a * a - wo * wo)
        zo.append(a + s)
    for zi in z:
        a = (0.5 * bw) / zi
        s = cmath.sqrt(a * a - wo * wo)
        zo.append(a - s)
    for pi in p:
        a = (0.5 * bw) / pi
        s = cmath.sqrt(a * a - wo * wo)
        po.append(a + s)
    for pi in p:
        a = (0.5 * bw) / pi
        s = cmath.sqrt(a * a - wo * wo)
        po.append(a - s)
    pz = 1.0
    pp = 1.0
    for zi in z:
        pz *= (-zi)
    for pi in p:
        pp *= (-pi)
    ko = k * (pz / pp).real
    for _ in range(deg):
        zo.append(complex(0.0, wo))
        zo.append(complex(0.0, -wo))
    return zo, po, ko


def tr_bilinear(z, p, k, fs):
    nz, np_ = len(z), len(p)
    deg = np_ - nz
    fs2 = 2.0 * fs
    pz = 1.0
    pp = 1.0
    for zi in z:
        pz *= (fs2 - zi)
    for pi in p:
        pp *= (fs2 - pi)
    k *= (pz / pp).real
    z = [(fs2 + zi) / (fs2 - zi) for zi in z]
    p = [(fs2 + pi) / (fs2 - pi) for pi in p]
    z += [-1.0 + 0.0j] * deg
    return z, p, k


# ----------------------------------------------------------------------
# zpk -> SOS (port of scipy-style nearest pairing)
# ----------------------------------------------------------------------

def _is_real(c):
    return abs(c.imag) <= 1e-12 * (1.0 + abs(c.real))


def _dist2(a, b):
    dr = a.real - b.real
    di = a.imag - b.imag
    return dr * dr + di * di


def _consume(arr, used, n, v):
    consumed = 0
    need_conj = not _is_real(v)
    target = 2 if need_conj else 1
    for i in range(n):
        if consumed >= target:
            break
        if used[i]:
            continue
        if _dist2(arr[i], v) < 1e-20:
            used[i] = 1
            consumed += 1
        elif need_conj and _dist2(arr[i], v.conjugate()) < 1e-20:
            used[i] = 1
            consumed += 1
    return consumed


def _nearest_zero(z, used, nz, p1, kind):
    best = nz
    best_d = 0.0
    for i in range(nz):
        if used[i]:
            continue
        if kind == 1 and not _is_real(z[i]):
            continue
        if kind == 2 and _is_real(z[i]):
            continue
        d = _dist2(z[i], p1)
        if best == nz or d < best_d:
            best = i
            best_d = d
    return best


def _section_make(zs, ps):
    b = [0.0, 0.0, 1.0]
    a = [0.0, 0.0, 1.0]
    if len(zs) == 1:
        b[0] = 0.0
        b[1] = 1.0
        b[2] = -zs[0].real
    elif len(zs) == 2:
        b[0] = 1.0
        b[1] = -(zs[0].real + zs[1].real)
        b[2] = zs[0].real * zs[1].real - zs[0].imag * zs[1].imag
    if len(ps) == 1:
        a[0] = 0.0
        a[1] = 1.0
        a[2] = -ps[0].real
    elif len(ps) == 2:
        a[0] = 1.0
        a[1] = -(ps[0].real + ps[1].real)
        a[2] = ps[0].real * ps[1].real - ps[0].imag * ps[1].imag
    return [b[0], b[1], b[2], a[1], a[2]]


def zpk2sos(z, p, k):
    """Port of the zpk->sos 'nearest' pairing loop. Returns list of
    [b0, b1, b2, a1, a2] sections."""
    dp = list(p)
    dz = list(z)
    npad = len(dp)
    nzpad = len(dz)
    while nzpad < npad:
        dz.append(0.0 + 0.0j)
        nzpad += 1
    while npad < nzpad:
        dp.append(0.0 + 0.0j)
        npad += 1
    if (npad & 1) != 0:
        dp.append(0.0 + 0.0j)
        dz.append(0.0 + 0.0j)
        npad += 1
        nzpad += 1
    ns = (npad + 1) // 2
    # BP/BS double the prototype order, so ns can reach `order` itself
    # (the old cap at (order+1)//2 wrongly rejected BP/BS orders > 16)
    if ns > MAX_IIR_ORDER:
        raise ValueError("too many sections")

    up = [0] * npad
    uz = [0] * nzpad
    np_rem = npad
    nz_rem = nzpad
    sos = [None] * ns

    for si in range(ns - 1, -1, -1):
        p1i = npad
        worst = 1e300
        for j in range(npad):
            if up[j]:
                continue
            d = abs(1.0 - abs(dp[j]))
            if p1i == npad or d < worst:
                p1i = j
                worst = d
        if p1i == npad:
            raise ValueError("no remaining pole")
        p1 = dp[p1i]
        np_rem -= _consume(dp, up, npad, p1)

        nreal_p = sum(1 for j in range(npad) if not up[j] and _is_real(dp[j]))
        nreal_z = sum(1 for j in range(nzpad) if not uz[j] and _is_real(dz[j]))

        if _is_real(p1) and nreal_p == 0:
            z1i = _nearest_zero(dz, uz, nzpad, p1, 1)
            if z1i < nzpad:
                nz_rem -= _consume(dz, uz, nzpad, dz[z1i])
            zs = [dz[z1i] if z1i < nzpad else (0.0 + 0.0j),
                  (0.0 + 0.0j)]
            ps = [p1, (0.0 + 0.0j)]
            sos[si] = _section_make(zs, ps)
        elif (np_rem + 1 == nz_rem) and (not _is_real(p1)) \
                and nreal_p == 1 and nreal_z == 1:
            z1i = _nearest_zero(dz, uz, nzpad, p1, 2)
            if z1i == nzpad:
                raise ValueError("pairing failed")
            nz_rem -= _consume(dz, uz, nzpad, dz[z1i])
            zs = [dz[z1i], dz[z1i].conjugate()]
            ps = [p1, p1.conjugate()]
            sos[si] = _section_make(zs, ps)
        else:
            if _is_real(p1):
                p2i = npad
                worst2 = 1e300
                for j in range(npad):
                    if up[j] or not _is_real(dp[j]):
                        continue
                    d = abs(1.0 - abs(dp[j]))
                    if p2i == npad or d < worst2:
                        p2i = j
                        worst2 = d
                if p2i == npad:
                    raise ValueError("no real pole 2")
                p2 = dp[p2i]
                np_rem -= _consume(dp, up, npad, p2)
            else:
                p2 = p1.conjugate()

            z1i = _nearest_zero(dz, uz, nzpad, p1, 0)
            if z1i == nzpad:
                ps = [p1, p2]
                sos[si] = _section_make([], ps)
                continue
            z1 = dz[z1i]
            nz_rem -= _consume(dz, uz, nzpad, z1)
            zs = []
            ps = []
            if not _is_real(z1):
                zs = [z1, z1.conjugate()]
            else:
                zs = [z1]
                z2i = _nearest_zero(dz, uz, nzpad, p1, 1)
                if z2i < nzpad:
                    zs.append(dz[z2i])
                    nz_rem -= _consume(dz, uz, nzpad, dz[z2i])
            ps = [p1, p2]
            sos[si] = _section_make(zs, ps)

    return sos, k


# ----------------------------------------------------------------------
# IIR auto order (port of fce_iir_auto_order)
# ----------------------------------------------------------------------

def _prewarp(f, fs):
    return 2.0 * fs * math.tan(PI * f / fs)


def _acosh(x):
    return math.log(x + math.sqrt(x * x - 1.0))


def _bs_obj(wp, ind, passb0, passb1, stopb0, stopb1, gpass, gstop,
            iir_family):
    """scipy band_stop_obj: fractional bandstop order as one passband
    edge is moved (the other stays fixed)."""
    p0 = wp if ind == 0 else passb0
    p1 = wp if ind == 1 else passb1
    n1 = abs(stopb0 * (p0 - p1) / (stopb0 * stopb0 - p0 * p1))
    n2 = abs(stopb1 * (p0 - p1) / (stopb1 * stopb1 - p0 * p1))
    nat = min(n1, n2)
    if iir_family == "butterworth":
        return (math.log10((gstop - 1.0) / (gpass - 1.0))
                / (2.0 * math.log10(nat)))
    if iir_family == "elliptic":
        arg1_sq = (gpass - 1.0) / (gstop - 1.0)
        arg0 = 1.0 / nat
        d00 = ellipk(arg0 * arg0)
        d01 = ellipkm1(arg0 * arg0)
        d10 = ellipk(arg1_sq)
        d11 = ellipkm1(arg1_sq)
        return d00 * d11 / (d01 * d10)
    # cheby 1 & 2
    return (_acosh(math.sqrt((gstop - 1.0) / (gpass - 1.0)))
            / _acosh(nat))


def _bs_fminbound(ind, x1, x2, passb0, passb1, stopb0, stopb1,
                  gpass, gstop, iir_family):
    """Bounded Brent minimization; 1:1 port of scipy
    _minimize_scalar_bounded (xatol=1e-5, maxfun=500)."""
    sqrt_eps = 1.4901161193847656e-08      # sqrt(2.2e-16)
    golden_mean = 0.38196601125010510      # 0.5*(3-sqrt(5))
    xatol = 1e-5
    a, b = x1, x2

    def f(x):
        return _bs_obj(x, ind, passb0, passb1, stopb0, stopb1,
                       gpass, gstop, iir_family)

    fulc = nfc = xf = a + golden_mean * (b - a)
    fx = f(xf)
    ffulc = fnfc = fx
    rat = e = 0.0
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * abs(xf) + xatol / 3.0
    tol2 = 2.0 * tol1
    num = 1
    while abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = True
        if abs(e) > tol1:
            golden = False
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = abs(q)
            r = e
            e = rat
            if (abs(p) < abs(0.5 * q * r)) and (p > q * (a - xf)) \
                    and (p < q * (b - xf)):
                rat = p / q
                x = xf + rat
                if (x - a) < tol2 or (b - x) < tol2:
                    si = 1.0 if xm >= xf else -1.0
                    rat = tol1 * si
            else:
                golden = True
        if golden:
            e = (a - xf) if xf >= xm else (b - xf)
            rat = golden_mean * e
        si = (1.0 if rat > 0.0 else (-1.0 if rat < 0.0 else 1.0))
        x = xf + si * max(abs(rat), tol1)
        fu = f(x)
        num += 1
        if num >= 500:
            break
        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu
        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * abs(xf) + xatol / 3.0
        tol2 = 2.0 * tol1
    return xf


def iir_auto_order(spec):
    """Returns dict(order, design_fc1, design_fc2). Raises ValueError."""
    fs = spec["fs"]
    iir_type = spec["iir_type"]
    iir_family = spec["iir_family"]
    gpass_db = spec.get("passband_ripple_db", 0.0)
    if not gpass_db > 0.0:
        gpass_db = 3.0
    gstop_db = spec.get("stopband_atten_db", 0.0)
    if not gstop_db > gpass_db:
        raise ValueError("need gstop > gpass")

    if iir_type == "lowpass":
        wp1 = spec["fc1"]
        ws1 = spec["edge1_hz"]
        if not ws1 > wp1:
            raise ValueError("lowpass auto-order edges")
        passb1 = _prewarp(wp1, fs)
        stopb1 = _prewarp(ws1, fs)
        nat = stopb1 / passb1
        design_fc1 = design_fc2 = 0.0
    elif iir_type == "highpass":
        wp1 = spec["fc1"]
        ws1 = spec["edge1_hz"]
        if not ws1 < wp1:
            raise ValueError("highpass auto-order edges")
        passb1 = _prewarp(wp1, fs)
        stopb1 = _prewarp(ws1, fs)
        nat = passb1 / stopb1
        design_fc1 = design_fc2 = 0.0
    elif iir_type == "bandpass":
        wp1, wp2 = spec["fc1"], spec["fc2"]
        ws1, ws2 = spec["edge1_hz"], spec["edge2_hz"]
        if not (ws1 < wp1 and wp2 < ws2):
            raise ValueError("bandpass auto-order edges")
        passb1 = _prewarp(wp1, fs)
        passb2 = _prewarp(wp2, fs)
        stopb1 = _prewarp(ws1, fs)
        stopb2 = _prewarp(ws2, fs)
        n1 = abs((stopb1 * stopb1 - passb1 * passb2) /
                 (stopb1 * (passb1 - passb2)))
        n2 = abs((stopb2 * stopb2 - passb1 * passb2) /
                 (stopb2 * (passb1 - passb2)))
        nat = min(n1, n2)
        design_fc1 = design_fc2 = 0.0
    elif iir_type == "bandstop":
        ws1, ws2 = spec["fc1"], spec["fc2"]
        wp1, wp2 = spec["edge1_hz"], spec["edge2_hz"]
        if not (wp1 < ws1 and ws2 < wp2):
            raise ValueError("bandstop auto-order edges")
        passb1 = _prewarp(wp1, fs)
        passb2 = _prewarp(wp2, fs)
        stopb1 = _prewarp(ws1, fs)
        stopb2 = _prewarp(ws2, fs)
        # scipy _find_nat_freq: each passband edge is optimized so the
        # fractional order is minimal (mirrors band_stop_obj/fminbound)
        gp = 10.0 ** (0.1 * gpass_db)
        gs = 10.0 ** (0.1 * gstop_db)
        opt0 = _bs_fminbound(0, passb1, stopb1 - 1e-12,
                             passb1, passb2, stopb1, stopb2,
                             gp, gs, iir_family)
        opt1 = _bs_fminbound(1, stopb2 + 1e-12, passb2,
                             passb1, passb2, stopb1, stopb2,
                             gp, gs, iir_family)
        passb1, passb2 = opt0, opt1
        n1 = abs((stopb1 * (passb1 - passb2)) /
                 (stopb1 * stopb1 - passb1 * passb2))
        n2 = abs((stopb2 * (passb1 - passb2)) /
                 (stopb2 * stopb2 - passb1 * passb2))
        nat = min(n1, n2)
        design_fc1 = design_fc2 = 0.0
    else:
        raise ValueError("bad iir_type")

    if not nat > 1.0:
        raise ValueError("bad nat")
    gpass = 10.0 ** (0.1 * gpass_db)
    gstop = 10.0 ** (0.1 * gstop_db)

    if iir_family == "butterworth":
        n = int(math.ceil(math.log10((gstop - 1.0) / (gpass - 1.0)) /
                          (2.0 * math.log10(nat))))
    elif iir_family in ("chebyshev1", "chebyshev2"):
        n = int(math.ceil(_acosh(math.sqrt((gstop - 1.0) / (gpass - 1.0))) /
                          _acosh(nat)))
    elif iir_family == "elliptic":
        arg1_sq = (gpass - 1.0) / (gstop - 1.0)
        arg0 = 1.0 / nat
        d00 = ellipk(arg0 * arg0)
        d01 = ellipkm1(arg0 * arg0)
        d10 = ellipk(arg1_sq)
        d11 = ellipkm1(arg1_sq)
        n = int(math.ceil(d00 * d11 / (d01 * d10)))
    elif iir_family == "bessel":
        raise ValueError("bessel has no auto-order rule")
    else:
        raise ValueError("bad iir_family")

    if n < 1:
        raise ValueError("bad order")
    clamped = 0
    if n > MAX_AUTO_ORDER:
        n = MAX_AUTO_ORDER
        clamped = 1  # surfaced like the C FCE_FLAG_ORDER_CLAMPED

    w0 = 1.0
    if iir_family == "butterworth":
        w0 = (gpass - 1.0) ** (-1.0 / (2.0 * float(n)))

    if iir_type == "lowpass":
        if iir_family == "chebyshev2":
            design_fc1 = spec["edge1_hz"]
        else:
            wc = passb1 * w0
            design_fc1 = fs * math.atan(wc / (2.0 * fs)) / PI
    elif iir_type == "highpass":
        if iir_family == "chebyshev2":
            design_fc1 = spec["edge1_hz"]
        else:
            wc = passb1 / w0
            design_fc1 = fs * math.atan(wc / (2.0 * fs)) / PI
    elif iir_type == "bandpass":
        bw = passb2 - passb1
        w0b = 1.0 if iir_family == "chebyshev2" else w0
        wn_hi = math.sqrt(w0b * w0b * bw * bw / 4.0 + passb1 * passb2) \
            + w0b * bw * 0.5
        wn_lo = math.sqrt(w0b * w0b * bw * bw / 4.0 + passb1 * passb2) \
            - w0b * bw * 0.5
        if iir_family == "chebyshev2":
            design_fc1 = spec["edge1_hz"]
            design_fc2 = spec["edge2_hz"]
        else:
            design_fc1 = fs * math.atan(wn_lo / (2.0 * fs)) / PI
            design_fc2 = fs * math.atan(wn_hi / (2.0 * fs)) / PI
    elif iir_type == "bandstop":
        if iir_family == "chebyshev2":
            # cheb2ord filter_type=3: -gpass frequency via new_freq
            v = _acosh(math.sqrt((gstop - 1.0) / (gpass - 1.0)))
            new_freq = 1.0 / math.cosh(v / float(n))
            nat0 = (new_freq * 0.5 * (passb1 - passb2)
                    + math.sqrt(new_freq * new_freq * (passb2 - passb1)
                                * (passb2 - passb1) * 0.25
                                + passb1 * passb2))
            nat1 = passb1 * passb2 / nat0
            lo, hi = min(nat0, nat1), max(nat0, nat1)
            design_fc1 = fs * math.atan(lo / (2.0 * fs)) / PI
            design_fc2 = fs * math.atan(hi / (2.0 * fs)) / PI
        else:
            # butter: discr formula; cheb1/ellip: the optimized edges
            bw = passb2 - passb1
            discr = math.sqrt(bw * bw + 4.0 * w0 * w0 * passb1 * passb2)
            wn_hi = (bw + discr) / (2.0 * w0)
            wn_lo = (discr - bw) / (2.0 * w0)
            design_fc1 = fs * math.atan(wn_lo / (2.0 * fs)) / PI
            design_fc2 = fs * math.atan(wn_hi / (2.0 * fs)) / PI

    return {"order": n, "design_fc1": design_fc1, "design_fc2": design_fc2,
            "clamped": clamped}


def iir_design(spec):
    """Port of fce_iir_design. Returns dict with sos, order, metadata."""
    fs = spec["fs"]
    iir_family = spec["iir_family"]
    iir_type = spec["iir_type"]
    n = spec["order"]

    # --- analog prototype ---
    if iir_family == "butterworth":
        pp, k = proto_butter(n)
        pz, nz = [], 0
    elif iir_family == "chebyshev1":
        pp, k = proto_cheb1(n, spec.get("passband_ripple_db", 0.0))
        pz, nz = [], 0
    elif iir_family == "chebyshev2":
        pp, pz, k = proto_cheb2(n, spec.get("stopband_atten_db", 0.0))
        nz = len(pz)
    elif iir_family == "elliptic":
        pp, pz, k = proto_ellip(n, spec.get("passband_ripple_db", 0.0),
                                spec.get("stopband_atten_db", 0.0))
        nz = len(pz)
    elif iir_family == "bessel":
        pp, k = proto_bessel(n)
        pz, nz = [], 0
    else:
        raise ValueError("bad iir_family")
    np_ = len(pp)

    # --- analog frequency transformation ---
    design_fc1 = spec.get("design_fc1", spec["fc1"])
    design_fc2 = spec.get("design_fc2", spec.get("fc2", 0.0))
    if iir_type == "lowpass":
        z, p, k = tr_lp2lp(pz, pp, k, _prewarp(design_fc1, fs))
    elif iir_type == "highpass":
        z, p, k = tr_lp2hp(pz, pp, k, _prewarp(design_fc1, fs))
    elif iir_type == "bandpass":
        w1 = _prewarp(design_fc1, fs)
        w2 = _prewarp(design_fc2, fs)
        wo = math.sqrt(w1 * w2)
        bw = w2 - w1
        z, p, k = tr_lp2bp(pz, pp, k, wo, bw)
    elif iir_type == "bandstop":
        w1 = _prewarp(design_fc1, fs)
        w2 = _prewarp(design_fc2, fs)
        wo = math.sqrt(w1 * w2)
        bw = w2 - w1
        z, p, k = tr_lp2bs(pz, pp, k, wo, bw)
    else:
        raise ValueError("bad iir_type")

    # --- bilinear ---
    z, p, k = tr_bilinear(z, p, k, fs)

    # --- zpk -> sos ---
    sos, k = zpk2sos(z, p, k)
    ns = len(sos)

    # --- section ordering + gain distribution ---
    sos_order = spec.get("sos_order", "default")
    peaks = [biquad_peak_gain(sos[i], 256) for i in range(ns)]
    if sos_order in ("pole_radius_desc", "internal_gain"):
        keys = []
        for i in range(ns):
            if sos_order == "internal_gain":
                keys.append(peaks[i])
            else:
                p1, p2 = biquad_poles(sos[i][3], sos[i][4])
                keys.append(max(abs(p1), abs(p2)))
        # simple selection sort, mirroring the C
        for i in range(ns):
            best = i
            for j in range(i + 1, ns):
                take = (keys[j] > keys[best]) \
                    if sos_order == "pole_radius_desc" \
                    else (keys[j] < keys[best])
                if take:
                    best = j
            if best != i:
                sos[i], sos[best] = sos[best], sos[i]
                keys[i], keys[best] = keys[best], keys[i]
                peaks[i], peaks[best] = peaks[best], peaks[i]

    for i in range(ns):
        g = k if i == 0 else 1.0
        for j in range(3):
            sos[i][j] *= g
        peaks[i] *= g

    return {"kind": "iir",
            "sos": sos,
            "order": n,
            "num_sections": ns,
            "design_fc1": design_fc1,
            "design_fc2": design_fc2,
            "norm_factor": 1.0,
            "section_gains": peaks}


# ======================================================================
# fixed-point quantization (port of fce_quant.c)
# ======================================================================

def _round_half(x):
    """round-half-away-from-zero (matches fce_round_half)."""
    return math.floor(x + 0.5) if x >= 0.0 else math.ceil(x - 0.5)


def quantize(coeffs, qformat, strategy="symmetric", sec_len=None,
             details=False):
    """Fixed-point conversion (port of fce_quant_core).

    Parameters
    ----------
    coeffs    : sequence of float coefficients (FIR taps, or flat SOS)
    qformat   : "q15" or "q31"
    strategy  : "symmetric" | "section_wise" | "coefficient_wise"
    sec_len   : coefficients per scaling unit for "section_wise"
                (5 for SOS; ignored otherwise). For "symmetric" on SOS,
                pass sec_len = 5*num_sections so the whole set shares a
                scale; for FIR it defaults to the full length.

    Returns
    -------
    q  : list of ints
    scale : float (symmetric) or None
    int_bits : int
    If details=True, also returns (max_abs, rms, max_rel).
    """
    n = len(coeffs)
    if n == 0:
        raise ValueError("empty coefficient set")
    frac = 15 if qformat == "q15" else 31
    qmax = (1 << frac) - 1
    qmax_d = float(qmax)

    if strategy == "symmetric":
        if sec_len is None:
            sec_len = n
        n_sec = (n + sec_len - 1) // sec_len
        mx = max(abs(c) for c in coeffs)
        scale = (qmax_d / mx) if mx > 0.0 else 1.0
        scales = [scale] * n
        sec_scales = [scale] * n_sec
        int_bits = int(math.ceil(math.log2(mx))) if mx > 1.0 else 0
    elif strategy == "section_wise":
        if sec_len is None:
            raise ValueError("section_wise needs sec_len")
        n_sec = (n + sec_len - 1) // sec_len
        scales = [0.0] * n
        sec_scales = []
        int_bits = 0
        scale = None
        for sec in range(n_sec):
            base = sec * sec_len
            cnt = sec_len if (base + sec_len <= n) else (n - base)
            mx = max(abs(coeffs[base + i]) for i in range(cnt)) \
                if cnt > 0 else 0.0
            s = (qmax_d / mx) if mx > 0.0 else 1.0
            for i in range(cnt):
                scales[base + i] = s
            sec_scales.append(s)
    elif strategy == "coefficient_wise":
        scales = []
        mx = 0.0
        sec_scales = None
        scale = None
        for c in coeffs:
            a = abs(c)
            if a > mx:
                mx = a
            scales.append((qmax_d / a) if a > 0.0 else 1.0)
        int_bits = int(math.ceil(math.log2(mx))) if mx > 1.0 else 0
    else:
        raise ValueError("unknown strategy: " + str(strategy))

    # quantize
    q = []
    sse = 0.0
    max_abs = 0.0
    max_rel = 0.0
    for i, c in enumerate(coeffs):
        s = scales[i]
        if s > 0.0:
            qv = _round_half(c * s)
        else:
            qv = 0
        qv = max(-qmax, min(qmax, qv))
        q.append(qv)
        vt = (qv / s) if s > 0.0 else 0.0
        err = abs(c - vt)
        if err > max_abs:
            max_abs = err
        sse += err * err
        if c != 0.0:
            rel = err / abs(c)
            if rel > max_rel:
                max_rel = rel
    rms = math.sqrt(sse / float(n))

    if details:
        return q, scale, sec_scales, int_bits, max_abs, rms, max_rel
    return q, scale, int_bits


# ======================================================================
# top-level design (port of fce_generate's design core)
# ======================================================================

def design(spec, details=False):
    """Design a filter from a spec dict (see module docstring). Returns
    a list of coefficients (FIR: h; IIR: flattened SOS) or, with
    details=True, a rich result dict."""
    s = dict(spec)
    kind = s.get("kind")
    if kind not in ("fir", "iir"):
        raise ValueError("spec['kind'] must be 'fir' or 'iir'")
    fs = s["fs"]
    if not fs > 0.0:
        raise ValueError("fs must be > 0")

    if kind == "fir":
        num_taps = s.get("num_taps", 0)
        taps_clamped = False
        if num_taps == 0:
            if s.get("window") != "kaiser" or not s.get("stopband_atten_db", 0) > 0 \
                    or not s.get("transition_hz", 0) > 0:
                raise ValueError("auto taps require Kaiser + atten + transition")
            num_taps = kaiser_taps(s["stopband_atten_db"], s["transition_hz"], fs)
            if num_taps > MAX_FIR_TAPS:
                num_taps = MAX_FIR_TAPS
                taps_clamped = True
        s["num_taps"] = num_taps
        res = fir_design(s)
        res["taps_clamped"] = taps_clamped
        if details:
            return res
        return res["h"]

    else:  # iir
        order = s.get("order", 0)
        clamped = 0
        if order == 0:
            auto = iir_auto_order(s)
            order = auto["order"]
            clamped = auto.get("clamped", 0)
            s["design_fc1"] = auto["design_fc1"]
            s["design_fc2"] = auto["design_fc2"]
        else:
            if not s.get("design_fc1"):
                s.setdefault("design_fc1", s["fc1"])
            if not s.get("design_fc2") and s.get("fc2"):
                s["design_fc2"] = s["fc2"]
        s["order"] = order
        res = iir_design(s)
        res["order_clamped"] = bool(clamped)
        if details:
            return res
        flat = [v for sec in res["sos"] for v in sec]
        return flat
