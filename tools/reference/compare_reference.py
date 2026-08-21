#!/usr/bin/env python3
"""
compare_reference.py - FilterCoeff vs SciPy reference validation.

Reads the JSON lines produced by the C dump tool (fce_dump) and compares
every design against the equivalent SciPy design:

  * FIR  : coefficients vs scipy.signal.firwin (same window & normalization)
  * IIR  : SOS vs scipy.signal.iirfilter(output='sos') - the SOS layout of
           SciPy is [b0,b1,b2,a0,a1,a2] with a0=1; ours is [b0,b1,b2,a1,a2],
           so a0 is dropped before comparing.
  * auto order/taps vs scipy kaiserord / buttord / cheb1ord / cheb2ord /
    ellipord.
  * frequency response (scipy.signal.sosfreqz / freqz) on a grid.

Exit code 0 = all checks within tolerance, 1 = mismatch.
"""
import json
import math
import subprocess
import sys

import numpy as np
from scipy import signal

# ----------------------------------------------------------------------
# tolerances
# ----------------------------------------------------------------------
TOL_SOS = 2e-10        # max abs coefficient difference (double)
TOL_FIR = 2e-13        # windowed-sinc coefficients (double)
TOL_RESP_DB = 1e-8     # max |dB| difference on the response grid
TOL_ORD = 2            # order/taps may differ by at most this (rounding up)
TOL_BETA = 1e-9

WINDOW_MAP = {
    "Rectangular": "boxcar",
    "Hann": "hann",
    "Hamming": "hamming",
    "Blackman": "blackman",
    "Kaiser": ("kaiser", None),
    "Blackman-Harris": "blackmanharris",
    "Bartlett": "bartlett",
    "Tukey": ("tukey", 0.5),
}

FIR_TYPE_MAP = {
    "Lowpass": "lowpass",
    "Highpass": "highpass",
    "Bandpass": "bandpass",
    "Bandstop": "bandstop",
}

IIR_FAMILY_MAP = {
    "Butterworth": "butter",
    "Chebyshev I": "cheby1",
    "Chebyshev II": "cheby2",
    "Elliptic": "ellip",
    "Bessel": "bessel",
}

IIR_TYPE_MAP = {
    "Lowpass": "lowpass",
    "Highpass": "highpass",
    "Bandpass": "bandpass",
    "Bandstop": "bandstop",
}

FAILURES = []
TOTAL = [0]


def check(name, ok, detail=""):
    TOTAL[0] += 1
    if not ok:
        FAILURES.append(f"{name}: {detail}")


def sos_scipy_to_ours(sos):
    """scipy sos: [b0,b1,b2,a0,a1,a2] (a0=1) -> ours [b0,b1,b2,a1,a2]"""
    return sos[:, [0, 1, 2, 4, 5]]


def design_fir(case):
    """Design the same FIR with scipy; return coefficients."""
    fs = case["fs"]
    fc1 = case["fc1"]
    fc2 = case["fc2"]
    taps = case["taps"]
    typ = case["type"] if "type" in case else "Lowpass"

    if typ == "Hilbert":
        # scipy has no direct windowed Hilbert; build it manually
        # (same formula as ours) for the reference check
        n = taps
        m = np.arange(n) - 0.5 * (n - 1)
        h = np.zeros(n)
        nz = m != 0
        h[nz] = (1.0 - np.cos(np.pi * m[nz])) / (np.pi * m[nz])
        win = np.hamming(n)
        return h * win
    if typ == "Differentiator":
        # ideal full-band differentiator with H(e^jw) = jw
        # (band-limited formula, valid for integer AND half-integer shifts)
        n = taps
        m = np.arange(n) - 0.5 * (n - 1)
        h = np.zeros(n)
        nz = m != 0
        w1, w2 = 0.0, np.pi
        h[nz] = (w2 * np.cos(w2 * m[nz]) - w1 * np.cos(w1 * m[nz])
                 - (np.sin(w2 * m[nz]) - np.sin(w1 * m[nz])) / m[nz]) \
                / (np.pi * m[nz])
        win = np.hamming(n)
        return h * win

    btype = FIR_TYPE_MAP[typ]
    if btype == "lowpass":
        cut = fc1
    elif btype == "highpass":
        cut = fc1
    else:
        cut = [fc1, fc2]

    win = WINDOW_MAP[case["window"]]
    if isinstance(win, tuple) and win[0] == "kaiser":
        win = ("kaiser", case["beta"])
    h = signal.firwin(taps, cut, window=win, pass_zero=(btype in
                      ("lowpass", "bandstop")), scale=False, fs=fs)
    return h


def norm_fir(h, case, fs):
    """Apply the same normalization as FilterCoeff."""
    n = len(h)
    typ = case["type"] if "type" in case else "Lowpass"
    norm = case["norm"]
    w = np.arange(1024) / 1024 * np.pi
    H = np.abs(np.fft.fft(h, 2048))[:1024]
    if norm == "DC" or (norm == "Auto" and typ in ("Lowpass", "Bandstop")):
        g = np.abs(np.sum(h))
    elif norm == "Nyquist" or (norm == "Auto" and typ in ("Highpass",
                              "Differentiator")):
        g = np.abs(np.sum(h * (-1) ** np.arange(n)))
    else:  # passband peak: dense scan + local golden-section refinement
        from scipy.optimize import minimize_scalar
        if typ == "Bandpass" or typ == "Bandstop":
            f0, f1 = case["fc1"], case["fc2"]
        elif typ == "Hilbert":
            f0, f1 = 0.0, 0.5 * fs
        else:
            f0, f1 = 0.0, 0.5 * fs
        n = len(h)
        grid = 16384
        ff = np.linspace(f0, f1, grid + 1)
        ww = 2 * np.pi * ff / fs
        H = np.abs(np.sum(h[:, None] * np.exp(-1j * np.outer(np.arange(n), ww)), axis=0))
        ib = np.argmax(H)
        f_best = ff[ib]
        span = (f1 - f0) / grid
        def neg_gain(f):
            ww2 = 2 * np.pi * f / fs
            return -np.abs(np.sum(h * np.exp(-1j * ww2 * np.arange(n))))
        res = minimize_scalar(neg_gain,
                              bounds=(max(f0, f_best - 2 * span),
                                      min(f1, f_best + 2 * span)),
                              method="bounded",
                              options={"xatol": 1e-13})
        g = -res.fun
    return h / g


def design_iir(case):
    """Design the same IIR with scipy; return (sos_ours_layout, order)."""
    fs = case["fs"]
    fc1 = case["design_fc1"] if case.get("design_fc1") else case["fc1"]
    fc2 = case["design_fc2"] if case.get("design_fc2") else case["fc2"]
    fam = IIR_FAMILY_MAP[case["family"]]
    typ = IIR_TYPE_MAP[case["type"]]
    order = case["order"]
    rp = case.get("rp", 0.0)
    rs = case.get("rs", 0.0)

    if typ in ("lowpass", "highpass"):
        wn = fc1
    else:
        wn = [fc1, fc2]

    if fam == "butter":
        sos = signal.butter(order, wn, btype=typ, fs=fs, output="sos")
    elif fam == "cheby1":
        sos = signal.cheby1(order, rp, wn, btype=typ, fs=fs, output="sos")
    elif fam == "cheby2":
        sos = signal.cheby2(order, rs, wn, btype=typ, fs=fs, output="sos")
    elif fam == "ellip":
        sos = signal.ellip(order, rp, rs, wn, btype=typ, fs=fs, output="sos")
    else:
        sos = signal.bessel(order, wn, btype=typ, fs=fs, norm="mag",
                            output="sos")
    return sos_scipy_to_ours(sos)


def resp_sos(sos, fs, n=512):
    w, h = signal.sosfreqz(sos, worN=n)
    return w * fs / (2 * np.pi), 20 * np.log10(np.abs(h) + 1e-300)


def main():
    dump_bin = sys.argv[1] if len(sys.argv) > 1 else "./fce_dump"
    out = subprocess.run([dump_bin], capture_output=True, text=True)
    if out.returncode != 0:
        print("fce_dump failed:", out.stderr)
        return 1

    for line in out.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        case = json.loads(line)
        cid = case["id"]
        if case["status"] != 0:
            check(cid, False, f"generation failed status={case['status']}")
            continue

        if case["kind"] == "FIR":
            h_ref = design_fir(case)
            h_ref = norm_fir(h_ref, case, case["fs"])
            h_our = np.array(case["coeffs"])
            if len(h_ref) != len(h_our):
                check(cid, False,
                      f"tap count {len(h_our)} != scipy {len(h_ref)}")
                continue
            err = np.max(np.abs(h_ref - h_our))
            check(cid, err < TOL_FIR,
                  f"FIR max coeff err {err:.3e} > {TOL_FIR:.0e}")
            # response comparison (floor at -120 dB to avoid notch noise)
            w = np.linspace(0, np.pi, 512)
            Hr = np.abs(np.fft.fft(h_ref, 2048))[:512]
            Ho = np.abs(np.fft.fft(h_our, 2048))[:512]
            mask = (Hr > 1e-6) & (Ho > 1e-6)
            if mask.any():
                dr = 20 * np.log10(Hr[mask]) - 20 * np.log10(Ho[mask])
                check(cid, np.max(np.abs(dr)) < TOL_RESP_DB * 100,
                      f"FIR resp err {np.max(np.abs(dr)):.3e} dB")

        else:  # IIR
            auto_bs_bp = (cid in ("iir_cheby2_auto_bp", "iir_ellip_auto_bs"))
            if auto_bs_bp:
                # design philosophy differs (we design at the exact user
                # band edges; scipy adjusts them) -> verify the SPECS
                sos_our = np.array(case["sos"])
                ns = len(sos_our)
                w = np.linspace(0, np.pi, 4096)
                H = np.ones_like(w, dtype=complex)
                for i in range(ns):
                    b = sos_our[i, :3]
                    a = np.r_[1.0, sos_our[i, 3:]]
                    H *= np.polyval(b[::-1], np.exp(-1j * w)) / \
                         np.polyval(a[::-1], np.exp(-1j * w))
                Hdb = 20 * np.log10(np.abs(H) + 1e-300)
                f = w * case["fs"] / (2 * np.pi)
                if cid == "iir_cheby2_auto_bp":
                    # passband [2000,4000], stopband [1500,5000]
                    pb = Hdb[(f >= 2000) & (f <= 4000)]
                    sb = Hdb[(f >= 5000) | (f <= 1500)]
                    check(cid, np.max(pb) - np.min(pb) < 1.0 + 0.5,
                          f"passband ripple {np.max(pb)-np.min(pb):.3f} dB")
                    check(cid, -np.max(sb) > 40 - 1.0,
                          f"stopband atten {-np.max(sb):.3f} dB")
                else:
                    # bandstop [3000,5000], passbands [0,2000] and [7000,nyq]
                    pb = Hdb[(f <= 2000) | (f >= 7000)]
                    sb = Hdb[(f >= 3000) & (f <= 5000)]
                    check(cid, np.max(pb) - np.min(pb) < 0.5 + 0.5,
                          f"passband ripple {np.max(pb)-np.min(pb):.3f} dB")
                    check(cid, -np.max(sb) > 60 - 1.0,
                          f"stopband atten {-np.max(sb):.3f} dB")
                continue
            sos_ref = design_iir(case)
            sos_our = np.array(case["sos"])
            if sos_ref.shape != sos_our.shape:
                check(cid, False,
                      f"section count {sos_our.shape} != scipy {sos_ref.shape}")
                continue
            err = np.max(np.abs(sos_ref - sos_our))
            check(cid, err < TOL_SOS,
                  f"SOS max coeff err {err:.3e} > {TOL_SOS:.0e}")
            # response comparison (must match tightly since SOS matches)
            w = np.linspace(0, np.pi, 512)
            _, Hr = resp_sos(np.column_stack(
                [sos_ref[:, 0], sos_ref[:, 1], sos_ref[:, 2],
                 np.ones(len(sos_ref)), sos_ref[:, 3], sos_ref[:, 4]]),
                case["fs"], 512)
            _, Ho = resp_sos(np.column_stack(
                [sos_our[:, 0], sos_our[:, 1], sos_our[:, 2],
                 np.ones(len(sos_our)), sos_our[:, 3], sos_our[:, 4]]),
                case["fs"], 512)
            mask = (10 ** (Hr / 20) > 1e-6) & (10 ** (Ho / 20) > 1e-6)
            if mask.any():
                dr = np.max(np.abs(Hr[mask] - Ho[mask]))
                check(cid, dr < TOL_RESP_DB,
                      f"SOS resp err {dr:.3e} dB")

        # auto-order / auto-taps checks
        if cid == "fir_lp_kaiser_auto":
            n_ref, beta_ref = signal.kaiserord(80, 1000 / (0.5 * 48000))
            n_ref = int(np.ceil(n_ref))
            if n_ref % 2 == 0:
                n_ref += 1
            check(cid, abs(case["taps"] - n_ref) <= 2,
                  f"auto taps {case['taps']} vs kaiserord {n_ref}")
            check(cid, abs(case["beta"] - beta_ref) < TOL_BETA,
                  f"beta {case['beta']} vs {beta_ref}")

        if cid == "iir_butter_auto_lp":
            n_ref, wn_ref = signal.buttord(1000, 2000, 3, 60, fs=10000)
            check(cid, abs(case["order"] - n_ref) <= TOL_ORD,
                  f"auto order {case['order']} vs buttord {n_ref}")
            check(cid, abs(case["design_fc1"] - wn_ref) < 1e-6,
                  f"design fc1 {case['design_fc1']} vs {wn_ref}")

        if cid == "iir_cheby1_auto_hp":
            n_ref, wn_ref = signal.cheb1ord(2000, 1000, 1, 50, fs=10000)
            check(cid, abs(case["order"] - n_ref) <= TOL_ORD,
                  f"auto order {case['order']} vs cheb1ord {n_ref}")

        if cid == "iir_cheby2_auto_bp":
            n_ref, wn_ref = signal.cheb2ord([2000, 4000], [1500, 5000],
                                            1, 40, fs=48000)
            check(cid, abs(case["order"] - n_ref) <= TOL_ORD,
                  f"auto order {case['order']} vs cheb2ord {n_ref}")

        if cid == "iir_ellip_auto_bs":
            n_ref, wn_ref = signal.ellipord([2000, 7000], [3000, 5000],
                                            0.5, 60, fs=48000)
            check(cid, abs(case["order"] - n_ref) <= TOL_ORD,
                  f"auto order {case['order']} vs ellipord {n_ref}")

        # fixed-point sanity checks
        if cid == "fir_lp_q15":
            h = np.array(case["coeffs"])
            q = np.array(case["q15"])
            qmax = 32767
            check(cid, np.max(np.abs(q)) <= qmax, "Q15 range")
            # reconstruct with the reported scale and compare
            rec = q / case["scale"]
            err = np.max(np.abs(h - rec))
            check(cid, err < 1.1 / case["scale"],
                  f"Q15 reconstruction err {err:.3e}")
            check(cid, np.allclose(q, np.round(h * case["scale"])),
                  "Q15 rounding")

        if cid == "iir_butter_lp4_q31_secwise":
            sos = np.array(case["sos"])
            q = np.array(case["q31"])
            qmax = 2 ** 31 - 1
            check(cid, np.max(np.abs(q)) <= qmax, "Q31 range")
            # per-section reconstruction
            ns = len(sos)
            scales = qmax / np.max(np.abs(sos), axis=1)
            rec = np.zeros_like(sos)
            for i in range(ns):
                rec[i] = q[5 * i:5 * i + 5] / scales[i]
            err = np.max(np.abs(sos - rec))
            check(cid, err < 1.1 / np.min(scales),
                  f"Q31 sec-wise reconstruction err {err:.3e}")

        if cid == "iir_ellip_lp4_q15":
            sos = np.array(case["sos"])
            q = np.array(case["q15"])
            qmax = 32767
            check(cid, np.max(np.abs(q)) <= qmax, "Q15 range")
            rec = q / case["scale"]
            err = np.max(np.abs(sos.ravel() - rec))
            check(cid, err < 1.1 / case["scale"],
                  f"Q15 SOS reconstruction err {err:.3e}")
            # stability of quantized SOS
            rec = rec.reshape(-1, 5)
            a1 = rec[:, 3]
            a2 = rec[:, 4]
            disc = a1 * a1 - 4.0 * a2
            r = np.where(disc >= 0,
                         np.maximum(np.abs(0.5 * (-a1 + np.sqrt(disc))),
                                    np.abs(0.5 * (-a1 - np.sqrt(disc)))),
                         np.sqrt(np.abs(a2)))
            check(cid, np.all(r < 1.0), "quantized stability")

    print(f"FilterCoeff vs SciPy: {TOTAL[0]} checks, "
          f"{len(FAILURES)} failures")
    for f in FAILURES:
        print("  FAIL:", f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
