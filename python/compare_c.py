#!/usr/bin/env python3
"""
compare_c.py - validate the pure-Python port against the C library.

Builds the C dump tool (tools/reference/fce_dump), which designs a fixed
matrix of FIR/IIR filters, then re-designs the exact same filters with the
Python port and compares the float64 coefficients.

    C (FilterCoeff)  <--->  python/filtercoeff.py

Usage:  python3 compare_c.py [path-to-fce_dump]
Exit 0 if every design matches to within the tolerance.
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import filtercoeff as fce  # noqa: E402

TOL = 1e-12  # max abs coefficient difference (both are float64)

NORM_MAP = {"DC": "dc", "Nyquist": "nyquist", "Peak": "passband_peak",
            "None": "none"}

WINDOW_MAP = {"Rectangular": "rectangular", "Hann": "hann",
              "Hamming": "hamming", "Blackman": "blackman",
              "Kaiser": "kaiser", "Blackman-Harris": "blackman_harris",
              "Bartlett": "bartlett", "Tukey": "tukey"}

FAMILY_MAP = {"Butterworth": "butterworth", "Chebyshev I": "chebyshev1",
              "Chebyshev II": "chebyshev2", "Elliptic": "elliptic",
              "Bessel": "bessel"}

FIR_TYPE_MAP = {"Lowpass": "lowpass", "Highpass": "highpass",
                "Bandpass": "bandpass", "Bandstop": "bandstop",
                "Hilbert": "hilbert", "Differentiator": "differentiator"}

IIR_TYPE_MAP = {"Lowpass": "lowpass", "Highpass": "highpass",
                "Bandpass": "bandpass", "Bandstop": "bandstop"}

FAILURES = []
TOTAL = [0]


def check(name, ok, detail=""):
    TOTAL[0] += 1
    if not ok:
        FAILURES.append(f"{name}: {detail}")


def fir_spec(case):
    sp = {
        "kind": "fir",
        "fs": case["fs"],
        "fc1": case["fc1"],
        "fc2": case.get("fc2", 0.0),
        "fir_type": FIR_TYPE_MAP[case["type"]],
        "window": WINDOW_MAP[case["window"]],
        "num_taps": case["taps"],
        "kaiser_beta": case["beta"],
        "normalization": NORM_MAP[case["norm"]],
    }
    return sp


def iir_spec(case, use_design_fc=True):
    sp = {
        "kind": "iir",
        "fs": case["fs"],
        "fc1": case["fc1"],
        "fc2": case.get("fc2", 0.0),
        "iir_family": FAMILY_MAP[case["family"]],
        "iir_type": IIR_TYPE_MAP[case["type"]],
        "order": case["order"],
        "passband_ripple_db": case["rp"],
        "stopband_atten_db": case["rs"],
    }
    if use_design_fc:
        sp["design_fc1"] = case.get("design_fc1", case["fc1"])
        sp["design_fc2"] = case.get("design_fc2", case.get("fc2", 0.0))
    return sp


def run(dump_bin):
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
            check(cid, False, f"C generation failed status={case['status']}")
            continue

        if case["kind"] == "FIR":
            # auto-taps Kaiser case: exercise the auto-tap path
            if cid == "fir_lp_kaiser_auto":
                sp = {"kind": "fir", "fs": case["fs"], "fc1": case["fc1"],
                      "fir_type": "lowpass", "window": "kaiser",
                      "num_taps": 0, "stopband_atten_db": 80,
                      "transition_hz": 1000, "normalization": "dc"}
                res = fce.design(sp, details=True)
                check(cid, res["num_taps"] == case["taps"],
                      f"auto taps {res['num_taps']} vs C {case['taps']}")
                h = res["h"]
            else:
                sp = fir_spec(case)
                res = fce.design(sp, details=True)
                h = res["h"]
            c = case["coeffs"]
            if len(h) != len(c):
                check(cid, False, f"len {len(h)} vs C {len(c)}")
                continue
            err = max(abs(a - b) for a, b in zip(h, c))
            check(cid, err < TOL, f"FIR max coeff err {err:.3e}")
            check(cid, abs(res["kaiser_beta"] - case["beta"]) < 1e-12,
                  f"beta {res['kaiser_beta']} vs C {case['beta']}")
            # fixed-point (fir_lp_q15: symmetric Q15)
            if "q15" in case:
                q, scale, _ = fce.quantize(h, "q15", "symmetric",
                                           sec_len=len(h))
                check(cid, list(q) == case["q15"],
                      f"Q15 mismatch max={max(abs(a-b) for a,b in zip(q, case['q15']))}")
                check(cid, abs(scale - case["scale"])
                      / max(1.0, abs(case["scale"])) < 1e-9,
                      f"Q15 scale {scale} vs C {case['scale']}")

        else:  # IIR
            if cid in ("iir_butter_auto_lp", "iir_cheby1_auto_hp",
                       "iir_cheby2_auto_bp", "iir_ellip_auto_bs"):
                # validate the auto-order path first (edge frequencies are
                # not in the dump JSON, so they are hardcoded to match
                # fce_dump.c)
                auto_spec = {
                    "iir_butter_auto_lp": {"fs": 10000, "fc1": 1000,
                                           "edge1_hz": 2000, "rp": 3.0,
                                           "rs": 60, "fam": "butterworth",
                                           "typ": "lowpass"},
                    "iir_cheby1_auto_hp": {"fs": 10000, "fc1": 2000,
                                           "edge1_hz": 1000, "rp": 1.0,
                                           "rs": 50, "fam": "chebyshev1",
                                           "typ": "highpass"},
                    "iir_cheby2_auto_bp": {"fs": 48000, "fc1": 2000,
                                           "fc2": 4000, "edge1_hz": 1500,
                                           "edge2_hz": 5000, "rp": 1.0,
                                           "rs": 40, "fam": "chebyshev2",
                                           "typ": "bandpass"},
                    "iir_ellip_auto_bs": {"fs": 48000, "fc1": 3000,
                                          "fc2": 5000, "edge1_hz": 2000,
                                          "edge2_hz": 7000, "rp": 0.5,
                                          "rs": 60, "fam": "elliptic",
                                          "typ": "bandstop"},
                }[cid]
                asp = {"kind": "iir", "fs": auto_spec["fs"],
                       "fc1": auto_spec["fc1"], "fc2": auto_spec.get("fc2", 0),
                       "iir_family": auto_spec["fam"],
                       "iir_type": auto_spec["typ"], "order": 0,
                       "passband_ripple_db": auto_spec["rp"],
                       "stopband_atten_db": auto_spec["rs"],
                       "edge1_hz": auto_spec["edge1_hz"],
                       "edge2_hz": auto_spec.get("edge2_hz", 0)}
                auto = fce.iir_auto_order(asp)
                check(cid, auto["order"] == case["order"],
                      f"auto order {auto['order']} vs C {case['order']}")
                check(cid, abs(auto["design_fc1"] - case["design_fc1"]) < 1e-9,
                      f"auto fc1 {auto['design_fc1']} vs C {case['design_fc1']}")
                # then compare the SOS using the C design frequencies
                sp = iir_spec(case, use_design_fc=True)
            else:
                sp = iir_spec(case, use_design_fc=True)
            res = fce.design(sp, details=True)
            c = case["sos"]
            if len(res["sos"]) != len(c):
                check(cid, False, f"sections {len(res['sos'])} vs C {len(c)}")
                continue
            err = max(abs(v - cc) for sec, csec in zip(res["sos"], c)
                      for v, cc in zip(sec, csec))
            check(cid, err < TOL, f"IIR max coeff err {err:.3e}")
            flat = [v for sec in res["sos"] for v in sec]
            if "q31" in case and cid == "iir_butter_lp4_q31_secwise":
                q, _, sec_scales, ib, *_ = fce.quantize(
                    flat, "q31", "section_wise", sec_len=5, details=True)
                check(cid, list(q) == case["q31"],
                      f"Q31 sec-wise mismatch "
                      f"max={max(abs(a-b) for a,b in zip(q, case['q31']))}")
            if "q15" in case and cid == "iir_ellip_lp4_q15":
                q, scale, _ = fce.quantize(flat, "q15", "symmetric",
                                           sec_len=len(flat))
                check(cid, list(q) == case["q15"],
                      f"Q15 SOS mismatch "
                      f"max={max(abs(a-b) for a,b in zip(q, case['q15']))}")
                check(cid, abs(scale - case["scale"])
                      / max(1.0, abs(case["scale"])) < 1e-9,
                      f"Q15 scale {scale} vs C {case['scale']}")

    print(f"C vs Python port: {TOTAL[0]} checks, {len(FAILURES)} failures")
    for f in FAILURES:
        print("  FAIL:", f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    dump_bin = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(__file__), "..", "tools",
                     "reference", "fce_dump")
    sys.exit(run(dump_bin))
