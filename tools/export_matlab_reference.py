#!/usr/bin/env python3
"""
export_matlab_reference.py - export reference coefficients from the C
library (tools/reference/fce_dump) as CSV files that the MATLAB validation
script (matlab/validate_filtercoeff.m) can read.

Usage: python3 tools/export_matlab_reference.py [dump_bin]
Writes one CSV per case into matlab/reference/ and a manifest listing the
selected case ids.
"""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "matlab", "reference")

# Representative cases spanning FIR windows and all IIR families.
SELECT = [
    "fir_lp_hann_odd",
    "fir_lp_kaiser",
    "fir_hp_hamming",
    "fir_bp_hamming",
    "fir_bs_kaiser",
    "fir_hilbert",
    "iir_butter_lp4",
    "iir_butter_hp8",
    "iir_butter_bp4",
    "iir_cheby1_lp5",
    "iir_cheby2_lp4",
    "iir_ellip_lp4",
    "iir_ellip_bp6",
    "iir_bessel_lp3",
    "fir_lp_kaiser_auto",
]


def main():
    dump_bin = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(ROOT, "tools", "reference", "fce_dump")
    out = subprocess.run([dump_bin], capture_output=True, text=True)
    if out.returncode != 0:
        print("fce_dump failed:", out.stderr)
        return 1

    os.makedirs(OUT, exist_ok=True)
    manifest = []
    for line in out.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        case = json.loads(line)
        if case["id"] not in SELECT:
            continue
        if case["status"] != 0:
            print("skip", case["id"], "status", case["status"])
            continue
        if case["kind"] == "FIR":
            coeffs = case["coeffs"]
        else:
            coeffs = [v for sec in case["sos"] for v in sec]
        path = os.path.join(OUT, case["id"] + ".csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["index", "coefficient"])
            for i, c in enumerate(coeffs):
                w.writerow([i, c])
        manifest.append(case["id"])
        print("exported", case["id"], "(", len(coeffs), "coeffs )")

    with open(os.path.join(OUT, "manifest.txt"), "w") as f:
        f.write("\n".join(manifest) + "\n")
    print("manifest written:", len(manifest), "cases")


if __name__ == "__main__":
    sys.exit(main())
