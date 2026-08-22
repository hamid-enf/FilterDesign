#!/usr/bin/env python3
"""Concrete input -> output examples for the pure-Python port.

    python3 python/example.py

Shows how to design FIR and IIR filters with filtercoeff.py and how to
convert the result to fixed point.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import filtercoeff as fce


def show_fir():
    # FIR lowpass, Hamming window, 101 taps, Fs = 48 kHz, Fc = 5 kHz
    spec = {
        "kind": "fir", "fir_type": "lowpass",
        "fs": 48000, "fc1": 5000, "num_taps": 101, "window": "hamming",
    }
    res = fce.design(spec, details=True)
    h = res["h"]
    print("=== FIR lowpass (Hamming) ===")
    print(f"  input : Fs={spec['fs']} Hz, Fc={spec['fc1']} Hz, taps={res['num_taps']}")
    print(f"  output: first 8 of {res['num_taps']} taps:")
    print("         " + " ".join(f"{v:.6f}" for v in h[:8]))
    print(f"  sum(h) (DC gain) ~ {sum(h):.4f}\n")


def show_iir():
    # IIR 4th-order Butterworth lowpass, Fs = 48 kHz, Fc = 5 kHz
    spec = {
        "kind": "iir", "iir_family": "butterworth", "iir_type": "lowpass",
        "fs": 48000, "fc1": 5000, "order": 4,
    }
    res = fce.design(spec, details=True)
    print("=== IIR Butterworth lowpass (order 4) ===")
    print(f"  input : Fs={spec['fs']} Hz, Fc={spec['fc1']} Hz, order={res['order']}")
    print(f"  output: {res['num_sections']} SOS sections, layout {{b0,b1,b2,a1,a2}}:")
    for sec in res["sos"]:
        print("         " + " ".join(f"{v:.8f}" for v in sec))
    print()


def show_iir_q15():
    # IIR 6th-order Elliptic bandpass + fixed-point Q15
    spec = {
        "kind": "iir", "iir_family": "elliptic", "iir_type": "bandpass",
        "fs": 48000, "fc1": 2000, "fc2": 4000, "order": 6,
        "passband_ripple_db": 1.0, "stopband_atten_db": 50,
    }
    res = fce.design(spec, details=True)
    flat = [v for sec in res["sos"] for v in sec]
    q, scale, intbits = fce.quantize(flat, "q15", "symmetric", sec_len=len(flat))
    print("=== IIR Elliptic bandpass (order 6) + Q15 ===")
    print(f"  input : Fs=48k, passband [{spec['fc1']},{spec['fc2']}] Hz, "
          f"rp={spec['passband_ripple_db']}, rs={spec['stopband_atten_db']}")
    print(f"  output: {res['num_sections']} sections (float + Q15)")
    print(f"          Q15 scale={scale:.4g}, int_bits={intbits}")
    print(f"          first section Q15: {q[:5]}\n")


def show_auto():
    # auto-order (no order given): Butterworth lowpass, buttord equivalent
    spec = {
        "kind": "iir", "iir_family": "butterworth", "iir_type": "lowpass",
        "fs": 10000, "fc1": 1000, "order": 0,
        "edge1_hz": 2000, "passband_ripple_db": 3.0, "stopband_atten_db": 60,
    }
    res = fce.design(spec, details=True)
    print("=== IIR auto-order (Butterworth lowpass) ===")
    print(f"  input : Fs=10k, passband 0..1000 Hz, stopband from 2000 Hz")
    print(f"  output: computed order = {res['order']}, "
          f"{res['num_sections']} sections (SciPy buttord gives the same)")


if __name__ == "__main__":
    show_fir()
    show_iir()
    show_iir_q15()
    show_auto()
