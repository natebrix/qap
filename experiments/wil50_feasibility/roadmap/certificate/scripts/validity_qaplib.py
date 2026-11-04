#!/usr/bin/env python3
"""Validity of the DNN_CERTX variants on small QAPLIB instances with a DNN gap.

Every variant at every check (DNN_CERTX=47, every check) must stay at or
below the instance's optimum; at the end the largest variant is compared with
the converged standard certificate (variants may exceed it only at
non-converged iterates, never the DNN value by more than its own residual).
Usage: validity_qaplib.py <dnntest> <qaplib dir> <logdir> name:opt ...
"""
import os, re, subprocess, sys

exe, qdir, logdir = sys.argv[1:4]
os.makedirs(logdir, exist_ok=True)
env = dict(os.environ, DNN_CERTX="47", DNN_CERTX_EVERY="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
ok = True
for spec in sys.argv[4:]:
    name, opt = spec.split(":")
    opt = float(opt)
    out = subprocess.run([exe, "-T", "1", "-g", "1.0", "-G", "1.0", "-H", "0", "-i", "1500",
                          os.path.join(qdir, name + ".dat")], env=env, capture_output=True, text=True).stdout
    open(os.path.join(logdir, f"valid_{name}.log"), "w").write(out)
    mx, last = -1e300, None
    for line in out.splitlines():
        if line.startswith("certx it"):
            v = [float(x) for x in re.findall(r"(?:std|exact|ls|avg\d+) (-?[\d.]+)", line)]
            v += [float(x) for x in re.findall(r" \d+:(-?[\d.]+)", line)]
            mx = max(mx, max(v)); last = line
    res = [l for l in out.splitlines() if l.startswith("RESULT")][0]
    lb = float(re.search(r"lb=(\S+)", res).group(1)); pr = float(re.search(r"primal=(\S+)", res).group(1))
    good = mx <= opt
    ok &= good
    print(f"{name}: optimum {opt:.0f} largest variant {mx:.4f} converged standard {lb:.4f} primal {pr:.4f} "
          f"{'OK' if good else 'VIOLATION'}")
print("PASS" if ok else "FAIL")
