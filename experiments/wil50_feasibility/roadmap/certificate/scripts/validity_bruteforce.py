#!/usr/bin/env python3
"""Validity check of the DNN_CERTX certificate variants against enumeration.

Random n = 8 QAPs (symmetric and asymmetric flows, integer data), optimum by
brute force over all 8! permutations; dnntest runs with every variant at every
check (DNN_CERTX=47, DNN_CERTX_EVERY=1) at the root and under a fixed node
(-x 0:3); every printed value must be <= the optimum of that (sub)problem.
Usage: validity_bruteforce.py <dnntest> <workdir> [count]
"""
import itertools, os, re, subprocess, sys
import numpy as np


def write_dat(path, A, B):
    n = len(A)
    with open(path, "w") as f:
        f.write(f"{n}\n\n")
        for M in (A, B):
            for row in M:
                f.write(" ".join(str(int(x)) for x in row) + "\n")
            f.write("\n")


def optimum(A, B, fixed=None):
    n = len(A)
    best = None
    for p in itertools.permutations(range(n)):
        if fixed and any(p[i] != j for i, j in fixed):
            continue
        v = sum(A[i][k] * B[p[i]][p[k]] for i in range(n) for k in range(n))
        best = v if best is None or v < best else best
    return best


def main():
    exe, wd = sys.argv[1], sys.argv[2]
    count = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(wd, exist_ok=True)
    rng = np.random.default_rng(20261104)
    env = dict(os.environ, DNN_CERTX="47", DNN_CERTX_EVERY="1", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1")
    worst = -1e300
    for t in range(count):
        n = 8
        A = rng.integers(0, 10, (n, n)); np.fill_diagonal(A, 0)
        if t % 2 == 0:
            A = np.triu(A) + np.triu(A).T
        D = rng.integers(0, 10, (n, n)); np.fill_diagonal(D, 0)
        B = D + D.T if t % 3 else D
        path = os.path.join(wd, f"r{t}.dat")
        write_dat(path, A, B)
        for fixes, fx in (("-", None), ("0:3", [(0, 3)])):
            opt = optimum(A.tolist(), B.tolist(), fx)
            out = subprocess.run([exe, "-T", "1", "-g", "1.0", "-G", "1.0", "-H", "0", "-i", "400",
                                  "-x", fixes, path], env=env, capture_output=True, text=True).stdout
            vals = []
            for line in out.splitlines():
                if line.startswith("certx it"):
                    vals += [float(x) for x in re.findall(r"(?:std|exact|ls|avg\d+|proj) (-?[\d.]+)", line)]
                    vals += [float(x) for x in re.findall(r" \d+:(-?[\d.]+)", line)]
            mx = max(vals)
            worst = max(worst, mx - opt)
            print(f"r{t} fixes {fixes:4s} opt {opt} max variant {mx:.6f} slack {opt - mx:.6f} "
                  f"({len(vals)} values) {'OK' if mx <= opt + 1e-9 else 'VIOLATION'}")
    print(f"worst (max variant - optimum) {worst:.6f}: {'PASS' if worst <= 1e-9 else 'FAIL'}")


main()
