#!/usr/bin/env python3
"""Validity check of DNN_RPRSM against enumeration (adapted from
roadmap/certificate/scripts/validity_bruteforce.py).

Random n = 8 QAPs (symmetric and asymmetric flows, integer data, and a linear
term on some), optimum by brute force over all 8! permutations (with the
node's fixed assignment and eliminations); dnntest certifies at every
iteration (-k 1) for 300 iterations with each DNN_RPRSM value and three
eigensolver modes; every certified value must be <= the optimum of that
(sub)problem. Usage: validity_rprsm.py <dnntest> <workdir> [count]
"""
import itertools, os, re, subprocess, sys
import numpy as np

BITS = ["0", "1", "2", "4", "10", "11", "15"]
MODES = [["-d", "3", "-L", "2", "-r", "40"], ["-e", "1"], ["-e", "0"]]


def write_dat(path, A, B, C=None):
    n = len(A)
    with open(path, "w") as f:
        f.write(f"{n}\n\n")
        for M in (A, B) + ((C,) if C is not None else ()):
            for row in M:
                f.write(" ".join(str(int(x)) for x in row) + "\n")
            f.write("\n")


def optimum(A, B, C, fixed=(), elim=()):
    n = len(A)
    best = None
    for p in itertools.permutations(range(n)):
        if any(p[i] != j for i, j in fixed) or any(p[i] == j for i, j in elim):
            continue
        v = sum(A[i][k] * B[p[i]][p[k]] for i in range(n) for k in range(n))
        if C is not None:
            v += sum(C[i][p[i]] for i in range(n))
        best = v if best is None or v < best else best
    return best


def main():
    exe, wd = sys.argv[1], sys.argv[2]
    count = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(wd, exist_ok=True)
    rng = np.random.default_rng(20261105)
    worst, nvals, nrun = -1e300, 0, 0
    for t in range(count):
        n = 8
        A = rng.integers(0, 10, (n, n)); np.fill_diagonal(A, 0)
        if t % 2 == 0:
            A = np.triu(A) + np.triu(A).T
        D = rng.integers(0, 10, (n, n)); np.fill_diagonal(D, 0)
        B = D + D.T if t % 3 else D
        C = rng.integers(0, 30, (n, n)) if t % 4 == 1 else None
        path = os.path.join(wd, f"r{t}.dat")
        write_dat(path, A, B, C)
        cases = (("-", (), "-", ()), ("0:3", ((0, 3),), "-", ()),
                 ("0:3", ((0, 3),), "1:1,2:5,4:4", ((1, 1), (2, 5), (4, 4))))
        for fixes, fx, elims, el in cases:
            opt = optimum(A.tolist(), B.tolist(), None if C is None else C.tolist(), fx, el)
            for mode in MODES:
                for bits in BITS:
                    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                               DNN_RPRSM=bits)
                    out = subprocess.run([exe, "-T", "1", "-g", "1.0", "-G", "1.0", "-H", "0", "-i", "300",
                                          "-k", "1", "-v"] + mode + ["-x", fixes, "-E", elims, path],
                                         env=env, capture_output=True, text=True).stdout
                    vals = [float(x) for x in re.findall(r"\(cert (-?[\d.]+)", out)]
                    if not vals:
                        print(f"r{t} {fixes} {elims} {' '.join(mode)} bits {bits}: no certificate\n{out[-500:]}")
                        continue
                    # the printed checks have 4 decimals; the RESULT lb (the best of
                    # them) has 17 digits and is what is compared
                    mres = re.search(r"^RESULT .*? lb=(\S+)", out, re.M)
                    mx = float(mres.group(1)) if mres else max(vals)
                    nvals += len(vals); nrun += 1
                    worst = max(worst, mx - opt)
                    print(f"r{t} fixes {fixes:4s} elim {elims:11s} {' '.join(mode):16s} bits {bits:2s} opt {opt} "
                          f"max cert {mx:.10f} slack {opt - mx:.3e} ({len(vals)}) "
                          f"{'OK' if mx <= opt + 1e-9 else 'VIOLATION'}", flush=True)
    print(f"{nrun} runs, {nvals} certified values; worst (max cert - optimum) {worst:.3e}: "
          f"{'PASS' if worst <= 1e-9 else 'FAIL'}")


main()
