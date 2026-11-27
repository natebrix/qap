#!/usr/bin/env python3
"""The step of a QAPLIB instance's objective grid, decided as the solver decides its
fathoming tolerance (latest/qap.C, set_fathom_tolerance, on the data as readQAPLIB
reads it): 2 when every objective value is an even integer, else 1. start.sh gives a
new run the best known + this step: the next value above the best known, so the run
must find the best known itself and prove that nothing is better.

  python3 infra/phase1/objstep.py qaplib/tho40t.dat
  2 even integers (the solver's fathoming tolerance 1.95)

Integers that are not all even: 1 (tolerance 0.95). Values that are not integers, or
too large to compute exactly in double: there is no grid (tolerance 0), and 1 is a
margin, not a step.
"""
import math
import sys


def read(path):
    """n, A, B and C as readQAPLIB reads them: C is the next n * n numbers if the file has any"""
    t = open(path).read().split()
    try:
        n = int(t[0])
        v = [float(x) for x in t[1:]]
    except (ValueError, IndexError):
        sys.exit(f"{path}: not a QAPLIB data file (n, then A and B)")
    if len(v) < 2 * n * n:
        sys.exit(f"{path}: {len(v)} numbers after n = {n}, not the {2 * n * n} of A and B")
    C = v[2 * n * n:3 * n * n]
    return n, v[:n * n], v[n * n:2 * n * n], C + [0.0] * (n * n - len(C))


def step(path):
    n, A, B, C = read(path)
    exact = sum(map(abs, A)) * max(map(abs, B)) + sum(map(abs, C)) < 2.0 ** 53
    integer = exact and all(x == math.floor(x) for M in (A, B, C) for x in M)
    sym = lambda M: all(M[i * n + j] == M[j * n + i] for i in range(n) for j in range(i + 1, n))
    even = integer and sym(A) and sym(B) and all(
        math.fmod(abs(A[i * n + i] * B[k * n + k]), 2.0) == 0 and math.fmod(abs(C[i * n + k]), 2.0) == 0
        for i in range(n) for k in range(n))
    if even:
        return 2, "even integers (the solver's fathoming tolerance 1.95)"
    if integer:
        return 1, "integers, not all even (the solver's fathoming tolerance 0.95)"
    return 1, ("not integers" if exact else "too large to compute exactly") + " (no grid: the solver's fathoming tolerance 0)"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    s, what = step(sys.argv[1])
    print(s, what)
