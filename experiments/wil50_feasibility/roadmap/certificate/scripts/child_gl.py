#!/usr/bin/env python3
"""Per-child bounds from a parent's certificate matrix (question 3).

Reads <prefix>.M.bin (M = L + scale Zc, the parent's lifted (m^2+1)^2 matrix,
row-major, unscaled) and <prefix>.meta (written by dnn.C under DNN_CERTX bit
16). For a child fixing facility a to location b (local indices), every
lifted permutation with sigma(a) = b satisfies

  <M, yy'> = M00 + lin_ab + sum_{k != a} [lin_{k s_k} + M_(k s_k),(ab) + M_(ab),(k s_k)
                                          + sum_{k' != a, k} M_(k s_k),(k' s_k')]

so the child's Gilmore-Lawler value with the row and column fixed exactly is
M00 + lin_ab + LAP_{m-1}(c^ab), c^ab_kl = lin_kl + M_kl,ab + M_ab,kl +
LAP_{m-2}(M_(kl),. without facilities a, k and locations b, l). The parent's
semidefinite term (m+1) scale lambda and margin are the same (the child's
permutations are the parent's), so child bound = value - margin + shift.
It is never below the reduced-cost bound GL + r_ab (dnn.C child_lb).
Prints the parent GL (a cross-check of dnn.C's exact value), and per child the
reduced-cost bound and the exact per-child bound; counts those above the cutoff.
Usage: child_gl.py <prefix> [all | line:row:<i> | line:col:<j>]
"""
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment


def lap(c):
    r, s = linear_sum_assignment(c)
    return c[r, s].sum()


def main():
    pre = sys.argv[1]
    which = sys.argv[2] if len(sys.argv) > 2 else "all"
    meta = {}
    for line in open(pre + ".meta"):
        k, v = line.split()
        meta[k] = float(v)
    m = int(meta["m"]); n1 = m * m + 1
    M = np.fromfile(pre + ".M.bin", dtype=np.float64).reshape(n1, n1)
    T = M[1:, 1:].reshape(m, m, m, m)
    lin = (M[0, 1:] + M[1:, 0] + np.diag(M)[1:]).reshape(m, m)
    M00 = M[0, 0]
    margin, shift, cutoff = meta["margin"], meta["shift"], meta["cutoff"]
    allidx = np.arange(m)
    # parent GL with exact inner LAPs
    C = np.empty((m, m))
    for i in range(m):
        ri = allidx[allidx != i]
        for j in range(m):
            cj = allidx[allidx != j]
            C[i, j] = lin[i, j] + lap(T[i, j][np.ix_(ri, cj)])
    gl = M00 + lap(C)
    print(f"parent m {m} it {int(meta['it'])} GL exact {gl - margin + shift:.4f} "
          f"(dnn.C exact {meta['exact'] + shift:.4f}) cutoff {cutoff:.4f} short {cutoff - (gl - margin + shift):.3f}")
    # exact restricted outer LAP: the best reduced-cost-type child bound
    if which == "all":
        kids = [(a, b) for a in range(m) for b in range(m)]
    else:
        _, kind, v = which.split(":")
        v = int(v)
        kids = [(v, b) for b in range(m)] if kind == "row" else [(a, v) for a in range(m)]
    over_red = over_exact = 0
    gains = []
    for a, b in kids:
        ra = allidx[allidx != a]; cb = allidx[allidx != b]
        red = M00 + C[a, b] + lap(C[np.ix_(ra, cb)])          # GL with the outer row fixed
        Cc = np.empty((m - 1, m - 1))
        for x, k in enumerate(ra):
            rk = allidx[(allidx != a) & (allidx != k)]
            for y, l in enumerate(cb):
                cl = allidx[(allidx != b) & (allidx != l)]
                Cc[x, y] = lin[k, l] + T[k, l, a, b] + T[a, b, k, l] + lap(T[k, l][np.ix_(rk, cl)])
        ex = M00 + lin[a, b] + lap(Cc)
        vr, ve = red - margin + shift, ex - 2 * margin + shift
        over_red += vr > cutoff; over_exact += ve > cutoff
        gains.append((ve - vr, ve, vr, a, b))
        print(f"child {a}:{b} outer-fixed {vr:.4f} exact-child {ve:.4f} gain {ve - vr:.4f}"
              f" {'FATHOM' if ve > cutoff else ''}", flush=True)
    g = np.array([x[0] for x in gains]); ve = np.array([x[1] for x in gains])
    print(f"children {len(kids)}: fathomed outer-fixed {over_red}, exact-child {over_exact}; "
          f"exact-child over parent GL: median {np.median(ve - (gl - margin + shift)):.3f} "
          f"max {np.max(ve - (gl - margin + shift)):.3f}; gain over outer-fixed median {np.median(g):.3f} max {g.max():.3f}")


main()
