#!/usr/bin/env python3
"""The best known tho40 solution (QAPSOLN/tho40.sln, 240,516) in tho40t.dat indices,
and node specs (-x i:j,...) that fix its first d assignments in a given order.
  python3 path.py            # verify the value, print the specs used in part 1
"""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../.."))

def read_dat(path):
    t = open(path).read().split()
    n = int(t[0]); v = np.array(t[1:1 + 2 * n * n], float)
    return n, v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)

def solution():
    t = open(os.path.join(ROOT, "QAPSOLN/tho40.sln")).read().split()
    n = int(t[0]); p = [int(x) - 1 for x in t[2:2 + n]]
    q = [0] * n
    for i, x in enumerate(p): q[x] = i          # inverted: tho40t swaps A and B
    return q

def image(q, fr, fc, rows=5, cols=8):
    def g(l):
        r, c = divmod(l, cols)
        return (rows - 1 - r if fr else r) * cols + (cols - 1 - c if fc else c)
    return [g(l) for l in q]

def value(A, B, q):
    q = np.array(q)
    return float((A * B[np.ix_(q, q)]).sum())

WARM = [26, 30, 2, 16, 35, 24, 34, 21, 33, 37, 39, 15]     # path_check warm dive's order
COLD = [26, 30, 39, 2, 16, 35, 24, 34, 33, 37, 21, 15]     # cold dive's order
INDEX = list(range(12))                                     # facilities 0..11 (not a branching order)
REV = list(range(39, 27, -1))                               # facilities 39..28

def spec(q, order, d):
    return ",".join("%d:%d" % (f, q[f]) for f in order[:d])

if __name__ == "__main__":
    n, A, B = read_dat(os.path.join(ROOT, "qaplib/tho40t.dat"))
    q0 = solution()
    for fr in (0, 1):
        for fc in (0, 1):
            q = image(q0, fr, fc)
            print("image rowflip %d colflip %d: value %.0f  facility 26 at %d" % (fr, fc, value(A, B, q), q[26]))
    q = [x for x in (image(q0, fr, fc) for fr in (0, 1) for fc in (0, 1)) if x[26] == 2][0]
    print("path image:", " ".join(map(str, q)))
    for name, order in (("warm", WARM), ("cold", COLD), ("index", INDEX), ("rev", REV)):
        for d in (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12):
            print(name, d, spec(q, order, d))
