#!/usr/bin/env python3
"""Find small instances/nodes with a DNN gap: brute force + DNN (SCS)."""
import sys, json, itertools
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from relax_small import *
out = []
def run(name, A, B, C, const, rows, cols, free0):
    m = len(A)
    v, p = brute(m, A, B, C)
    L = lifted_cost(m, A, B, C)
    d, st, t, _ = solve(m, L, ["dnn"], None, solver="SCS")
    rec = dict(name=name, m=m, opt=v + const, dnn=d + const, gap=v - d, rel=(v - d) / (v + const))
    print(json.dumps(rec), flush=True)
    return rec
mode = sys.argv[1]
if mode == "nug12":
    n, A, B, C = read_qaplib('/home/user/qap/qaplib/nug12.dat')
    # nug12: first matrix is the flow, the second the 3x4 grid distance: put
    # the grid on side 0 by swapping (A <-> B transposes the assignment)
    A, B = B, A
    rng = np.random.default_rng(int(sys.argv[2]))
    for t in range(int(sys.argv[3])):
        k = int(sys.argv[4])
        I = rng.choice(12, k, replace=False); J = rng.choice(12, k, replace=False)
        fixes = list(zip(I.tolist(), J.tolist()))
        Ar, Br, Cr, const, f0, f1 = reduce_instance(A, B, None, fixes)
        run("nug12:" + ",".join("%d:%d" % f for f in fixes), Ar, Br, Cr, const, 3, 4, f0)
