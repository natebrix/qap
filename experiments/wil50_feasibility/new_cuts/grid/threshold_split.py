#!/usr/bin/env python3
"""Split a node's objective by the 13 threshold regions (the exact identity
d = sum_t |1[a in S_t] - 1[b in S_t]|) and compare, region by region, the DNN
iterate's value of that threshold's cut term with the heuristic (SA) minimum
of the same term alone (a capacitated min-bisection with the node's linear
terms). A threshold whose DNN value lies below its own bisection optimum is
one where the DNN's projection is outside the bisection hull in the
objective's direction; the sum over thresholds of the SA minima is an upper
bound on the fully decoupled bound sum_t MinBis_t.

  python3 threshold_split.py SPEC YFILE [--runs 10 --sweeps 300]
"""
import argparse, itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from decomp_ub import write_sa, run_sa


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y', nargs='?')
    ap.add_argument('--runs', type=int, default=10); ap.add_argument('--sweeps', type=int, default=300)
    a = ap.parse_args()
    nd = gl.Node(a.spec)
    n = nd.n
    regs_full = [('row<%d' % t, np.array([c // 10 < t for c in range(n)])) for t in range(1, 5)] + \
                [('col<%d' % t, np.array([c % 10 < t for c in range(n)])) for t in range(1, 10)]
    if a.Y:
        Y = nd.load_Y(a.Y); Y4, x = nd.blocks(Y); dg = np.einsum('aiai->ai', Y4)
    L = np.array(nd.free_loc); Fc = np.array(nd.free_fac)
    tot_y = tot_sa = 0
    print('region  cap  DNN-iterate  SA-min  (iterate - SA)')
    for name, ind in regs_full:
        Dt = np.abs(ind[:, None].astype(float) - ind[None, :].astype(float))
        D = Dt[np.ix_(L, L)]
        lin = np.zeros((nd.m, nd.m)); c0 = 0.0
        for x_, y_ in nd.fixes:
            lin += 2 * np.outer(Dt[L, x_], nd.B[Fc, y_])
        for (x1, y1), (x2, y2) in itertools.product(nd.fixes, repeat=2):
            c0 += Dt[x1, x2] * nd.B[y1, y2]
        p = f'{gl.ROOT}/experiments/wil50_feasibility/new_cuts/grid/runs/sa_thr.txt'
        write_sa(p, D, nd.F, lin)
        best, _ = run_sa(p, a.runs, a.sweeps, 3)
        best += c0
        vy = np.nan
        if a.Y:
            vy = np.einsum('ab,aibj,ij->', D, Y4, nd.F, optimize=True) + (lin * dg).sum() + c0
            tot_y += vy
        tot_sa += best
        cap = int(ind[L].sum())
        print(f'{name:6s} {cap:4d} {vy:11.2f} {best:9.1f} {vy - best:9.2f}', flush=True)
    print(f'total  {tot_y:11.2f} {tot_sa:9.1f} {tot_y - tot_sa:9.2f}')


if __name__ == '__main__':
    main()
