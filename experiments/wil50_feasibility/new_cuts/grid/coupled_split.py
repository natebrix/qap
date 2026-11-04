#!/usr/bin/env python3
"""Coupled threshold subsets. For a set K of threshold regions, the part of
the objective carried by K is sum_{t in K} cut_t (with the node's linear
terms); minimizing it alone over permutations couples the thresholds in K
through the joint cell capacities (for K = {row<t, col<s}: a 4-way quadrant
partition with fixed quadrant capacities). Compare, for each K:
  DNN iterate's value of the K-part | SA optimum of the K-part alone (coupled)
  | sum of the SA optima of each threshold alone (decoupled).
(iterate - coupled) < 0: the iterate is outside the coupled hull of K in the
objective direction; (coupled - decoupled): the value of coupling K.

  python3 coupled_split.py SPEC YFILE [--pairs] [--sets "row<3,col<5;row<2,col<4,col<5"]
"""
import argparse, itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from decomp_ub import write_sa, run_sa

HERE = os.path.dirname(os.path.abspath(__file__))


def thr_matrix(name, n=50):
    t = int(name[4:]); ind = np.array([(c // 10 if name.startswith('row') else c % 10) < t for c in range(n)])
    return np.abs(ind[:, None].astype(float) - ind[None, :].astype(float))


def part_value(nd, Y4, dg, Dfull):
    L = np.array(nd.free_loc); Fc = np.array(nd.free_fac)
    D = Dfull[np.ix_(L, L)]
    lin = np.zeros((nd.m, nd.m)); c0 = 0.0
    for x_, y_ in nd.fixes:
        lin += 2 * np.outer(Dfull[L, x_], nd.B[Fc, y_])
    for (x1, y1), (x2, y2) in itertools.product(nd.fixes, repeat=2):
        c0 += Dfull[x1, x2] * nd.B[y1, y2]
    vy = np.einsum('ab,aibj,ij->', D, Y4, nd.F, optimize=True) + (lin * dg).sum() + c0 if Y4 is not None else np.nan
    p = HERE + '/runs/sa_cpl.txt'
    write_sa(p, D, nd.F, lin)
    best, _ = run_sa(p, 8, 200, 11)
    return vy, best + c0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--pairs', action='store_true'); ap.add_argument('--sets', default='')
    a = ap.parse_args()
    nd = gl.Node(a.spec); Y = nd.load_Y(a.Y); Y4, _ = nd.blocks(Y); dg = np.einsum('aiai->ai', Y4)
    names = [f'row<{t}' for t in range(1, 5)] + [f'col<{t}' for t in range(1, 10)]
    single = {}
    for nm in names:
        single[nm] = part_value(nd, Y4, dg, thr_matrix(nm))
    sets = []
    if a.pairs:
        sets += [[r, c] for r in names[:4] for c in names[4:]]
    sets += [s.split(',') for s in a.sets.split(';') if s]
    print('K | iterate  coupled  decoupled | iterate-coupled  coupled-decoupled')
    for K in sets:
        Dfull = sum(thr_matrix(nm) for nm in K)
        vy, cpl = part_value(nd, Y4, dg, Dfull)
        dec = sum(single[nm][1] for nm in K)
        print(f'{"+".join(K):30s} | {vy:9.1f} {cpl:8.0f} {dec:8.0f} | {vy - cpl:8.1f} {cpl - dec:7.0f}', flush=True)


if __name__ == '__main__':
    main()
