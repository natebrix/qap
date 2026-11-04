#!/usr/bin/env python3
"""Flow-weighted dispersion / sub-QAP inequalities at a dumped iterate.

For a set T of free facilities, every permutation satisfies
   sum_{i,j in T} F_ij d(loc_i, loc_j) [+ sum_{i in T} lin_i(loc_i)]  >=  Q(T),
Q(T) = the optimum of the |T|-facility sub-QAP on the node's free cells (with
or without the fixed facilities' linear terms). The left side is linear in Y.
Here Q(T) is estimated by SA (an upper bound on Q, so the printed violation is
optimistic: if it is negative the inequality is certainly slack).

T: greedy cliques of high flow (seeded at every facility), and the most
clustered sets of E[d].
  python3 subqap_viol.py SPEC YFILE [--kmin 4 --kmax 12]
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from decomp_ub import write_sa, run_sa

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--kmin', type=int, default=4); ap.add_argument('--kmax', type=int, default=12)
    ap.add_argument('--seeds', type=int, default=8)
    a = ap.parse_args()
    nd = gl.Node(a.spec); Y = nd.load_Y(a.Y); Y4, x = nd.blocks(Y); dg = np.einsum('aiai->ai', Y4)
    Er, Ec = nd.expected_dist(Y4); E = Er + Ec
    D = nd.Dr + nd.Dc; lin = nd.lin_r + nd.lin_c; F = nd.F; m = nd.m
    rng = np.random.default_rng(1)
    seeds = rng.choice(m, a.seeds, replace=False)
    print('k  seed  withlin  LHS(Y)  Q_SA  viol(optimistic)')
    res = []
    for k in range(a.kmin, a.kmax + 1, 2):
        for s0 in seeds:
            T = [int(s0)]
            while len(T) < k:
                rest = [j for j in range(m) if j not in T]
                T.append(max(rest, key=lambda j: F[j, T].sum()))
            for withlin in (0, 1):
                FT = np.zeros_like(F); FT[np.ix_(T, T)] = F[np.ix_(T, T)]
                LT = np.zeros_like(lin)
                if withlin:
                    LT[:, T] = lin[:, T]
                lhs = (FT * E).sum() + (LT * dg).sum()
                p = HERE + '/runs/sa_sub.txt'
                write_sa(p, D, FT, LT)
                q, _ = run_sa(p, 6, 100, 5)
                res.append((k, withlin, lhs, q))
                print(f'{k:2d} {s0:3d} {withlin} {lhs:9.2f} {q:8.1f} {q - lhs:8.2f}', flush=True)


if __name__ == '__main__':
    main()
