#!/usr/bin/env python3
"""Heavier SA on a coupled threshold set K at a node, and the K-part of the
node's best full SA solution.  python3 check_coupled.py SPEC "row<2,row<3,..." RUNS SWEEPS"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from decomp_ub import write_sa, run_sa
from coupled_split import thr_matrix
HERE = os.path.dirname(os.path.abspath(__file__))
spec, K, runs, sweeps = sys.argv[1], sys.argv[2].split(','), int(sys.argv[3]), int(sys.argv[4])
nd = gl.Node(spec); L = np.array(nd.free_loc); Fc = np.array(nd.free_fac)
def parts(Dfull):
    D = Dfull[np.ix_(L, L)]; lin = np.zeros((nd.m, nd.m)); c0 = 0.0
    for x_, y_ in nd.fixes: lin += 2 * np.outer(Dfull[L, x_], nd.B[Fc, y_])
    for (x1, y1), (x2, y2) in itertools.product(nd.fixes, repeat=2): c0 += Dfull[x1, x2] * nd.B[y1, y2]
    return D, lin, c0
DK = sum(thr_matrix(nm) for nm in K)
D, lin, c0 = parts(DK)
write_sa(HERE + '/runs/sa_chk.txt', D, nd.F, lin)
for seed in (21, 22):
    best, perm = run_sa(HERE + '/runs/sa_chk.txt', runs, sweeps, seed)
    print(f'K-part SA optimum (seed {seed}): {best + c0:.0f}', flush=True)
# best full solution and its K part
Dr, Dc = gl.grid_dist(); Df, linf, cf = parts(Dr + Dc)
write_sa(HERE + '/runs/sa_chkfull.txt', Df, nd.F, linf)
bf, pf = run_sa(HERE + '/runs/sa_chkfull.txt', runs, sweeps, 23)
pf = np.array(pf)
kval = (nd.F * D[np.ix_(pf, pf)]).sum() + lin[pf, np.arange(nd.m)].sum() + c0
print(f'full SA {bf + cf:.0f}; its K part {kval:.0f}, rest {bf + cf - kval:.0f}')
