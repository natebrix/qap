#!/usr/bin/env python3
"""Where the DNN 'cheats' on the facility side: decompose (value of a good
permutation) - (DNN primal) into facility pairs: g_jl = B_jl (D*_jl - Delta_jl),
plus the linear part. python3 gapdecomp.py SPEC Y.bin PERM.txt"""
import sys
import numpy as np
from ynode import Node
nd = Node(sys.argv[1]).load(sys.argv[2]); m = nd.m
full = list(map(int, open(sys.argv[3]).read().split()))
pos = {f: a for a, f in enumerate(full)}          # facility (orig) -> location (orig)
locl = {a: i for i, a in enumerate(nd.loc)}
pi_loc = np.array([locl[pos[f]] for f in nd.fac])  # local facility j -> local location
Dst = nd.A[np.ix_(pi_loc, pi_loc)]
D = nd.delta()
G = nd.B * (Dst - D)
xs = np.zeros((m, m)); xs[pi_loc, np.arange(m)] = 1
lin = float((nd.C * (xs - nd.x)).sum())
val = float((nd.B * Dst).sum() + (nd.C * xs).sum() + nd.shift)
print('permutation value %.0f, DNN primal %.1f, gap %.1f = quadratic %.1f + linear %.1f' % (val, nd.objective(), val - nd.objective(), G.sum(), lin))
iu = np.triu_indices(m, 1)
g = 2 * G[iu]; order = np.argsort(-g)
print('pairs: positive part %.1f, negative part %.1f; top 20 pairs carry %.1f; pairs with B>=8 carry %.1f of %.1f' % (
    g[g > 0].sum(), g[g < 0].sum(), g[order[:20]].sum(), g[nd.B[iu] >= 8].sum(), g.sum()))
per = G.sum(1) * 1.0
deg = nd.B.sum(1)
print('per facility (sum over its pairs): min %.1f max %.1f; corr with flow degree %.2f; top 5 facilities carry %.1f' % (
    per.min(), per.max(), np.corrcoef(per, deg)[0, 1], np.sort(per)[-5:].sum()))
# distance profile: for heavy pairs, solution distance vs DNN expected distance
for lo, hi in ((0, 2), (3, 5), (6, 7), (8, 9)):
    sel = (nd.B[iu] >= lo) & (nd.B[iu] <= hi)
    print('  flows %d-%d: %4d pairs, solution mean distance %.3f, DNN %.3f' % (lo, hi, sel.sum(), Dst[iu][sel].mean(), D[iu][sel].mean()))
