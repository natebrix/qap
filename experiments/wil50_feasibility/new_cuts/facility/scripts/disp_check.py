#!/usr/bin/env python3
"""Dispersion inequalities on facility subsets: for k facilities, the sum of
their pairwise distances is at least sigma_k = min over k distinct free
locations of the sum of pairwise distances (exact by enumeration of compact
k-sets). Implied by the parity triples for k = 4 (sigma_4 = 8), not for k >= 5
(sigma_5 = 16 > 14). Finds the k-subsets with the smallest sum of expected
distances Delta (greedy + swaps) and reports the slack.
  python3 disp_check.py SPEC Y.bin"""
import sys, itertools
import numpy as np
from ynode import Node
nd = Node(sys.argv[1]).load(sys.argv[2]); m = nd.m; A = nd.A; D = nd.delta()
rng = np.random.default_rng(0)
def sigma(k):
    best = np.inf
    for a0 in range(m):
        near = [b for b in range(m) if b != a0 and A[a0, b] <= 3]
        for rest in itertools.combinations(near, k - 1):
            S = (a0,) + rest
            best = min(best, A[np.ix_(S, S)].sum() / 2)
    return best
for k in (4, 5, 6, 7):
    sg = sigma(k)
    best = (np.inf, None)
    for start in range(300):
        S = list(rng.choice(m, k, replace=False))
        cur = D[np.ix_(S, S)].sum() / 2; improved = True
        while improved:
            improved = False
            for i in range(k):
                for j in range(m):
                    if j in S: continue
                    T = S.copy(); T[i] = j; v = D[np.ix_(T, T)].sum() / 2
                    if v < cur - 1e-12: S, cur, improved = T, v, True
        if cur < best[0]: best = (cur, sorted(S))
    print('k=%d: sigma_k (free locations) = %g; smallest DNN sum of expected distances %.4f (slack %.4f) on %s' % (k, sg, best[0], best[0] - sg, best[1]))
