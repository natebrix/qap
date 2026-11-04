#!/usr/bin/env python3
"""Star inequalities: for a hub facility j and weights w_l = B_jl, sum_l w_l
d(pi_j, pi_l) >= min_a LAP_a with LAP_a = min over assignments of the other
facilities to the other free locations of sum_l w_l A_(a, loc(l)). Implied by
the DNN (the conditional block Y_(aj),(.)/x_aj is doubly stochastic), so the
slack below must be >= 0 up to the iterate's infeasibility. Also reports the
slack against the conditional bound sum_a x_aj LAP_a (tight when the DNN puts
each hub position's conditional on an optimal assignment)."""
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from ynode import Node
nd = Node(sys.argv[1]).load(sys.argv[2]); m = nd.m; A = nd.A; B = nd.B; x = nd.x
D = nd.delta()
sl_min, sl_cond = [], []
for j in range(m):
    w = B[j].copy(); w[j] = 0
    others = [l for l in range(m) if l != j]
    lap = np.zeros(m)
    for a in range(m):
        locs = [b for b in range(m) if b != a]
        C = np.outer(A[a, locs], w[others])          # cost of facility l at location b
        r, c = linear_sum_assignment(C); lap[a] = C[r, c].sum()
    val = (w * D[j]).sum()
    sl_min.append(val - lap.min()); sl_cond.append(val - (x[:, j] * lap).sum())
sl_min, sl_cond = np.array(sl_min), np.array(sl_cond)
print('star (hub j, weights B_j.): DNN value - min_a LAP_a: min %.3f median %.1f; DNN value - sum_a x_aj LAP_a: min %.4f median %.3f max %.3f' % (
    sl_min.min(), np.median(sl_min), sl_cond.min(), np.median(sl_cond), sl_cond.max()))
