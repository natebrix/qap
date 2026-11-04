#!/usr/bin/env python3
"""Pick a count split for the disjunction test: over the 13 threshold regions,
search signed facility sets s (support 4..15) for the most violated box
inequality (Var(c) < f(1 - f)), i.e. the count whose degree-2 split cut the
iterate violates most; print region, Q+, Q-, beta = floor(E c), E c, Var c.
  python3 pick_split.py SPEC Y.bin"""
import sys
import numpy as np
from ynode import Node
from box_sep import viol
nd = Node(sys.argv[1]).load(sys.argv[2]); m = nd.m; x = nd.x; T = nd.T
rng = np.random.default_rng(5)
best = (-1, None)
regs = [('row<%d' % t, nd.rc[:, 0] < t) for t in range(1, 5)] + [('col<%d' % t, nd.rc[:, 1] < t) for t in range(1, 10)]
if len(sys.argv) > 3: regs = [r for r in regs if r[0] in sys.argv[3].split(',')]   # optional region filter
for name, mask in regs:
    idx = np.where(mask)[0]
    p = x[idx].sum(0); W = T[np.ix_(idx, range(m), idx, range(m))].sum(axis=(0, 2))
    Wd = W.copy(); np.fill_diagonal(Wd, p); Sig = Wd - np.outer(p, p)
    for st in range(80):
        k = rng.integers(4, 10); s = np.zeros(m); I = rng.choice(m, k, replace=False); s[I] = rng.choice([-1, 1], k)
        cur = viol(s, p, Sig); imp = True
        while imp:
            imp = False
            for i in rng.permutation(m):
                for v in (-1.0, 0.0, 1.0):
                    if v == s[i]: continue
                    o = s[i]; s[i] = v; nz = np.count_nonzero(s)
                    if nz < 4 or nz > 15: s[i] = o; continue
                    nv = viol(s, p, Sig)
                    if nv > cur + 1e-12: cur, imp = nv, True
                    else: s[i] = o
        if cur > best[0]:
            mu = s @ p + (s < 0).sum(); var = s @ Sig @ s
            best = (cur, name, np.where(s > 0)[0].tolist(), np.where(s < 0)[0].tolist(), int(np.floor(mu)), mu, var)
v, name, qp, qm, beta, mu, var = best
print('best split: region %s  Q+ %s  Q- %s  beta %d  E c %.4f  Var c %.4f  box violation %.4f (triangle units %.4f)' % (name, qp, qm, beta, mu, var, v, v / 2))
print('ARGS --region %s --plus %s%s --beta %d' % (name, ','.join(map(str, qp)), (' --minus ' + ','.join(map(str, qm))) if qm else '', beta))
