#!/usr/bin/env python3
"""An external objective term (dnncuts -Ladd L.bin:CONST) from the K most
violated box inequalities on the INTERIOR lines of the grid (rows 1..3,
columns 1..8; the edge lines are the complements of the edge thresholds, whose
boxes new_cuts/box_pool killed), separated at a dumped iterate, at a fixed
multiplier. Same convention as new_cuts/literature/wil50/make_box_ladd.py so
the two screens compare: per cut <G, Y> <= h with G = 2b on (a,i),(a,i) for a
in S, i in Q and -1 on (a,i),(c,j) for a != c in S, i != j in Q (gangster
entries left at 0), h = b (b + 1); L = mu sum_k G_k, CONST = mu sum_k h_k.
Every cut is checked on random permutations of the node and on permutations
realising each count before it is written (a violated one aborts).

  python3 make_line_ladd.py SPEC Y.bin --mu 0.3 --top 200 --out L.bin
"""
import argparse

import numpy as np

from nlc_analyze import Node, ROWS, COLS
from line_boxes import moments, separate_boxes

ap = argparse.ArgumentParser()
ap.add_argument('spec'); ap.add_argument('ybin')
ap.add_argument('--mu', type=float, default=0.3)
ap.add_argument('--top', type=int, default=200)
ap.add_argument('--out', required=True)
ap.add_argument('--nperm', type=int, default=300)
ap.add_argument('--starts', type=int, default=40)
ap.add_argument('--all-lines', action='store_true', help='edge lines too')
a = ap.parse_args()

nd = Node(a.spec).load(a.ybin)
m = nd.m; N = m * m + 1
lines = {}
for r in range(ROWS):
    if a.all_lines or 0 < r < ROWS - 1:
        lines[f'row={r}'] = [c for c in range(m) if nd.r[c] == r]
for c_ in range(COLS):
    if a.all_lines or 0 < c_ < COLS - 1:
        lines[f'col={c_}'] = [c for c in range(m) if nd.c[c] == c_]
cands = []
for name, cells in lines.items():
    if not 0 < len(cells) < m:
        continue
    p, W = moments(nd, cells)
    found = separate_boxes(p, W, m, a.starts, 24, 7, seed=hash(name) & 0xffff)
    for (b, Q), v in found.items():
        cands.append((v, name, list(Q), b, cells))
cands.sort(key=lambda t: -t[0])
cuts = cands[:a.top]
print(f'{len(cands)} violated interior-line boxes; taking {len(cuts)}: violations {cuts[-1][0]:.4f}..{cuts[0][0]:.4f}')
from collections import Counter
print('  per line:', dict(Counter(c[1] for c in cuts)))
print('  |Q| median', np.median([len(c[2]) for c in cuts]), ' b values', dict(Counter(c[3] for c in cuts)))

rng = np.random.default_rng(7)
perms = [rng.permutation(m) for _ in range(a.nperm)]   # perm[a] = facility at local cell a
L = np.zeros((N, N)); const = 0.0; tight = 0
for (v, name, Q, b, S) in cuts:
    idx = np.array([1 + s * m + j for s in S for j in Q])
    G = np.zeros((N, N))
    G[np.ix_(idx, idx)] = -1.0
    for s in S:                       # same cell, different facilities: gangster zero
        for j in Q:
            for k in Q:
                if j != k:
                    G[1 + s * m + j, 1 + s * m + k] = 0.0
    for j in Q:                       # different cells, same facility: gangster zero
        for s in S:
            for c in S:
                if s != c:
                    G[1 + s * m + j, 1 + c * m + j] = 0.0
    for e in idx:
        G[e, e] = 2.0 * b
    h = float(b * (b + 1))
    Sl = list(S); Ql = list(Q); notS = [x for x in range(m) if x not in S]; notQ = [j for j in range(m) if j not in Q]
    tests = list(perms)
    for t in range(0, min(len(Sl), len(Ql)) + 1):
        if len(Ql) - t > len(notS) or len(Sl) - t > len(notQ):
            continue
        pi = -np.ones(m, dtype=int)
        locs_in = rng.permutation(Sl); locs_out = rng.permutation(notS); facs = rng.permutation(Ql)
        for i in range(t):
            pi[locs_in[i]] = facs[i]
        for i in range(len(Ql) - t):
            pi[locs_out[i]] = facs[t + i]
        rest_f = list(rng.permutation(notQ)); free_l = [x for x in range(m) if pi[x] < 0]
        for x, f in zip(free_l, rest_f):
            pi[x] = f
        tests.append(pi)
    worst = -np.inf
    for pi in tests:
        idxp = 1 + np.arange(m) * m + pi
        worst = max(worst, G[np.ix_(idxp, idxp)].sum() - h)
    assert worst <= 1e-9, ('INVALID cut', name, Q, b, worst)
    tight += abs(worst) < 1e-9
    # the violation at the iterate in this convention (= 2 x line_boxes.py's)
    L += a.mu * G; const += a.mu * h
L.tofile(a.out)
print(f'{len(cuts)} cuts ({tight} tight on a test permutation; all checked on {a.nperm} random permutations), '
      f'mu {a.mu}, nonzeros {np.count_nonzero(L)}, |L| {np.abs(L).sum():.3f}, const {const:.6f} -> use: -Ladd {a.out}:{const:.9g}')
