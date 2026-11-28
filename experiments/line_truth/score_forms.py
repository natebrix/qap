#!/usr/bin/env python3
"""score_forms.py SHARDDIR: which function of a line's children's shortfalls tracks the work under the line, in
wil50-261126's own tree. For every branched node whose subtree is closed: the work under its children (sum of
their subCpu) against scores formed from the children's shortfalls after their first DNN solve (what a firing
sees; 0 for a child fathomed there): the plain sum, sums of powers, sums of exp(k s) - 1 (a child with no shortfall adds nothing), the largest. Across nodes,
not across lines at one node: but a score that cannot order nodes by work cannot order lines by it either.
A second table holds the plain sum nearly fixed (pairs of nodes within 5% of each other's sum) and asks how
often each score puts the pair in the right order: the case a line choice turns on."""
import sys, math, itertools, collections
import numpy as np
import tree
ix, nodes, kids, closed, meta = tree.load(sys.argv[1]); G = lambda r, c: r[ix[c]]
def rank(a):
    r = np.empty(len(a)); r[np.argsort(a)] = np.arange(len(a)); return r
forms = [('sum s', lambda s: s.sum()), ('sum s^2', lambda s: (s**2).sum()), ('sum s^4', lambda s: (s**4).sum()),
         ('sum s^8', lambda s: (s**8).sum()), ('sum s^1.5', lambda s: (s**1.5).sum()), ('sum s^3', lambda s: (s**3).sum()),
         ('sum e^.01s - 1', lambda s: np.expm1(.01*s).sum()), ('sum e^.02s - 1', lambda s: np.expm1(.02*s).sum()),
         ('sum e^.03s - 1', lambda s: np.expm1(.03*s).sum()), ('sum e^.04s - 1', lambda s: np.expm1(.04*s).sum()),
         ('max s', lambda s: s.max()), ('open children', lambda s: (s > 0).sum())]
print(f"{meta['run']} as of {meta['generatedAt']}")
for depth in (2, 3, 4, 5):
    X, y = [], []
    for r in nodes.values():
        if G(r, 'depth') != depth or G(r, 'outcome') != 'B' or G(r, 'id') not in closed: continue
        ch = [nodes[c] for c in kids.get(G(r, 'id'), [])]
        s = np.array([max(G(c, 'cutoff') - G(c, 'dnn'), 0.0) if G(c, 'dnn') is not None else 0.0 for c in ch])
        w = sum((G(c, 'subCpu') or 0) for c in ch) / 3600
        if w <= 0 or s.max() <= 0: continue
        X.append([f(s) for _, f in forms]); y.append(w)
    X, y = np.array(X), np.array(y)
    if len(y) < 30: continue
    ly = np.log(y)
    print(f"\ndepth {depth}: {len(y)} closed branched nodes, {y.sum():.0f} CPU h under their children (median {np.median(y):.2f} h, max {y.max():.1f})")
    print("  score            Spearman   r2 of ln work on ln score   right order among pairs within 5% in sum s")
    ss = X[:, 0]; order = np.argsort(ss); pairs = []
    for a in range(len(order)):
        for b in range(a + 1, min(a + 40, len(order))):
            i, j = order[a], order[b]
            if ss[j] > 1.05 * ss[i]: break
            if abs(ly[i] - ly[j]) > math.log(1.5): pairs.append((i, j))
    for c, (nm, _) in enumerate(forms):
        x = X[:, c]; rho = np.corrcoef(rank(x), rank(y))[0, 1]
        lx = np.log(np.maximum(x, 1e-9)); r2 = np.corrcoef(lx, ly)[0, 1] ** 2
        ok = sum((x[i] - x[j]) * (y[i] - y[j]) > 0 for i, j in pairs); 
        print(f"  {nm:15s} {rho:8.3f} {r2:14.3f} {'':14s} {100*ok/max(len(pairs),1):5.1f}% of {len(pairs)}")
