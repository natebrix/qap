#!/usr/bin/env python3
"""where.py SHARDDIR: where in the tree a strong-branching firing would pay.
(1) The size law by depth: among branched nodes, ln(subtree CPU) against the node's own shortfall when it branched
    (cutoff - its final bound), closed subtrees and all.
(2) A gate on that shortfall: at each depth, the branched nodes at least S short: how many, what a firing at each
    would cost (FIRE CPU h, a flat figure), and how much of the run's CPU so far lies under them."""
import sys, math, numpy as np
from tree import load
FIRE = 4.0   # CPU h a firing: 8 lines x 48 children x 65 s x 0.57 (sb.early), wil50 at depth 2
ix, nodes, kids, closed, meta = load(sys.argv[1])
root = [r for r in nodes.values() if r[ix['depth']] == 0][0]
TOT = root[ix['subCpu']] / 3600
print(f"{meta['run']} at {meta['generatedAt']}: {len(nodes)} nodes, {TOT:.0f} CPU h under the root, cutoff {root[ix['cutoff']]}")
B = {}
for r in nodes.values():
    if r[ix['outcome']] == 'B' and r[ix['depth']] >= 1 and r[ix['lb']] is not None:
        B.setdefault(r[ix['depth']], []).append((r[ix['cutoff']] - r[ix['lb']], r[ix['subCpu']] / 3600, r[ix['id']] in closed, r[ix['size']]))
print("\n(1) ln(subtree CPU h) = a + k s over branched nodes (s: the node's shortfall when it branched)")
print("depth  branched  closed | closed only: k, x per 100, rms factor, R2 | all: k, x per 100, R2 | median s, median CPU h")
for d in sorted(B):
    a = np.array([(s, w, c) for s, w, c, n in B[d]]); 
    if len(a) < 12: continue
    out = []
    for sel in (a[:, 2] == 1, np.ones(len(a), bool)):
        x, y = a[sel, 0], np.log(np.maximum(a[sel, 1], 1e-4))
        if len(x) < 8: out.append("      -"); continue
        k, c0 = np.polyfit(x, y, 1); res = y - (c0 + k * x)
        out.append(f"k {k:.4f} x{math.exp(100*k):6.1f} rms x{math.exp(res.std()):.2f} R2 {1 - res.var()/y.var():.2f}")
    print(f"{d:5d} {len(a):8d} {int(a[:,2].sum()):7d} | {out[0]} | {out[1]} | {np.median(a[:,0]):5.0f} {np.median(a[:,1]):7.2f}")
print(f"\n(2) a firing at every branched node at least S short, {FIRE:.0f} CPU h each: nodes, firing CPU h, CPU under them (h, share of the run so far), work per firing hour")
for d in (2, 3, 4):
    a = sorted(B.get(d, []), key=lambda t: -t[0])
    print(f"depth {d}: {len(a)} branched nodes, {sum(t[1] for t in a):.0f} CPU h under them ({100*sum(t[1] for t in a)/TOT:.0f}%)")
    for S in (400, 350, 300, 250, 200, 150, 100):
        g = [t for t in a if t[0] >= S]
        if not g: continue
        w = sum(t[1] for t in g)
        small = sum(1 for t in g if t[1] < 2 * FIRE)
        print(f"   S >= {S}: {len(g):5d} nodes, firing {len(g)*FIRE:6.0f} h, under them {w:6.0f} h ({100*w/TOT:4.1f}%), {w/(len(g)*FIRE):6.1f} h of work per firing hour; {small} of them under {2*FIRE:.0f} h so far")
