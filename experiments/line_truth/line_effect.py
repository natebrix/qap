#!/usr/bin/env python3
"""line_effect.py SHARDDIR: does a child's shortfall say the same about the work under it whichever line it came
from? lines_price.py prices a line the run did not take with the size law fitted on the line it did take; that
holds only if the law does not depend on the line. Where the run took many lines (depth 3 down, rule 9), fit
ln(CPU h under a child) = a + k (its first-solve shortfall) over all children at a depth whose subtrees are closed
and that branched, then look at the residual by the parent's line: rows against columns, and line by line."""
import sys, math, collections
import numpy as np
import tree
ix, nodes, kids, closed, meta = tree.load(sys.argv[1]); G = lambda r, c: r[ix[c]]
print(f"{meta['run']} as of {meta['generatedAt']}: {len(nodes)} nodes")
for d in (3, 4, 5, 6):
    rows = []
    for r in nodes.values():
        if G(r, 'depth') != d or G(r, 'outcome') != 'B' or G(r, 'id') not in closed or G(r, 'dnn') is None: continue
        p = nodes.get(G(r, 'parent'));  s = G(r, 'cutoff') - G(r, 'dnn'); h = (G(r, 'subCpu') or 0) / 3600
        if p is None or s <= 0 or h <= 0.02: continue
        rows.append((s, math.log(h), G(p, 'line'), G(p, 'id')))
    if len(rows) < 50: continue
    s = np.array([x[0] for x in rows]); y = np.array([x[1] for x in rows]); ln = [x[2] for x in rows]
    k, a = np.polyfit(s, y, 1); res = y - (a + k * s)
    print(f"\nchildren at depth {d} (closed, branched, over 0.02 h): {len(rows)} under {len(set(x[3] for x in rows))} parents; ln h = {a:.2f} + {k:.4f} s, rms factor {math.exp(res.std()):.2f}")
    isrow = np.array([l[0] == 'r' for l in ln])
    for nm, m in (("row", isrow), ("column", ~isrow)):
        if m.sum() > 5:
            k2, a2 = np.polyfit(s[m], y[m], 1)
            print(f"  parent took a {nm}: {m.sum()} children, mean residual {res[m].mean():+.3f} (work x{math.exp(res[m].mean()):.2f} at equal shortfall), own slope {k2:.4f}")
    by = collections.defaultdict(list); bp = collections.defaultdict(set)
    for (ss, yy, l, pid), e in zip(rows, res): by[l].append(e); bp[l].add(pid)
    big = sorted((l for l in by if len(by[l]) >= 20), key=lambda l: -len(by[l]))
    print(f"  by line (20 children or more): {len(big)} lines of {len(by)}")
    print("  line   children  parents  mean residual  work at equal shortfall  s.e.")
    m = []
    for l in big:
        e = np.array(by[l]); m.append(e.mean())
        print(f"  {l:5s} {len(e):8d} {len(bp[l]):8d} {e.mean():+13.3f} {math.exp(e.mean()):12.2f}x {e.std()/math.sqrt(len(e)):14.3f}")
    if m: print(f"  spread of the lines' mean residuals: s.d. {np.std(m):.3f} (x{math.exp(np.std(m)):.2f}), from x{math.exp(min(m)):.2f} to x{math.exp(max(m)):.2f}")
