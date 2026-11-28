#!/usr/bin/env python3
"""lines_price.py SHARDDIR [d1|d2]: every line a firing solved in full at wil50's depth-1 nodes (../wil50_sb/runs)
or at four heavy depth-2 nodes (d2/runs), priced by the size law of wil50-261126's own tree.

The law is fitted where both sides are known: the children of the line the run took (r45 at depth 1, r31 at
depth 2) have a shortfall in the firing (cutoff less the 150-iteration bound) and a subtree in the run; ln(CPU h
under the child so far) = a + k shortfall over those children, pooled over the fired nodes of that depth. A
line's price is the sum of exp(a + k s) over its open children. Then: the line each score would take (sum of
(s/smax)^P for P = 1, 2, 4, 8; the lowest price itself), and its price next to the run's line. Open subtrees
are still growing and the heaviest are the most open, so k is if anything too small."""
import sys, os, glob, math, collections
import numpy as np
import tree
H = os.path.dirname(os.path.abspath(__file__))
ix, nodes, kids, closed, meta = tree.load(sys.argv[1]); G = lambda r, c: r[ix[c]]
which = sys.argv[2] if len(sys.argv) > 2 else 'd1'
pat, depth = (('../wil50_sb/runs/21-*.progress.txt', 1) if which == 'd1' else ('d2/runs/*.progress.txt', 2))
def path(r):
    f = []
    while r is not None and G(r, 'depth') > 0: f.append(G(r, 'fix')); r = nodes.get(G(r, 'parent'))
    return ','.join(reversed(f))
bypath = {path(r): r for r in nodes.values() if G(r, 'depth') == depth}
fired = {}
for f in sorted(glob.glob(os.path.join(H, pat))):
    node = os.path.basename(f).split('.')[0].replace('-', ':').replace('_', ',')
    cells = {}
    for ln in open(f):
        t = ln.split()
        if t[0] != 'sb-child': continue
        a, b = t[4].split(':'); fath = t[12] == '1'
        cells[(int(a), int(b))] = 0.0 if fath else max(float(t[8]) - float(t[10]), 0.0)
    if node in bypath: fired[node] = cells
print(f"{meta['run']} as of {meta['generatedAt']}; firings at depth {depth}: {', '.join(fired)}")
# the law, from the run's line
X, Y = [], []; runline = {}
for node, cells in fired.items():
    p = bypath[node]; live = {}
    for c in kids.get(G(p, 'id'), []):
        r = nodes[c]; a, b = G(r, 'fix').split(':'); live[(int(a), int(b))] = r
    runline[node] = live
    for k, r in live.items():
        if k in cells and cells[k] > 0 and (G(r, 'subCpu') or 0) / 3600 > 0.05: X.append(cells[k]); Y.append(math.log(G(r, 'subCpu') / 3600))
X, Y = np.array(X), np.array(Y)
if len(X) < 10: sys.exit("the run's line is not solved in these firings yet: no law")
k, a = np.polyfit(X, Y, 1); res = Y - (a + k * X)
print(f"law from {len(X)} children of the run's line: ln h = {a:.2f} + {k:.4f} s (x{math.exp(100*k):.1f} per 100 units), rms factor {math.exp(res.std()):.2f}, r2 {1 - res.var()/Y.var():.2f}")
price = lambda s: float(np.exp(a + k * s[s > 0]).sum())
tot = collections.defaultdict(float)
for node, cells in fired.items():
    p = bypath[node]; live = runline[node]
    nfac = collections.Counter(x for x, y in cells); nloc = collections.Counter(y for x, y in cells)
    full = max(max(nfac.values()), max(nloc.values()))
    lines = {('r%d' % x): np.array([cells[c] for c in cells if c[0] == x]) for x in nfac if nfac[x] == full}
    lines.update({('c%d' % y): np.array([cells[c] for c in cells if c[1] == y]) for y in nloc if nloc[y] == full})
    took = G(p, 'line'); act = sum((G(r, 'subCpu') or 0) for r in live.values()) / 3600
    nopen = sum(G(r, 'id') not in closed for r in live.values())
    smax = max(v.max() for v in lines.values()); pr = {l: price(v) for l, v in lines.items()}
    print(f"\n{node}: the run took {took}; {act:.0f} CPU h under its children so far ({nopen} of {len(live)} subtrees open); {len(lines)} lines solved in full, largest shortfall {smax:.0f}")
    if took not in lines: print(f"  {took} is not solved in full in this firing yet"); continue
    print("  line   open  largest   sum s   price h   price / run's line")
    for l in sorted(lines, key=lambda l: pr[l]):
        v = lines[l]; print(f"  {l:5s} {int((v > 0).sum()):5d} {v.max():8.1f} {v.sum():7.0f} {pr[l]:9.1f} {pr[l]/pr[took]:8.2f}{'   <- the run' if l == took else ''}")
    picks = []
    for P in (1, 2, 4, 8):
        sc = {l: ((v / smax) ** P).sum() for l, v in lines.items()}; b = min(sc, key=sc.get)
        picks.append(f"P = {P}: {b} ({pr[b]/pr[took]:.2f})"); tot[P] += pr[b]
    b = min(pr, key=pr.get); tot['price'] += pr[b]; tot['run'] += pr[took]; tot['act'] += act
    print("  a score's line (its price / the run's line's): " + "; ".join(picks) + f"; lowest price: {b} ({pr[b]/pr[took]:.2f})")
if tot['run']:
    print(f"\nover these nodes: the run's lines price at {tot['run']:.0f} h (work under them so far {tot['act']:.0f} h); the lines of P = 1: {tot[1]/tot['run']:.2f}, "
          f"P = 2: {tot[2]/tot['run']:.2f}, P = 4: {tot[4]/tot['run']:.2f}, P = 8: {tot[8]/tot['run']:.2f}, lowest price: {tot['price']/tot['run']:.2f} of that")
