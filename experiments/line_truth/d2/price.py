#!/usr/bin/env python3
"""price.py SHARDDIR: what the lines tried at wil50's heavy depth-2 nodes (d2/runs) would cost, next to row 31,
the line wil50-261126 took there.

For each firing: every line with all its children solved (from the progress file: a row is a facility with every
free location solved, a column a location with every free facility), its children's shortfalls (cutoff less the
150-iteration bound; 0 when fathomed), the scores strong branching can form from them, and a price in CPU hours,
sum over the open children of exp(a + k shortfall). a and k are fitted here on the run's own depth-3 nodes: the
work under each (subCpu) against its shortfall after its first DNN solve, the quantity a firing sees. Open
subtrees are still growing, so the fit is given three ways: closed subtrees only, all, and all under the nodes
fired here. Row 31's children are in the run's tree, so its price has a direct check: their subtree CPU so far."""
import sys, os, glob, math, collections
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import tree
D = os.path.dirname(os.path.abspath(__file__))
ix, nodes, kids, closed, meta = tree.load(sys.argv[1])
G = lambda r, c: r[ix[c]]
print(f"{meta['run']} as of {meta['generatedAt']}: {len(nodes)} nodes")

def path(r):
    f = []
    while r is not None and G(r, 'depth') > 0:
        f.append(G(r, 'fix')); r = nodes.get(G(r, 'parent'))
    return ','.join(reversed(f))

# the size law at depth 3, on the shortfall after the first DNN solve
d3 = [r for r in nodes.values() if G(r, 'depth') == 3 and G(r, 'dnn') is not None and G(r, 'cutoff') is not None]
def fit(rs, label):
    x = np.array([G(r, 'cutoff') - G(r, 'dnn') for r in rs]); y = np.array([(G(r, 'subCpu') or 0) / 3600 for r in rs])
    m = (x > 0) & (y > 0.02)
    if m.sum() < 10: print(f"  {label}: too few"); return None
    k, a = np.polyfit(x[m], np.log(y[m]), 1); res = np.log(y[m]) - (a + k * x[m])
    print(f"  {label}: {m.sum()} nodes, ln h = {a:.2f} + {k:.4f} s (x{math.exp(100*k):.1f} per 100), rms factor {math.exp(res.std()):.2f}, "
          f"r2 {1 - res.var()/np.log(y[m]).var():.2f}")
    return a, k
print("size law at depth 3 (work under a node against its shortfall after its first DNN solve):")
laws = {}
br = [r for r in d3 if G(r, 'outcome') == 'B']
laws['closed'] = fit([r for r in br if G(r, 'id') in closed], 'branched, subtree closed')
laws['all'] = fit(br, 'branched, all (open ones still growing)')

fired = {}
for f in sorted(glob.glob(os.path.join(D, 'runs', '*.progress.txt'))):
    name = os.path.basename(f).split('.')[0]; node = name.replace('-', ':').replace('_', ',')
    cells = {}
    for ln in open(f):
        t = ln.split()
        if t[0] != 'sb-child': continue
        a, b = t[4].split(':'); cut = float(t[8]); lb = float(t[10]); fath = t[12] == '1'
        cells[(int(a), int(b))] = (0.0 if fath else max(cut - lb, 0.0), fath)
    fired[node] = cells
under = [r for r in br if path(nodes[G(r, 'parent')]) in fired]
laws['fired'] = fit(under, 'branched, under the four nodes fired here')

def price(sh, law):
    a, k = law; return sum(math.exp(a + k * s) for s, f in sh if not f)
for node, cells in fired.items():
    par = [r for r in nodes.values() if G(r, 'depth') == 2 and path(r) == node]
    live = {}
    if par:
        p = par[0]
        for c in kids.get(G(p, 'id'), []):
            r = nodes[c]; a, b = G(r, 'fix').split(':'); live[(int(a), int(b))] = r
        print(f"\n{node}: in the run, line {G(p, 'line')} (rule {G(p, 'rule')}), bound {G(p, 'lb'):.1f}, shortfall {G(p, 'cutoff') - G(p, 'lb'):.1f}, "
              f"{G(p, 'kids')} children, {(G(p, 'subCpu') or 0)/3600:.1f} CPU h under it so far, {'closed' if G(p, 'id') in closed else 'open'}")
    else: print(f"\n{node}: not in the tree")
    nfac = collections.Counter(a for a, b in cells); nloc = collections.Counter(b for a, b in cells)
    full = max(max(nfac.values()), max(nloc.values()))
    lines = [('r%d' % a, [cells[c] for c in cells if c[0] == a]) for a in nfac if nfac[a] == full]
    lines += [('c%d' % b, [cells[c] for c in cells if c[1] == b]) for b in nloc if nloc[b] == full and full > 8]
    print(f"  {len(cells)} children solved; complete lines ({full} children each):")
    print("  line  open  largest  sum s   sum s^8 (rel)   price h: closed-fit  all-fit  fired-fit")
    out = []
    for nm, sh in lines:
        s = np.array([x for x, f in sh]); smax = max(s.max(), 1e-9)
        pr = [price(sh, laws[w]) if laws[w] else float('nan') for w in ('closed', 'all', 'fired')]
        out.append((pr[2] if pr[2] == pr[2] else pr[1], nm, int(sum(not f for x, f in sh)), s.max(), s.sum(), float(((s / 400.0) ** 8).sum()), pr))
    for _, nm, op, mx, sm, p8, pr in sorted(out):
        print(f"  {nm:5s} {op:4d} {mx:8.1f} {sm:7.0f} {p8:12.4f}        {pr[0]:10.1f} {pr[1]:9.1f} {pr[2]:9.1f}")
    # row 31 in the run: its children's first DNN shortfalls priced the same way, and the work under them so far
    l31 = [(k, r) for k, r in live.items()]
    if l31:
        sh = []
        for k, r in l31:
            dn = G(r, 'dnn'); fa = G(r, 'outcome') in tree.FATH and (G(r, 'kids') or 0) == 0
            s = max(G(r, 'cutoff') - dn, 0.0) if dn is not None else 0.0
            sh.append((s, dn is None or s <= 0))
        act = sum((G(r, 'subCpu') or 0) for k, r in l31) / 3600
        nopen = sum(G(r, 'id') not in closed for k, r in l31)
        pr = [price(sh, laws[w]) if laws[w] else float('nan') for w in ('closed', 'all', 'fired')]
        s = np.array([x for x, f in sh])
        print(f"  run's line from the tree ({len(l31)} children, first-solve shortfalls): open {int(sum(not f for x, f in sh))}, largest {s.max():.1f}, "
              f"sum {s.sum():.0f}; price {pr[0]:.1f} / {pr[1]:.1f} / {pr[2]:.1f} h; work under them so far {act:.1f} h, {nopen} subtrees open")
        both = [(cells[k][0], max(G(r, 'cutoff') - G(r, 'dnn'), 0.0)) for k, r in l31 if k in cells and G(r, 'dnn') is not None]
        if both:
            b = np.array(both); d = b[:, 0] - b[:, 1]
            print(f"  same children here and in the run ({len(both)}): shortfall here less the run's first solve: mean {d.mean():+.1f}, rms {math.sqrt((d**2).mean()):.1f}")
