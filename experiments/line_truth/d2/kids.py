#!/usr/bin/env python3
"""kids.py SHARDDIR: the children solved to the end by kids.sh at wil50's node 21:1,45:0, against the law that
priced their lines. For each: its line, its shortfall in the firing, what the law says (lines_price.py d2's fit
on r31's children in the run), and the CPU h and nodes of its subtree here (cold start, one process). r31's
children are the control: the run solved them too (from a warm start), so the cold start's own effect shows.
Then: at equal shortfall, do the children of c9 and c49 cost what r31's do?"""
import sys, os, math, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..'))
import tree
ix, nodes, kids, closed, meta = tree.load(sys.argv[1]); G = lambda r, c: r[ix[c]]
A, K = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (-3.27, 0.0262)
def path(r):
    f = []
    while r is not None and G(r, 'depth') > 0: f.append(G(r, 'fix')); r = nodes.get(G(r, 'parent'))
    return ','.join(reversed(f))
d3 = {path(r): r for r in nodes.values() if G(r, 'depth') == 3 and path(r).startswith('21:1,45:0,')}
done = {}
for ln in open(os.path.join(D, 'kids', 'done.txt')):
    t = ln.split()
    if t[0] == 'all': continue
    out = open(os.path.join(D, 'kids', t[0] + '.out.txt')).read()
    m = re.search(r'^\s*Nodes: (\d+)', out, re.M); tt = re.search(r'^\s*Total Time: ([\d.]+)', out, re.M)
    done[t[0]] = (int(t[2]), int(m.group(1)) if m else None, float(tt.group(1)) / 3600 if tt else None)
print(f"law: ln h = {A} + {K} s; the run's tree as of {meta['generatedAt']}")
print("child              line  shortfall  law h   here h   nodes   here/law   in the run h")
rows = []
for ln in open(os.path.join(D, 'kids', 'list.txt')):
    node, line, _, s, *_ = ln.split(); s = float(s); tag = node.replace(':', '-').replace(',', '_')
    if tag not in done: continue
    ex, nn, h = done[tag]; law = math.exp(A + K * s); r = d3.get(node)
    run = f"{(G(r, 'subCpu') or 0)/3600:.2f}{'' if G(r, 'id') in closed else ' (open)'}" if r is not None else ''
    if h is None: print(f"{node:18s} {line:4s} {s:9.1f} {law:6.2f}   no result (exit {ex})"); continue
    cap = ' (capped: at least)' if ex == 3 else ''
    print(f"{node:18s} {line:4s} {s:9.1f} {law:6.2f} {h:8.2f} {nn:7d} {h/law:9.2f}   {run}{cap}")
    rows.append((line, s, h, law, (G(r, 'subCpu') or 0) / 3600 if r is not None and G(r, 'id') in closed else None, ex))
print()
for line in ('r31', 'c9', 'c49', 'c9+c49'):
    v = [x for x in rows if x[0] in line.split('+')]
    if len(v) < 3: continue
    s = np.array([x[1] for x in v]); h = np.array([x[2] for x in v]); law = np.array([x[3] for x in v]); e = np.log(h / law)
    k, a = np.polyfit(s, np.log(h), 1)
    print(f"{line}: {len(v)} children ({sum(x[5] == 3 for x in v)} capped), {h.sum():.1f} CPU h here against the law's {law.sum():.1f} ({h.sum()/law.sum():.2f}); "
          f"ln(here / law) mean {e.mean():+.2f} (x{math.exp(e.mean()):.2f}), s.e. {e.std(ddof=1)/math.sqrt(len(v)):.2f}; own fit ln h = {a:.2f} + {k:.4f} s")
    w = [x for x in v if x[4]]
    if w:
        c = np.log(np.array([x[2] for x in w]) / np.array([x[4] for x in w]))
        print(f"  cold start here / the run's warm subtree, {len(w)} children: x{math.exp(c.mean()):.2f} (from x{math.exp(c.min()):.2f} to x{math.exp(c.max()):.2f}); "
              f"total {sum(x[2] for x in w):.2f} h here, {sum(x[4] for x in w):.2f} h in the run")
ctl = [x for x in rows if x[0] == 'r31']; alt = [x for x in rows if x[0] != 'r31']
if len(ctl) >= 3 and len(alt) >= 3:
    ec = np.log([x[2] / x[3] for x in ctl]); ea = np.log([x[2] / x[3] for x in alt])
    d = ea.mean() - ec.mean(); se = math.sqrt(ea.var(ddof=1) / len(ea) + ec.var(ddof=1) / len(ec))
    print(f"\nat equal shortfall a child of c9 or c49 costs x{math.exp(d):.2f} what a child of r31 does here (ln {d:+.2f}, s.e. {se:.2f}: from x{math.exp(d - 2*se):.2f} to x{math.exp(d + 2*se):.2f})")
