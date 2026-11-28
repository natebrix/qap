#!/usr/bin/env python3
"""first_look.py SHARDDIR: the live run's depth-2 nodes under the depth-1 nodes fired in ../wil50_sb (all on r45):
each child's shortfall after its 150-iteration strong-branching solve there against the work under it in the run
so far (subtree CPU; open subtrees are still growing)."""
import json, sys, os, glob, math
import numpy as np
S = sys.argv[1]
d = json.load(open(os.path.join(S, '0.json'))); ix = {c: i for i, c in enumerate(d['columns'])}
rows = d['rows']; byid = {r[ix['id']]: r for r in rows}
kids = {}
for r in rows:
    kids.setdefault(r[ix['parent']], []).append(r)
def rank(a):
    r = np.empty(len(a)); r[np.argsort(a)] = np.arange(len(a)); return r
allx, ally = [], []
print("node   kids  hours  top5 share  Spearman(shortfall, CPU)  fit log CPU = a + k s (children over 0.05 h): k, factor per 100, rms factor")
for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../wil50_sb/runs/21-*.progress.txt'))):
    node = os.path.basename(f).split('.')[0].replace('-', ':')
    sb = {}
    for ln in open(f):
        t = ln.split(); a, b = t[4].split(':')
        if a == '45': sb[b] = max(float(t[8]) - float(t[10]), 0.0)
    par = [r for r in rows if r[ix['depth']] == 1 and r[ix['fix']] == node]
    if not par: print(f"{node:6s} not in the run's tree yet"); continue
    ch = kids.get(par[0][ix['id']], [])
    x, y, n = [], [], []
    for r in ch:
        c = r[ix['fix']].split(':')[1]
        if c not in sb: continue
        x.append(sb[c]); y.append(r[ix['subCpu']] / 3600); n.append(r[ix['size']])
    x, y = np.array(x), np.array(y)
    if len(x) < 5: continue
    rho = np.corrcoef(rank(x), rank(y))[0, 1]
    m = y > 0.05
    k, a = np.polyfit(x[m], np.log(y[m]), 1); res = np.log(y[m]) - (a + k * x[m])
    top5 = np.sort(y)[-5:].sum() / y.sum()
    print(f"{node:6s} {len(x):3d} {y.sum():7.1f}   {100*top5:4.0f}%        {rho:5.2f}                 k {k:.4f}  x{math.exp(100*k):6.1f}  x{math.exp(res.std()):.2f}  ({m.sum()} children)")
    allx += list(x[m]); ally += list(y[m])
allx, ally = np.array(allx), np.array(ally)
k, a = np.polyfit(allx, np.log(ally), 1); res = np.log(ally) - (a + k * allx)
print(f"all {len(allx)} children over 0.05 h: log CPU h = {a:.2f} + {k:.4f} s: x{math.exp(100*k):.1f} per 100 units of shortfall, rms a factor {math.exp(res.std()):.2f}, r2 {1 - res.var()/np.log(ally).var():.2f}")
for lo, hi in ((0, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 400)):
    m = (allx >= lo) & (allx < hi)
    if m.sum(): print(f"  shortfall {lo:3d}-{hi:3d}: {m.sum():3d} children, CPU h median {np.median(ally[m]):8.2f}, mean {ally[m].mean():8.2f}, max {ally[m].max():8.1f}")
