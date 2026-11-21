#!/usr/bin/env python3
"""siblings.py: sibling groups in the c47 descent (depths 8-18: runs/sc47, runs/d9_*, runs/down). For every parent with
two or more branching children: each child's line (rule 9's choice) and its own branching children. A group where the
same line is chosen at several siblings, each with many branching children, is "one bad choice repeated". Writes
out/siblings.txt."""
import glob, collections
CUT = 608215054.05
T = {}
for p in ['runs/sc47/nodes.txt'] + glob.glob('runs/d9_*/nodes.txt') + glob.glob('runs/down/d*/nodes/*.txt'):
    for l in open(p):
        f = l.split()
        if len(f) >= 12 and f[0] == 'node': T[f[2]] = f
kids = collections.defaultdict(list)
for k, f in T.items():
    if ',' in k: kids[k.rsplit(',', 1)[0]].append(f)
lb = lambda f: float(f[5] if f[5] != '-' else f[4] if f[4] != '-' else f[3])
rows = []
for par, K in kids.items():
    B = [f for f in K if f[7] == 'B']
    if len(B) < 2: continue
    rec = []
    for f in sorted(B, key=lambda f: lb(f)):
        nb = sum(g[7] == 'B' for g in kids.get(f[2], []))
        nk = len(kids.get(f[2], []))
        rec.append((f[2].rsplit(',', 1)[1], f[9], nb, nk, CUT - lb(f)))
    lines = collections.Counter(r[1] for r in rec)
    rows.append((int(T[par][1]) + 1 if par in T else par.count(',') + 2, par, rec, lines))
out = []
for dep, par, rec, lines in sorted(rows):
    same = [ln for ln, c in lines.items() if c >= 2]
    tot = sum(r[2] for r in rec)
    out.append(f"depth {dep} parent ...{','.join(par.split(',')[-3:])}: {len(rec)} branching siblings, their branching children {tot}; "
               f"lines repeated: {same or '-'}")
    for r in rec: out.append(f"    {r[0]:6s} line {r[1]:4s} shortfall {r[4]/1e3:6,.0f}k -> {r[2]} branching of {r[3]} bounded")
open('out/siblings.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
