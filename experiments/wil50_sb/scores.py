#!/usr/bin/env python3
"""Per-line strong-branching scores from the firings' progress files (runs/*.progress.txt):
a line's children are consecutive; a row's share their first index, a column's their second.
Score as bnb.C: the sum over children of shortfall^P, P = 8 at the root (printed as its 8th root), 1 below."""
import glob, os, sys
d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs')
for f in sorted(glob.glob(os.path.join(d, '*.progress.txt'))):
    recs = []
    for ln in open(f):
        t = ln.split(); a, b = t[4].split(':')
        recs.append((a, b, max(float(t[8]) - float(t[10]), 0.0), float(t[t.index('cpu') + 1])))
    lines = []; rows = set()   # [kind, index, shortfalls, cpu]
    i = 0
    while i < len(recs):
        a, b = recs[i][0], recs[i][1]
        kind = 'r' if i + 1 < len(recs) and recs[i + 1][0] == a else 'c' if i + 1 < len(recs) and recs[i + 1][1] == b else '?'
        key = a if kind != 'c' else b
        j = i
        while j < len(recs) and (recs[j][0] if kind != 'c' else recs[j][1]) == key: j += 1
        x = [r[2] for r in recs[i:j]]
        if kind == 'c':   # a column's children on rows already tried were solved there and are not solved again
            x += [r[2] for r in recs[:i] if r[1] == key and ('r' + r[0]) in rows]
        else: rows.add('r' + key)
        lines.append((kind + key, x, sum(r[3] for r in recs[i:j])))
        i = j
    root = os.path.basename(f).startswith('root')
    print(f"{os.path.basename(f)}: {len(recs)} children, {sum(l[2] for l in lines)/3600:.1f} CPU h")
    for name, x, cpu in lines:
        sc = sum(s ** 8 for s in x) ** .125 if root else sum(x)
        print(f"  {name:>4} n={len(x):3d} score={sc:9.1f} fathomed={sum(s == 0 for s in x)}")
