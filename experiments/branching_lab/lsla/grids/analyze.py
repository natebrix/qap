#!/usr/bin/env python3
"""lsla/grids: rank every root row by the CAP 150 lookahead. Usage: analyze.py INST KIDS PRODROW > runs/INST/rank.md
Reads runs/INST/children.txt (one line per child solve) and runs/INST/p*/la.txt (the solver's own per-row scores).
A child the QPB test drops is not solved or logged: it counts as fathomed (KIDS - logged)."""
import glob, re, sys
inst, kids, prod = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rows = {}
for ln in open(f'runs/{inst}/children.txt'):
    m = re.match(r'child depth 0 cap 150 cell (\d+):(\d+) lb (\S+) shortfall (\S+) iters (\d+)', ln)
    if not m: continue
    r, c, sh, it = int(m[1]), int(m[2]), float(m[4]), int(m[5])
    rows.setdefault(r, {})[c] = (sh, it)
la = {}
for f in glob.glob(f'runs/{inst}/p*/la.txt'):
    for ln in open(f):
        for r, s, o in re.findall(r' r(\d+):[-\d.e]+:([\d.e+]+):(\d+)', ln.split(';', 1)[1] if ';' in ln else ''):
            la[int(r)] = (float(s), int(o))
out = []
for r, ch in rows.items():
    s = sum(v[0] for v in ch.values()); o = sum(v[0] > 0 for v in ch.values())
    fe = sum(v[0] == 0 for v in ch.values()) + kids - len(ch)
    out.append((s, o, fe, r, len(ch)))
out.sort()
best = out[0][0]
print(f'| rank | row | score (shortfall sum) | vs best | open | fathomed within cap | logged | solver la.txt |')
print('|---|---|---|---|---|---|---|---|')
for i, (s, o, fe, r, n) in enumerate(out, 1):
    mark = '**' if r == prod else ''
    l = la.get(r); lt = f'{l[0]:.1f}/{l[1]}' if l else '-'
    print(f'| {i} | {mark}r{r}{mark} | {s:.1f} | {s/best:.3f}x | {o} | {fe} | {n} | {lt} |')
done = [x for x in out if x[4] + 0 >= 0]
print(f'\nrows scored {len(out)}; production r{prod} rank ' +
      str(next((i for i, x in enumerate(out, 1) if x[3] == prod), '-')))
if len(out) > 1:
    print('spread: ' + ', '.join(f'#{k} {out[k-1][0]/best:.3f}x' for k in (2, 3, 5, 10, len(out)) if k <= len(out)))
