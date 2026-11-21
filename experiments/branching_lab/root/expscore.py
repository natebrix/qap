#!/usr/bin/env python3
"""expscore.py CHILDREN.txt TRACE.gz ROW: a root score that predicts tree size. From a whole tree branched on ROW
at the root (TRACE: nodes.txt.gz), the nodes under each depth-1 child; with the root lookahead's shortfall for
the same children (CHILDREN: QAP_LA_PROGRESS lines) fit ln(nodes) = a + b * shortfall; then for every row
scored at the root, the predicted tree = sum over its children of exp(a + b * shortfall) (a child that fathoms
in the lookahead counts 1). Prints the fit, and the rows ranked by the prediction beside the linear score (the
sum of shortfalls). A child missing from CHILDREN (QPB dropped it) counts 1."""
import sys, gzip, math, collections, statistics
ch, tr, row = sys.argv[1], sys.argv[2], sys.argv[3]
S = collections.defaultdict(dict)
for l in open(ch):
    t = l.split()
    if t[:1] == ['child'] and t[2] == '0':
        r, c = t[6].split(':'); S[r][c] = float(t[10])
sub = collections.Counter()
for l in gzip.open(tr, 'rt'):
    t = l.split()
    if t[:1] != ['node'] or t[1] == '0': continue
    cells = t[2].split(',')
    top = [x for x in cells if x.split(':')[0] == row]
    if top: sub[top[0].split(':')[1]] += 1
xs, ys = [], []
print(f'row r{row}: depth-1 child, lookahead shortfall, nodes under it in the tree')
for c in sorted(sub, key=lambda c: -sub[c]):
    s = S[row].get(c)
    print(f'  {row}:{c}  {s if s is not None else "-":>10}  {sub[c]}')
    if s is not None and s > 0: xs.append(s); ys.append(math.log(sub[c]))
mx, my = statistics.fmean(xs), statistics.fmean(ys)
b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
a = my - b * mx
res = [y - (a + b * x) for x, y in zip(xs, ys)]
r = 1 - sum(e * e for e in res) / sum((y - my) ** 2 for y in ys)
print(f'fit over {len(xs)} children: ln(nodes) = {a:.3f} + {b:.6f} * shortfall, R^2 {r:.2f}, residual sd {statistics.pstdev(res):.2f} (a factor {math.exp(statistics.pstdev(res)):.2f})')
nch = max(len(v) for v in S.values())
rows = []
for rr, v in S.items():
    pred = sum(math.exp(a + b * s) if s > 0 else 1.0 for s in v.values()) + (nch - len(v))
    rows.append((pred, rr, sum(v.values()), max(v.values())))
lin = {rr: i + 1 for i, (_, rr) in enumerate(sorted((x[2], x[1]) for x in rows))}
print('| rank | row | predicted nodes | linear score | its rank | largest child shortfall |'); print('|---|---|---|---|---|---|')
for i, (p, rr, l, m) in enumerate(sorted(rows)):
    print(f'| {i + 1} | r{rr} | {p:,.0f} | {l:,.1f} | {lin[rr]} | {m:,.1f} |')
