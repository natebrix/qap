#!/usr/bin/env python3
"""rootfit.py CHILDREN.txt ROW=TRACE.gz [ROW=TRACE.gz ...] [--total ROW=NODES ...]
The root score as a size model, pooled over every traced tree. CHILDREN: the root lookahead's QAP_LA_PROGRESS
lines (every row's children and their shortfalls). ROW=TRACE: a whole tree branched on ROW at the root with its
node trace; the nodes under each depth-1 child are counted from it. Fits ln(nodes under a child) = a + b *
shortfall over all those children, then scores every row by the tree it predicts: 1 + the sum over its children
of exp(a + b * shortfall) (a child the lookahead fathoms, or QPB drops, counts 1).
Prints the fit; for each traced row the prediction from the fit WITHOUT that row's children (leave one row out)
beside the actual tree; and every row ranked. --total ROW=NODES adds an untraced tree's node count to the
comparison (not to the fit)."""
import sys, gzip, math, collections, statistics
args = sys.argv[1:]
ch = args[0]; traced = {}; totals = {}
it = iter(args[1:])
for a in it:
    if a == '--total': continue
    r, v = a.split('=', 1)
    if v.isdigit(): totals[r] = int(v)
    else: traced[r] = v
S = collections.defaultdict(dict)
for l in open(ch):
    t = l.split()
    if t[:1] == ['child'] and t[2] == '0':
        r, c = t[6].split(':'); S[r][c] = float(t[10])
nch = max(len(v) for v in S.values())
sub, actual = {}, dict(totals)
for row, tr in traced.items():
    cnt = collections.Counter(); n = 0
    for l in gzip.open(tr, 'rt'):
        t = l.split()
        if t[:1] != ['node']: continue
        n += 1
        if t[1] == '0': continue
        top = [x for x in t[2].split(',') if x.split(':')[0] == row]
        if top: cnt[top[0].split(':')[1]] += 1
    sub[row] = cnt; actual[row] = n
def fit(rows):
    xs, ys = [], []
    for r in rows:
        for c, s in S[r].items():
            if s > 0 and sub[r].get(c, 0) > 0: xs.append(s); ys.append(math.log(sub[r][c]))
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    res = [y - (a + b * x) for x, y in zip(xs, ys)]
    r2 = 1 - sum(e * e for e in res) / sum((y - my) ** 2 for y in ys)
    return a, b, r2, math.exp(statistics.pstdev(res)), len(xs), min(xs), max(xs)
def pred(r, a, b):
    return 1 + sum(math.exp(a + b * s) if s > 0 else 1.0 for s in S[r].values()) + (nch - len(S[r]))
a, b, r2, fac, n, lo, hi = fit(list(sub))
print(f'fit over {n} children of rows {" ".join("r" + r for r in sub)}: ln(nodes) = {a:.3f} + {b:.6f} * shortfall, '
      f'R^2 {r2:.2f}, residual a factor {fac:.2f}; shortfalls {lo:.0f} to {hi:.0f}')
print('\n| row | actual nodes | predicted, its own children left out of the fit | ratio | predicted, all in | largest child shortfall |')
print('|---|---|---|---|---|---|')
for r in sorted(actual, key=lambda r: actual[r]):
    if r in sub and len(sub) > 1:
        a1, b1 = fit([q for q in sub if q != r])[:2]; p1 = pred(r, a1, b1); note = ''
    else:
        p1 = pred(r, a, b); note = ' (untraced: not in any fit)' if r not in sub else ' (the only traced row: in sample)'
    print(f'| r{r} | {actual[r]:,} | {p1:,.0f}{note} | {p1 / actual[r]:.2f} | {pred(r, a, b):,.0f} | {max(S[r].values()):,.1f} |')
rows = sorted((pred(r, a, b), r) for r in S)
lin = {r: i + 1 for i, (_, r) in enumerate(sorted((sum(v.values()), r) for r, v in S.items()))}
if len(actual) > 2:
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); o = [0] * len(v)
        for i, j in enumerate(s): o[j] = i + 1
        return o
    def corr(x, y):
        mx, my = statistics.fmean(x), statistics.fmean(y)
        return sum((p - mx) * (q - my) for p, q in zip(x, y)) / (sum((p - mx) ** 2 for p in x) * sum((q - my) ** 2 for q in y)) ** 0.5
    ks = list(actual)
    print(f'\nover the {len(ks)} rows run: Spearman(predicted size, actual nodes) = {corr(rank([pred(r, a, b) for r in ks]), rank([actual[r] for r in ks])):+.2f}; '
          f'(sum of shortfalls, actual nodes) = {corr(rank([sum(S[r].values()) for r in ks]), rank([actual[r] for r in ks])):+.2f}')
print('\n| rank | row | predicted nodes | rank by the sum of shortfalls | largest child shortfall | actual |'); print('|---|---|---|---|---|---|')
for i, (p, r) in enumerate(rows):
    print(f'| {i + 1} | r{r} | {p:,.0f} | {lin[r]} | {max(S[r].values()):,.1f} | {actual.get(r, "")} |')
