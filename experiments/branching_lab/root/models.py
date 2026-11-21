#!/usr/bin/env python3
"""models.py CHILDREN.txt ROW=TRACE.gz ... [--total ROW=NODES ...]: which form of the size model? Over the traced
rows' depth-1 children (shortfall s after the lookahead, nodes n under the child in the tree):
  exp:   ln n = a + b s          power: ln n = a + p ln s          quad: ln n = a + b s + c s^2
For each: the fit, and each row's tree predicted with that row's children left out, against the actual tree
(untraced --total rows: predicted from all traced rows). Children with s = 0 or 1 node are left out of the fits
and count 1 node in a prediction."""
import sys, gzip, math, collections, statistics
args = sys.argv[1:]; ch = args[0]; traced, totals = {}, {}
for a in args[1:]:
    if a == '--total': continue
    r, v = a.split('=', 1)
    (totals if v.isdigit() else traced)[r] = int(v) if v.isdigit() else v
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
def lsq(X, y):          # normal equations, tiny
    k = len(X[0]); A = [[sum(r[i] * r[j] for r in X) for j in range(k)] for i in range(k)]; v = [sum(r[i] * t for r, t in zip(X, y)) for i in range(k)]
    for i in range(k):
        p = max(range(i, k), key=lambda q: abs(A[q][i])); A[i], A[p] = A[p], A[i]; v[i], v[p] = v[p], v[i]
        for q in range(i + 1, k):
            f = A[q][i] / A[i][i]
            for j in range(i, k): A[q][j] -= f * A[i][j]
            v[q] -= f * v[i]
    w = [0.0] * k
    for i in reversed(range(k)): w[i] = (v[i] - sum(A[i][j] * w[j] for j in range(i + 1, k))) / A[i][i]
    return w
FEAT = {'exp': lambda s: [1.0, s], 'power': lambda s: [1.0, math.log(s)], 'quad': lambda s: [1.0, s, s * s]}
def pts(rows):
    return [(S[r][c], math.log(sub[r][c])) for r in rows for c in S[r] if S[r][c] > 0 and sub[r].get(c, 0) > 1]
def fit(model, rows):
    P = pts(rows); f = FEAT[model]; w = lsq([f(s) for s, _ in P], [y for _, y in P])
    res = [y - sum(a * b for a, b in zip(w, f(s))) for s, y in P]; my = statistics.fmean(y for _, y in P)
    return w, 1 - sum(e * e for e in res) / sum((y - my) ** 2 for _, y in P), math.exp(statistics.pstdev(res)), len(P)
def pred(model, w, r):
    f = FEAT[model]
    return 1 + sum(max(1.0, math.exp(sum(a * b for a, b in zip(w, f(s))))) if s > 0 else 1.0 for s in S[r].values()) + (nch - len(S[r]))
for model in FEAT:
    w, r2, fac, n = fit(model, list(sub))
    print(f'{model}: coefficients {" ".join(f"{x:.5g}" for x in w)}; R^2 {r2:.2f}, residual a factor {fac:.2f}, {n} children')
    errs = []
    for r in sorted(actual, key=lambda r: actual[r]):
        w1 = fit(model, [q for q in sub if q != r])[0] if (r in sub and len(sub) > 1) else w
        p = pred(model, w1, r); errs.append(abs(math.log(p / actual[r])))
        print(f'   r{r}: actual {actual[r]:,}, predicted {"(left out) " if r in sub else "(untraced) "}{p:,.0f}, ratio {p / actual[r]:.2f}')
    print(f'   mean |ln ratio| {statistics.fmean(errs):.2f} (a factor {math.exp(statistics.fmean(errs)):.2f})')
    rows = sorted((pred(model, w, r), r) for r in S)
    print('   ranking: ' + ' '.join(f'r{r}:{p:,.0f}' for p, r in rows[:6]) + ' ... ' + ' '.join(f'r{r}:{p:,.0f}' for p, r in rows[-3:]))
