#!/usr/bin/env python3
"""Does the firing's score predict subtree work? Calibrate a size model on the estimate run's
shares of the tree by depth-1 node (wil50-est-261124, the recipe's lines), then apply it to the
lines strong branching scored at the same nodes.
Model: a child's subtree work = A * exp(k * shortfall) (shape of experiments/estimates/wil50/refit.txt,
k = 0.0205 there on frontier shortfalls); a fathomed child costs nothing beyond its own solve."""
import glob, math, os, sys
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs')
def lines(f):
    recs = []
    for ln in open(f):
        t = ln.split(); a, b = t[4].split(':')
        recs.append((a, b, max(float(t[8]) - float(t[10]), 0.0)))
    out = {}; rows = set(); i = 0
    while i < len(recs):
        a, b = recs[i][0], recs[i][1]
        kind = 'r' if i + 1 < len(recs) and recs[i + 1][0] == a else 'c' if i + 1 < len(recs) and recs[i + 1][1] == b else '?'
        key = a if kind != 'c' else b
        j = i
        while j < len(recs) and (recs[j][0] if kind != 'c' else recs[j][1]) == key: j += 1
        x = [r[2] for r in recs[i:j]]
        if kind == 'c': x += [r[2] for r in recs[:i] if r[1] == key and ('r' + r[0]) in rows]
        else: rows.add('r' + key)
        out[kind + key] = x; i = j
    return out
L = {os.path.basename(f).split('.')[0]: lines(f) for f in glob.glob(os.path.join(d, '21-*.progress.txt'))}
# the estimate's share of the rest (10,145 CPU h) and the line the recipe took there
est = {'21-0': ('r45', .38), '21-10': ('c0', .37), '21-1': ('c0', .08), '21-20': ('r45', .07), '21-11': ('r45', .04)}
def work(x, k): return sum(math.exp(k * s) for s in x if s > 0)
def fit(k):   # least squares in logs with a free scale
    pts = [(math.log(work(L[n][l], k)), math.log(sh)) for n, (l, sh) in est.items() if n in L and l in L[n]]
    c = sum(y - x for x, y in pts) / len(pts)
    return c, math.sqrt(sum((y - x - c) ** 2 for x, y in pts) / len(pts))
ks = [i / 1000 for i in range(0, 41)]
best = min(ks, key=lambda k: fit(k)[1])
for k in (0.0, 0.01, best, 0.0205, 0.03):
    c, rms = fit(k)
    print(f"\nk = {k:.4f}{' (best fit)' if k == best else ''}: rms log error {rms:.2f} (a factor {math.exp(rms):.2f})")
    for n in sorted(L):
        tag = est.get(n, ('', None))
        row = []
        for l, x in L[n].items():
            if len(x) < 29: continue
            p = math.exp(c) * work(x, k)
            row.append(f"{l} {100*p:5.1f}%{'*' if l == tag[0] else ''}")
        print(f"  {n:6s} est {'%4.0f%%' % (100*tag[1]) if tag[1] else '    -'} | " + '  '.join(row))
print("\n* the recipe's line (the estimate's share is for it). Shares are of the estimate's 10,145 CPU h.")
for n in sorted(L):
    for l, x in L[n].items():
        if len(x) < 29: continue
        s = sorted(x)
        print(f"{n:6s} {l:>4} n={len(x)} fathomed={sum(v==0 for v in x)} sum={sum(x):7.0f} median={s[len(s)//2]:5.0f} top5={' '.join('%.0f' % v for v in s[-5:])}")
