#!/usr/bin/env python3
"""cop5_ls.py Y.bin [SEEDS] [SEED]: local search for a violated Horn cut on a dumped iterate (see cop5.py).
From random 5-tuples and from 5-cycles of the strongest scaled entries, replace one index at a time by the index
that lowers min_d d'(H o Yh_S)d the most, until no replacement helps. Prints the least values found."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location('cop5', __import__('os').path.join(sys.path[0], 'cop5.py'))
src = open(spec.origin).read().replace('\nmain()\n', '\n'); ns = {}; exec(compile(src, 'cop5', 'exec'), ns)
cop_min, HORN = ns['cop_min'], ns['HORN']

path = sys.argv[1]; nseed = int(sys.argv[2]) if len(sys.argv) > 2 else 200; rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
Y = np.fromfile(path); N = int(round(len(Y) ** .5)); Y = Y.reshape(N, N); Y = (Y + Y.T) / 2; m = int(round((N - 1) ** .5))
d = np.diag(Y).copy(); live = np.where(d > 1e-7)[0]; s = np.sqrt(d[live]); Yh = np.maximum(Y[np.ix_(live, live)], 0) / np.outer(s, s); n = len(live)
off = Yh.copy(); np.fill_diagonal(off, 0)
def val(S): return cop_min(Yh[S[:, :, None], S[:, None, :]] * HORN)
def lab(S): return ' '.join('const' if live[t] == 0 else f"{(live[t]-1)//m}:{(live[t]-1)%m}" for t in S)
seeds = []
for _ in range(nseed):                                   # greedy cycles: follow strong entries
    S = [int(rng.integers(n))]
    while len(S) < 5:
        w = off[S[-1]].copy(); w[S] = 0; p = w ** 4; S.append(int(rng.choice(n, p=p / p.sum())))
    seeds.append(S)
results = []
for S in seeds:
    S = np.array(S); cur = val(S[None])[0]
    while True:
        best, bS = cur, None
        for pos in range(5):
            cand = np.repeat(S[None], n, 0); cand[:, pos] = np.arange(n)
            ok = np.array([len(set(r)) == 5 for r in cand]); v = np.where(ok, val(cand), np.inf); i = v.argmin()
            if v[i] < best - 1e-9: best, bS = v[i], cand[i].copy()
        if bS is None: break
        S, cur = bS, best
    results.append((cur, lab(S)))
results.sort(); seen = set()
print(f"{path}: {nseed} starts; least Horn values (negative = violated cut):")
for v, l in results:
    k = frozenset(l.split())
    if k in seen: continue
    seen.add(k); print(f"  {v:+.5f}  {l}")
    if len(seen) >= 8: break
print(f"  starts ending below 0: {sum(v < -1e-7 for v, _ in results)} of {len(results)}")
