#!/usr/bin/env python3
"""harvest.py KICKS_RANDOM: the signal at tai50b's root with a second, incumbent-free search per cell, and
every local optimum found credited to all n cells it satisfies (T[r][v] = best cost seen with p[r] = v).
Per cell (r, v): (a) the incumbent with r:v forced by one exchange, descent, 30 kicks (the signal's base search);
(b) a random start with r:v forced, descent, KICKS_RANDOM kicks. Prints each line's predicted survivors
(T - best <= slack) under (a) alone, and under (a)+(b)+harvest, and checks c26 and r31 child by child
against the measured outcomes (root_children.tsv)."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/qap/experiments/branching_lab/landscape")
import lsq
K = int(sys.argv[1]); I = lsq.Inst("tai50b"); n = I.n; slack = 3372320.9; t0 = time.time()
Ta = np.full((n, n), np.inf); T = np.full((n, n), np.inf); idx = np.arange(n)
def credit(c, p):
    T[idx, p] = np.minimum(T[idx, p], c)
def ils(p, fix, free, kicks, rng, harvest):
    c, p = I.descend(p, fix)
    if harvest: credit(c, p)
    for _ in range(kicks):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix)
        if harvest: credit(cq, q)
        if cq < c: c, p = cq, q
    return c
rng = np.random.default_rng(7)
for r in range(n):
    free = [i for i in range(n) if i != r]
    for v in range(n):
        fix = {r: v}
        Ta[r, v] = ils(I.adjust(I.p0, fix), fix, free, 30, rng, True)
        T[r, v] = min(T[r, v], Ta[r, v])
        c = ils(I.adjust(rng.permutation(n).astype(np.int32), fix), fix, free, K, rng, True)
        T[r, v] = min(T[r, v], c)
    print(f"row {r} done {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
np.save(f"harvest_T_k{K}.npy", T); np.save(f"harvest_Ta_k{K}.npy", Ta)
xa = (Ta - I.best) / slack; x = (T - I.best) / slack
def lines(X):
    out = [(f"r{r}", int((X[r] <= 1).sum())) for r in range(n)] + [(f"c{v}", int((X[:, v] <= 1).sum())) for v in range(n)]
    return sorted(out, key=lambda t: t[1])
la, lh = lines(xa), lines(x)
print("base search alone, fewest survivors:", la[:8]); print("with random starts + harvest, fewest:", lh[:12])
rank = {nm: i + 1 for i, (nm, c) in enumerate(lh)}; cnt = dict(lh); cnta = dict(la)
truth = {}
for l in open("root_children.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "line": truth[(f[0], f[1])] = f[4]
for nm in ("c26", "r31"):
    kids = [(k, o) for (l, k), o in truth.items() if l == nm]
    tp = fp = fn = tn = 0
    for k, o in kids:
        r, v = map(int, k.split(":")); pred = x[r, v] <= 1; br = o == "B"
        tp += pred and br; fp += pred and not br; fn += (not pred) and br; tn += (not pred) and (not br)
    print(f"{nm}: base count {cnta[nm]}, harvested count {cnt[nm]} (rank {rank[nm]} of {2 * n}); measured branched {sum(o == 'B' for k, o in kids)}; tp {tp} fp {fp} fn {fn} tn {tn}")
print(f"CPU {time.time() - t0:.0f}s, K {K}")
