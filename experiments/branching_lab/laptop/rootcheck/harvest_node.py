#!/usr/bin/env python3
"""harvest_node.py INST FIXES SLACK KICKS_RANDOM TAG [LINE=truthB ...]: harvest.py at any node. FIXES r:v,...
(p[r] = v) or -. Per free cell (r, v): (a) the incumbent with the fixes and r:v forced, descent, 30 kicks;
(b) a random completion of the fixes with r:v forced, descent, KICKS_RANDOM kicks; every local optimum is
credited to every cell it satisfies. Prints each line's predicted survivors (excess <= SLACK) under (a) alone
and under (a)+(b)+harvest, and for the named lines the counts beside the measured branching children."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/qap/experiments/branching_lab/landscape")
import lsq
inst, fx, slack, K, tag = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
truth = dict(a.split("=") for a in sys.argv[6:])
I = lsq.Inst(inst); n = I.n; t0 = time.time()
fix0 = {} if fx == "-" else {int(a.split(":")[0]): int(a.split(":")[1]) for a in fx.split(",")}
rows = [r for r in range(n) if r not in fix0]; cols = [v for v in range(n) if v not in fix0.values()]
Ta = np.full((n, n), np.inf); T = np.full((n, n), np.inf); idx = np.arange(n)
def credit(c, p): T[idx, p] = np.minimum(T[idx, p], c)
def ils(p, fix, free, kicks, rng):
    c, p = I.descend(p, fix); credit(c, p)
    for _ in range(kicks):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix); credit(cq, q)
        if cq < c: c, p = cq, q
    return c
def randstart(rng, fix):
    p = np.empty(n, dtype=np.int32); fr = [r for r in range(n) if r not in fix]; fv = [v for v in range(n) if v not in fix.values()]
    for r, v in fix.items(): p[r] = v
    p[fr] = rng.permutation(fv); return p
rng = np.random.default_rng(7)
for r in rows:
    for v in cols:
        fix = dict(fix0); fix[r] = v; free = [i for i in rows if i != r]
        Ta[r, v] = ils(I.adjust(I.adjust(I.p0, fix0), fix), fix, free, 30, rng); T[r, v] = min(T[r, v], Ta[r, v])
        T[r, v] = min(T[r, v], ils(randstart(rng, fix), fix, free, K, rng))
np.save(f"harvest_{tag}_k{K}.npy", T)
def lines(X):
    out = [(f"r{r}", int(((X[r, cols] - I.best) <= slack).sum())) for r in rows] + [(f"c{v}", int(((X[rows, v] - I.best) <= slack).sum())) for v in cols]
    return sorted(out, key=lambda t: t[1])
la, lh = lines(Ta), lines(T); rank = {nm: i + 1 for i, (nm, c) in enumerate(lh)}; ranka = {nm: i + 1 for i, (nm, c) in enumerate(la)}
print(f"{tag}: {inst} m {len(rows)} slack {slack:.0f} K {K} CPU {time.time() - t0:.0f}s")
print(" base search alone, fewest survivors:", la[:8]); print(" with random starts + harvest, fewest:", lh[:10])
for nm, b in truth.items():
    print(f" {nm}: measured branching {b}; base count {dict(la)[nm]} (rank {ranka[nm]}), harvested count {dict(lh)[nm]} (rank {rank[nm]} of {len(lh)})")
