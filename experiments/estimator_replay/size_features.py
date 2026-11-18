"""Which features of a surviving child (known once it is bounded) predict its subtree, within its siblings:
within-parent least squares of log(subtree CPU) on
  x1 shortfall / the depth's median, x2 the cut phase's gain / the depth's median (0 without a phase),
  x3 log of the node's own CPU, x4 QPB's shortfall / the depth's median.
    python3 -I size_features.py TREE [TREE ...]"""
import math, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay

def fit(X, y):
    k = len(X[0]); A = [[sum(r[a] * r[b] for r in X) for b in range(k)] for a in range(k)]
    v = [sum(r[a] * t for r, t in zip(X, y)) for a in range(k)]
    # solve A w = v (Gauss)
    M = [A[a][:] + [v[a]] for a in range(k)]
    for c in range(k):
        p = max(range(c, k), key=lambda r: abs(M[r][c])); M[c], M[p] = M[p], M[c]
        if abs(M[c][c]) < 1e-12: continue
        for r in range(k):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y_ for x, y_ in zip(M[r], M[c])]
    w = [M[a][k] / M[a][a] if abs(M[a][a]) > 1e-12 else 0.0 for a in range(k)]
    pred = [sum(wi * xi for wi, xi in zip(w, r)) for r in X]
    sse = sum((t - p) ** 2 for t, p in zip(y, pred)); sst = sum(t * t for t in y)
    return w, 1 - sse / sst if sst else float("nan")

for path in sys.argv[1:]:
    T = replay.Tree(replay.load_shards(path) if os.path.isdir(path) else replay.load_trace(path))
    nodes = None
    X, y = [], []
    for i in range(len(T.depth)):
        if not T.branched[i]: continue
        surv = [j for j in T.kids[i] if T.branched[j] and T.C[j] > 0]
        if len(surv) < 2: continue
        rows = []
        for j in surv:
            m = T.med_short(T.depth[j])
            dnn = T.best[j]
            x1 = (T.cutoff[j] - T.best[j]) / m
            x4 = (T.cutoff[j] - T.qpb[j]) / m if not math.isnan(T.qpb[j]) else 0.0
            rows.append([x1, 0.0, math.log(max(T.cpu[j], 1e-3)), x4, math.log(T.C[j])])
        mu = [statistics.fmean(r[c] for r in rows) for c in range(5)]
        for r in rows:
            X.append([r[c] - mu[c] for c in range(4)]); y.append(r[4] - mu[4])
    name = os.path.basename(path.rstrip('/'))
    for cols, label in (([0], "shortfall"), ([2], "own CPU"), ([3], "QPB shortfall"), ([0, 2], "shortfall + own CPU"), ([0, 2, 3], "all three")):
        w, r2 = fit([[r[c] for c in cols] for r in X], y)
        print(f"{name[:30]:30s} {len(y):6d} survivors  {label:22s} R^2 {r2:.2f}  w {[round(v, 2) for v in w]}")
