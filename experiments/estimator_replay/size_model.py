"""How a surviving child's subtree size depends on its shortfall, within its siblings: for every branched
node with >= 2 children that branch, regress log(subtree CPU) on x = shortfall / the median shortfall of
the branched nodes at the child's depth, demeaned by parent (the slope a survivor weight exp(b x) needs).
    python3 -I size_model.py TREE [TREE ...]"""
import math, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay
for path in sys.argv[1:]:
    T = replay.Tree(replay.load_shards(path) if os.path.isdir(path) else replay.load_trace(path))
    med = {}
    by = {}
    for i in range(len(T.depth)):
        if T.branched[i]: by.setdefault(T.depth[i], []).append(T.cutoff[i] - T.best[i])
    med = {d: statistics.median(v) for d, v in by.items()}
    xs, ys = [], []
    for i in range(len(T.depth)):
        if not T.branched[i]: continue
        surv = [j for j in T.kids[i] if T.branched[j] and T.C[j] > 0 and med.get(T.depth[j], 0) > 0]
        if len(surv) < 2: continue
        x = [(T.cutoff[j] - T.best[j]) / med[T.depth[j]] for j in surv]; y = [math.log(T.C[j]) for j in surv]
        mx, my = statistics.fmean(x), statistics.fmean(y)
        xs += [v - mx for v in x]; ys += [v - my for v in y]
    sxx = sum(v * v for v in xs); sxy = sum(a * b for a, b in zip(xs, ys)); syy = sum(v * v for v in ys)
    b = sxy / sxx if sxx else float("nan"); r2 = sxy * sxy / (sxx * syy) if sxx and syy else float("nan")
    print(f"{os.path.basename(path.rstrip('/'))[:40]:40s} survivors {len(xs):6d}  slope b {b:6.2f}  within-parent R^2 {r2:.2f}")
