import sys, random, statistics
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import replay
T = replay.Tree(replay.load_trace(sys.argv[1]) if sys.argv[1].endswith(".gz") else replay.load_shards(sys.argv[1]))
N, C = T.N[T.root], T.C[T.root]
rng = random.Random(5)
for mode, q in (("surv", 0), ("surv", 1), ("knuth", 0), ("knuth", 1.5)):
    xs = []; costs = []; depths = []
    for _ in range(int(sys.argv[2])):
        paid = set()
        en, ec, c, last = replay.dive(T, T.root, 1.0, mode, q, rng, paid)
        xs.append(ec / C); costs.append(c / C); depths.append(last)
    xs.sort()
    print(mode, q, "mean", round(statistics.fmean(xs), 3), "median", round(xs[len(xs)//2], 4), "p90", round(xs[int(.9*len(xs))], 3), "p99", round(xs[int(.99*len(xs))], 2), "max", round(xs[-1], 1), "cost/dive", round(statistics.fmean(costs), 4), "depth mean", round(statistics.fmean(depths), 1), "max", max(depths))
