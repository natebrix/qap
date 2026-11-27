"""Knuth's dives against Knuth's dives that solve the rest exactly once they reach depth D, at the same
CPU, on a complete node trace (the dives' variance below D is what carried tho40-est's estimate).

A dive from the root moves to one child a level (uniformly: QPB's gaps barely separate siblings under
the DNN bound) and, in the "full D" rule, on reaching a node at depth D adds that node's whole subtree,
weighted by the path's 1/p, and stops. Its cost is the CPU of the nodes it bounds: the path's, plus the
subtree's for "full D". Each design gets the same CPU (a share of the tree's, the root's own bound included once): as many
dives as fit.

    python3 -I sim_full.py TRACE[.gz] [--budget 0.05,0.2] [--reps 300] [--full 4,5]
"""
import argparse, os, random, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tho40b_analysis"))
from ntrace import load  # noqa: E402


def dive(t, rng, root, full):
    """one dive below the root (the stratum: bounded once, outside the dives, as latest/estimate.C does)"""
    ch = t.children[root]
    if t.outcome[root] != "B" or not ch:
        return 0.0, 0.0, 0.0
    i, w, n, c, cost = ch[rng.randrange(len(ch))], float(len(ch)), 0.0, 0.0, 0.0
    while True:
        if full is not None and t.depth[i] >= full:
            n += w * t.size[i]
            c += w * t.sub_cpu[i]
            cost += t.sub_cpu[i]
            return n, c, cost
        n += w
        c += w * t.cpu[i]
        cost += t.cpu[i]
        ch = t.children[i]
        if t.outcome[i] != "B" or not ch:
            return n, c, cost
        k = rng.randrange(len(ch))
        w *= len(ch)
        i = ch[k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trace")
    ap.add_argument("--budget", default="0.02,0.05,0.2")
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--full", default="4,5")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t = load(a.trace)
    root = t.roots[0]
    N, C = t.size[root], t.sub_cpu[root]
    print(f"{a.trace}: root depth {t.depth[root]}, {N:,} nodes, {C:,.0f} CPU s; by depth: " +
          ", ".join(f"{d}: {sum(1 for x in t.depth if x == d)}" for d in sorted(set(t.depth))))
    print("| rule | CPU budget (share of the tree) | dives | nodes: median (10%-90%) | CPU: median (10%-90%) | within 20% (CPU) |")
    print("|---|---|---|---|---|---|")
    rng = random.Random(a.seed)
    rules = [("Knuth, uniform", None)] + [(f"full below depth {d}", int(d)) for d in a.full.split(",")]
    for name, full in rules:
        for bud in (float(x) for x in a.budget.split(",")):
            rn, rc, nd = [], [], []
            for _ in range(a.reps):
                sn = sc = spent = 0.0
                k = 0
                while spent < bud * C - t.cpu[root] or k == 0:
                    n, c, cost = dive(t, rng, root, full)
                    sn += n; sc += c; spent += cost; k += 1
                rn.append((1 + sn / k) / N); rc.append((t.cpu[root] + sc / k) / C); nd.append(k)
            rn.sort(); rc.sort()
            p = lambda xs, f: xs[min(int(f * len(xs)), len(xs) - 1)]
            w20 = sum(0.8 <= x <= 1.2 for x in rc) / len(rc)
            print(f"| {name} | {bud:g} | {statistics.median(nd):,.0f} | {p(rn, .5):.2f} ({p(rn, .1):.2f}-{p(rn, .9):.2f}) | "
                  f"{p(rc, .5):.2f} ({p(rc, .1):.2f}-{p(rc, .9):.2f}) | {w20:.0%} |", flush=True)


if __name__ == "__main__":
    main()
