"""Knuth's tree-size estimator, simulated on a complete node trace: how many
dives, and which way of choosing children, give an estimate within a factor of
the true size (docs/tho40b_results_plan.md, Better estimates).

A dive starts at the root and, at each branched node, moves to one child,
child c with probability p(c); it ends at a node that is not branched. With w
the product of 1/p along the path so far, the dive's estimate of the tree's
nodes is the sum of w over the nodes it visits (its CPU: the sum of w x CPU).
For a fixed incumbent the tree does not depend on the search order, and each
dive's estimate is unbiased; what the sampling rule changes is the variance.
On a complete trace every node's children are traced, so a simulated dive
sees what the solver's estimator (bnbdrv -e, latest/estimate.C) would, without
running the solver: it costs nothing once the trace exists.

Sampling rules:
  uniform   p(c) = 1/kids (Knuth 1975)
  gap       p(c) proportional to 10^(b x gap%(c)), the child's own bound's gap
            (b = 0.7 per %, the slope of log10 subtree size in the Nov 8 fits):
            children likely to grow big are sampled more and weighted less.
            Needs every child's bound at each step of the dive: a QPB bound
            each, not a DNN bound
  oracle    p(c) proportional to the child's true subtree size: zero variance,
            the limit of a perfect predictor (a check of the code)

usage: python3 knuth_sim.py TRACE[.gz] [--dives 30,100,300,1000] [--reps 200] [--b 0.7] [--seed 1]
"""
import argparse
import math
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ntrace import load  # noqa: E402


def qpb_gap_pct(t, i):
    v = t.qpb[i] if not math.isnan(t.qpb[i]) else t.lb(i)
    return max(0.0, t.cutoff[i] - v) / t.cutoff[i] * 100


def dive(t, rng, rule, b):
    i = t.roots[0]
    w = 1.0
    nodes = cpu = 0.0
    while True:
        nodes += w
        cpu += w * t.cpu[i]
        ch = t.children[i]
        if t.outcome[i] != "B" or not ch:
            return nodes, cpu
        if rule == "uniform":
            c = rng.randrange(len(ch))
            p = 1.0 / len(ch)
        else:
            ws = [10 ** (b * qpb_gap_pct(t, k)) for k in ch] if rule == "gap" else [t.size[k] for k in ch]
            tot = sum(ws)
            x = rng.random() * tot
            for c, wc in enumerate(ws):
                x -= wc
                if x <= 0:
                    break
            p = ws[c] / tot
        w /= p
        i = ch[c]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace")
    ap.add_argument("--dives", default="30,100,300,1000")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--b", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t = load(a.trace)
    if len(t.roots) != 1:
        sys.exit(f"{a.trace}: {len(t.roots)} roots; a dive needs one")
    N, C = len(t), t.sub_cpu[t.roots[0]]
    rng = random.Random(a.seed)
    print(f"{a.trace}: {N:,} nodes, {C:,.0f} CPU s. Each cell: the median and the 10th-90th percentile of the "
          f"estimate over the truth, over {a.reps} repetitions of n dives")
    print(f"{'rule':>8} {'dives':>6}  {'nodes: median (10%-90%)':>28}  {'CPU: median (10%-90%)':>28}  "
          f"{'one dive: largest share of the sum':>34}")
    for rule in ("uniform", "gap", "oracle"):
        for n in (int(x) for x in a.dives.split(",")):
            rn, rc, top = [], [], []
            for _ in range(a.reps):
                d = [dive(t, rng, rule, a.b) for _ in range(n)]
                rn.append(statistics.fmean(x for x, _ in d) / N)
                rc.append(statistics.fmean(y for _, y in d) / C)
                top.append(max(x for x, _ in d) / sum(x for x, _ in d))
            q = lambda v: statistics.quantiles(v, n=10)
            qn, qc = q(rn), q(rc)
            print(f"{rule:>8} {n:>6}  {statistics.median(rn):>8.2f} ({qn[0]:.2f}-{qn[-1]:.2f}){'':>9}"
                  f"  {statistics.median(rc):>8.2f} ({qc[0]:.2f}-{qc[-1]:.2f}){'':>9}  {statistics.median(top):>34.2f}")


if __name__ == "__main__":
    main()
