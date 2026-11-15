"""What a better bound would save, from a complete node trace
(docs/tho40b_results_plan.md, (b) 1).

A branched node whose best bound fell short of the cutoff by s would have been
fathomed by a bound d >= s higher, saving its whole subtree but itself. For each
d, the nodes saved are the subtrees of the topmost branched nodes with
shortfall <= d (a node under another such node is already counted); the CPU
saved is their subtree CPU less their own. --max-depth applies the better
bound only down to that depth (a costly bound used only near the top).

The near-miss list ranks branched nodes by nodes saved per unit of bound
(subtree size / shortfall): the frozen nodes worth a bound lab.

usage: python3 bound_value.py TRACE[.gz] [--max-depth D] [--near-misses K]
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ntrace import load  # noqa: E402

DELTAS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]


def saved(t, delta, max_depth=None):
    """(nodes saved, CPU saved, nodes saved by the depth of the new fathom)"""
    nodes = cpu = 0.0
    by_depth = collections.Counter()
    stack = [(r, False) for r in t.roots]
    while stack:
        i, covered = stack.pop()
        if covered:
            continue
        hit = t.outcome[i] == "B" and t.shortfall(i) <= delta and (max_depth is None or t.depth[i] <= max_depth)
        if hit:
            nodes += t.size[i] - 1
            cpu += t.sub_cpu[i] - t.cpu[i]
            by_depth[t.depth[i]] += t.size[i] - 1
            continue
        stack.extend((c, False) for c in t.children[i])
    return nodes, cpu, by_depth


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace")
    ap.add_argument("--max-depth", type=int)
    ap.add_argument("--near-misses", type=int, default=10)
    a = ap.parse_args()
    t = load(a.trace)
    N = len(t)
    total_cpu = sum(t.sub_cpu[r] for r in t.roots)
    cutoff = t.cutoff[t.roots[0]]
    branched = [i for i in range(N) if t.outcome[i] == "B"]
    print(f"{a.trace}: {N} nodes ({len(branched)} branched), {total_cpu:,.0f} CPU s, "
          f"roots at depth {sorted(set(t.depth[r] for r in t.roots))}, cutoff {cutoff:,.2f}"
          f"{'' if a.max_depth is None else f'; the better bound only at depths <= {a.max_depth}'}")
    print(f"{'d':>6} {'d % of cutoff':>14} {'nodes saved':>12} {'%':>6} {'CPU saved':>12} {'%':>6}  where (depth: nodes)")
    for d in DELTAS:
        n, c, by = saved(t, d, a.max_depth)
        where = ", ".join(f"{k}: {v:,.0f}" for k, v in sorted(by.items()) if v)
        print(f"{d:>6} {d / cutoff * 100:>13.3f}% {n:>12,.0f} {n / N * 100:>5.1f}% {c:>12,.0f} {c / total_cpu * 100:>5.1f}%  {where}")
    if a.near_misses:
        roots = set(t.roots)
        cand = sorted((i for i in branched if i not in roots), key=lambda i: -(t.size[i] - 1) / max(t.shortfall(i), 1e-9))
        print(f"\nnear misses: the {a.near_misses} branched nodes with the most nodes saved per unit of bound")
        print(f"{'depth':>5} {'shortfall':>10} {'gap %':>7} {'subtree':>8} {'CPU s':>9}  -F")
        for i in cand[:a.near_misses]:
            print(f"{t.depth[i]:>5} {t.shortfall(i):>10.1f} {t.gap_pct(i):>7.3f} {t.size[i]:>8,} {t.sub_cpu[i]:>9,.0f}  {t.fixes_flag(i)}")


if __name__ == "__main__":
    main()
