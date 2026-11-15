"""How well a run part-way through predicts its own final tree, scored on a
complete node trace (docs/tho40b_results_plan.md, (a) 1).

The trace's order is the order the run merged its nodes in, so its first
fraction f of lines is the tree the run had built at that point (what the live
tree showed). Two estimators predict the whole tree from that prefix:

  gap fit  (experiments/tho40_estimate/estimate.py's model) per depth,
           log10(subtree size) = a + b * gap% fitted on the prefix's complete
           branched subtrees, a pooled fit with depth as a covariate where a
           depth has fewer than 3; an incomplete node's missing children share
           its model size
  rates    the fathom rate and the children per branched node at each depth,
           S(d) = 1 + branched(d)/bounded(d) * kids(d) * S(d + 1); a missing child
           at depth d counts S(d)

Both keep the prefix's complete subtrees as they are. The score is the
predicted total over the actual.

usage: python3 estimator_score.py TRACE[.gz] [--fractions 0.05,0.1,...]
"""
import argparse
import collections
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ntrace import load, open_trace  # noqa: E402


def gap_fit(t):
    byd = collections.defaultdict(list)
    for i in range(len(t)):
        if t.outcome[i] == "B" and t.complete[i] and t.size[i] > 1:
            byd[t.depth[i]].append((t.gap_pct(i), math.log10(t.size[i])))
    model = {}
    for d, pts in byd.items():
        if len(pts) >= 3 and len(set(g for g, _ in pts)) > 1:
            b, a = np.polyfit([g for g, _ in pts], [y for _, y in pts], 1)
            model[d] = (a, b)
    xs = [(g, d, y) for d, pts in byd.items() for g, y in pts]
    coef = None
    if len(xs) >= 3:
        A = np.c_[np.ones(len(xs)), [g for g, _, _ in xs], [d for _, d, _ in xs]]
        coef, *_ = np.linalg.lstsq(A, np.array([y for _, _, y in xs]), rcond=None)

    def size(d, g):
        if d in model:
            a, b = model[d]
            return 10 ** (a + b * g)
        if coef is not None:
            return 10 ** (coef[0] + coef[1] * g + coef[2] * d)
        return 1.0
    return size


def rates(t):
    bounded, branched, kids = collections.Counter(), collections.Counter(), collections.Counter()
    for i in range(len(t)):
        d = t.depth[i]
        bounded[d] += 1
        if t.outcome[i] == "B":
            branched[d] += 1
            kids[d] += t.kids[i]
    S = {}
    for d in range(max(bounded) + 1, min(bounded) - 1, -1):
        if bounded[d] == 0:
            S[d] = 1.0
            continue
        k = kids[d] / branched[d] if branched[d] else 0.0
        S[d] = 1.0 + branched[d] / bounded[d] * k * S.get(d + 1, 1.0)
    return lambda d: S.get(d, 1.0)


def predict(t, child_est):
    """the prefix's whole-tree estimate: complete subtrees as they are, each
    missing child of an incomplete node at child_est(node)"""
    est = [0.0] * len(t)
    for i in sorted(range(len(t)), key=lambda i: -t.depth[i]):
        if t.complete[i]:
            est[i] = t.size[i]
            continue
        present = t.children[i]
        missing = max(0, t.kids[i] - len(present))
        est[i] = 1 + sum(est[c] for c in present) + missing * child_est(i)
    return sum(est[r] for r in t.roots)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace")
    ap.add_argument("--fractions", default="0.02,0.05,0.1,0.2,0.3,0.5,0.75")
    a = ap.parse_args()
    with open_trace(a.trace) as f:
        total_lines = sum(1 for ln in f if ln.startswith("node "))
    full = load(a.trace)
    N = len(full)
    print(f"{a.trace}: {N} nodes; prefix = the first fraction of the trace (merge order)")
    print(f"{'prefix':>7} {'nodes':>8} {'complete':>9} {'gap fit':>10} {'ratio':>6} {'rates':>10} {'ratio':>6}")
    for fr in (float(x) for x in a.fractions.split(",")):
        t = load(a.trace, limit=max(1, int(total_lines * fr)))
        size = gap_fit(t)
        g = predict(t, lambda i: max(0.0, size(t.depth[i], t.gap_pct(i)) - 1) / max(1, t.kids[i]))
        S = rates(t)
        r = predict(t, lambda i: S(t.depth[i] + 1))
        comp = sum(1 for i in range(len(t)) if t.complete[i] and t.outcome[i] == "B")
        print(f"{fr:>7.2f} {len(t):>8,} {comp:>9,} {g:>10,.0f} {g / N:>6.2f} {r:>10,.0f} {r / N:>6.2f}")


if __name__ == "__main__":
    main()
