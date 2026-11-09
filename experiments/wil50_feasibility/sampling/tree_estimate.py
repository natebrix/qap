#!/usr/bin/env python3
"""The tree estimate from stages A and B (README.md, "Stage B on wil50t"):

  total DNN nodes = (the nodes at depths 0-2) + (frontier size) x (mean subtree of a
  frontier node), the mean over stage B's 400 sampled nodes of size(shortfall),
  per depth-3 orbit = total / 29,580.

size(s) is the DNN nodes in the subtree below a depth-3 node s short at
convergence (the converged shortfall for the 80 converged nodes, the give-up
shortfall, an upper bound on it, for the others; 1 for a node that fathoms in
its first solve), log-linear between measured subtrees in the production
configuration or the nearest one measured:

  s = 0        1 node (it fathoms)
  s = 158.2    94 DNN nodes: node 10, cpcutt (pipeline_ab/runs_cpt/cpcutt, Nov 7)
  s = 315.7    1,198: node 5, cpcutt, matched probe estimate (tail_probe/runs_node5_cpt)
  s = 404.8    13,250 [e]: node 1, 14,262 under the Nov 5 `full` setting
               (tail_probe/runs/full), times node 5's cpcutt/full ratio 1,198/1,290
  beyond       47^(s/130) per unit of shortfall from the last anchor

and, as the alternative, the per-level model alone, size(s) = 47^(s/130).
Shortfalls are uncut (1,500 iterations for the anchors, 4,000 here). A 95%
bootstrap over the 400 nodes gives the interval; the anchors' own error (node 5's
SE is about 35%) is not in it. Writes results/stageB/tree_estimate.md.

  tree_estimate.py [--stage-a runs/stageA] [--stage-b runs/stageB]
"""
import argparse
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import D  # noqa: E402

ANCHORS = [(0.0, 1.0), (158.2, 94.0), (315.7, 1198.0), (404.8, 14262.0 * 1198.0 / 1290.0)]
ORBITS = 29580
CPU_PER_NODE = (2021 / 94 + 28095 / 1198) / 2      # node 10 and node 5 under cpcutt, CPU s per DNN node


def size_anchor(s):
    if s <= 0:
        return 1.0
    for (s0, n0), (s1, n1) in zip(ANCHORS, ANCHORS[1:]):
        if s <= s1:
            return math.exp(math.log(n0) + (math.log(n1) - math.log(n0)) * (s - s0) / (s1 - s0))
    s0, n0 = ANCHORS[-1]
    return n0 * 47 ** ((s - s0) / 130)


def size_lines(s):
    return 47 ** (max(s, 0.0) / 130)


def klass(s):
    if s <= 0:
        return "fathomed"
    return "near (< 150)" if s < 150 else "150-210" if s < 210 else "210-300" if s < 300 else "tail (>= 300)"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage-a", default=os.path.join(D, "runs", "stageA"))
    ap.add_argument("--stage-b", default=os.path.join(D, "runs", "stageB"))
    ap.add_argument("--converge", type=int, default=80)
    ap.add_argument("--out", default=os.path.join(D, "results", "stageB", "tree_estimate.md"))
    a = ap.parse_args()
    frontier = sum(1 for _ in open(os.path.join(a.stage_a, "frontier.tsv"))) - 1
    top = sum(1 for t in os.listdir(os.path.join(a.stage_a, "tasks"))
              if os.path.exists(os.path.join(a.stage_a, "tasks", t, "done.json")))
    sample = [l.rstrip("\n").split("\t") for l in open(os.path.join(a.stage_b, "sample.tsv")) if l[0].isdigit()]
    load = lambda key: (json.load(open(p)) if os.path.exists(p := os.path.join(a.stage_b, "tasks", key, "result.json"))
                        else None)
    rows = []
    for r in sample:
        sp = r[4]
        k = sp.replace(":", "_").replace(",", "-")
        g = load("g_" + k)
        if g is None:
            continue
        t = g["trace"]
        rows.append({"spec": sp, "g": t, "v": None, "k": k})
    miss = sorted((x for x in rows if x["g"] and not x["g"]["fathomed"]), key=lambda x: -x["g"]["short"])
    far = {x["spec"] for x in miss[:a.converge]}
    for x in rows:
        if x["spec"] in far:
            v = load("v_" + x["k"])
            x["v"] = v["trace"] if v else None
    nconv = sum(1 for x in rows if x["v"])
    for x in rows:
        if x["g"] is None or x["g"]["fathomed"]:
            x["s"] = -1.0
        elif x["v"] is not None:
            x["s"] = -1.0 if x["v"]["fathomed"] else x["v"]["short"]
        else:
            x["s"] = x["g"]["short"]
    n = len(rows)
    L = [f"# Tree estimate from stages A and B ({n} of {len(sample)} sampled frontier nodes; {nconv} converged)", "",
         f"Frontier {frontier:,} depth-3 nodes (stage A), {top} nodes bounded at depths 0-2; {ORBITS:,} depth-3 orbits "
         f"(Burnside), so {frontier / ORBITS:.3f} frontier nodes per orbit.", ""]
    rng = random.Random(1)
    res = {}
    for name, fn in (("anchors (central)", size_anchor), ("per-level model 47^(s/130)", size_lines)):
        sizes = [fn(x["s"]) for x in rows]
        mean = sum(sizes) / n
        tail = sum(z for z, x in zip(sizes, rows) if x["s"] >= 300)
        boots, tshares = [], []
        for _ in range(4000):
            idx = [rng.randrange(n) for _ in range(n)]
            m = sum(sizes[i] for i in idx) / n
            boots.append(m)
            tot = sum(sizes[i] for i in idx)
            tshares.append(sum(sizes[i] for i in idx if rows[i]["s"] >= 300) / tot if tot else 0)
        boots.sort()
        tshares.sort()
        lo, hi = boots[100], boots[3899]
        total = top + frontier * mean
        res[name] = (mean, lo, hi, total, tail / sum(sizes), tshares[100], tshares[3899])
    L += ["| size model | mean subtree per frontier node (95%) | total DNN nodes (95%) | per depth-3 orbit | "
          "tail (>= 300) share (95%) | CPU s at 22 s per DNN node [e] |", "|---|---|---|---|---|---|"]
    for name, (mean, lo, hi, total, ts, tlo, thi) in res.items():
        L.append(f"| {name} | {mean:,.1f} ({lo:,.1f}-{hi:,.1f}) | {total:,.0f} ({top + frontier * lo:,.0f}-"
                 f"{top + frontier * hi:,.0f}) | {total / ORBITS:,.1f} | {ts:.2f} ({tlo:.2f}-{thi:.2f}) | "
                 f"{total * CPU_PER_NODE:.2e} |")
    # by class
    L += ["", "By class (anchors):", "", "| class | sampled nodes | shortfall range | mean subtree | share of the total |",
          "|---|---|---|---|---|"]
    sizes = {id(x): size_anchor(x["s"]) for x in rows}
    tot = sum(sizes.values())
    for c in ("fathomed", "near (< 150)", "150-210", "210-300", "tail (>= 300)"):
        xs = [x for x in rows if klass(x["s"]) == c]
        if not xs:
            L.append(f"| {c} | 0 | | | 0 |")
            continue
        ss = [x["s"] for x in xs]
        L.append(f"| {c} | {len(xs)} | {min(ss):.1f}-{max(ss):.1f} | {sum(sizes[id(x)] for x in xs) / len(xs):,.1f} | "
                 f"{sum(sizes[id(x)] for x in xs) / tot:.3f} |")
    L += ["", "Anchors (shortfall, DNN nodes): " + ", ".join(f"({s:.1f}, {v:,.0f})" for s, v in ANCHORS) +
          f"; CPU per DNN node {CPU_PER_NODE:.1f} s (node 10 and node 5 under cpcutt).",
          "Against roadmap/tree_size.md: 1.0e4-6e4 DNN nodes per depth-3 orbit, 3e8-1.8e9 in all, 85-91% from the "
          "tail (the Oct 29 settings, the random-line 16-node sample)."]
    txt = "\n".join(L) + "\n"
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
