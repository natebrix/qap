#!/usr/bin/env python3
"""The wil50 tree estimate of sampling/tree_estimate.py, read from the committed stage-B
table (sampling/results/stageB/summary.md) instead of the run directories, with more
measured subtrees ("anchors") than the record's three:

  refit.py [SUBTREE_DIR=SHORTFALL ...]
  refit.py --table STAGEB/summary.md --stage-a STAGEA --no-record [SUBTREE_DIR=SHORTFALL ...]   (another instance)

Each SUBTREE_DIR is a subtree.sh output whose node's converged shortfall is SHORTFALL
(the stage-B table's "converged short"); its DNN node count is read from its trace.
Always included: the record's anchors (node 10, 94 at 158.2; node 5, 1,198 at 315.7;
node 1, 13,245 [e] at 404.8) and the Nov 10 laptop subtree of 13:1,21:0,45:2 (3,374 at
284.5, zone 0.005, ../../wil50_feasibility/pipeline_ab/runs_zone005/).

Size models for a frontier node s units short (s from the table as tree_estimate.py
takes it: the converged shortfall for the 80 converged nodes, the give-up shortfall
otherwise, size 1 for a first-solve fathom):
  record anchors     log-linear between the record's anchors (tree_estimate.py)
  per-level          47^(s/130)
  fit                log(size) = a + b s, least squares over every measured subtree
  fit, mean          the fit times exp(sigma^2/2), sigma the fit's residual spread: the
                     mean of a subtree at s when sizes scatter log-normally about the fit
Totals: 651 nodes at depths 0-2 plus 14,192 frontier nodes times the mean subtree; CPU at
22.5 s per DNN node (the record's, cloud Xeon). 95% bootstrap over the 400 nodes."""
import math, os, random, re, sys
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.abspath(os.path.join(H, "../../.."))
TABLE = os.path.join(R, "experiments/wil50_feasibility/sampling/results/stageB/summary.md")
FRONT, TOP, CPU = 14192, 651, 22.5
RECORD = [(158.2, 94.0, "node 10"), (315.7, 1198.0, "node 5 [e]"), (404.8, 14262.0 * 1198.0 / 1290.0, "node 1 [e]")]
MEASURED = [(284.5, 3374.0, "13:1,21:0,45:2 (laptop, Nov 10)")]

def nodes_in(d):
    return sum(1 for ln in open(os.path.join(d, "trace.txt")) if ln.startswith("pid "))

def shortfalls():
    s = []
    for ln in open(TABLE):
        if not ln.startswith("| ") or ln.startswith("| node") or ln.startswith("|---"): continue
        c = [x.strip() for x in ln.strip().strip("|").split("|")]
        if len(c) < 15: continue
        try: short = float(c[4])
        except ValueError: continue
        if c[5] == "fathomed": s.append(-1.0); continue
        try: s.append(float(c[12]))
        except ValueError: s.append(short)
    return s

def interp(anchors):
    a = [(0.0, 1.0)] + sorted((x, n) for x, n, _ in anchors)
    def f(s):
        if s <= 0: return 1.0
        for (s0, n0), (s1, n1) in zip(a, a[1:]):
            if s <= s1: return math.exp(math.log(n0) + (math.log(n1) - math.log(n0)) * (s - s0) / (s1 - s0))
        s0, n0 = a[-1]; return n0 * 47 ** ((s - s0) / 130)
    return f

def fit(anchors):
    xs = [x for x, _, _ in anchors]; ys = [math.log(n) for _, n, _ in anchors]; k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs); a = my - b * mx
    res = [y - (a + b * x) for x, y in zip(xs, ys)]
    s2 = sum(r * r for r in res) / max(1, k - 2)
    return a, b, s2

def main():
    global TABLE, FRONT, TOP, CPU
    args, extra, record = sys.argv[1:], [], True
    while args and args[0].startswith("--"):
        o = args.pop(0)
        if o == "--table": TABLE = args.pop(0)
        elif o == "--stage-a":   # the frontier and the nodes above it from a stage-A run directory
            A = args.pop(0)
            FRONT = sum(1 for _ in open(os.path.join(A, "frontier.tsv"))) - 1
            TOP = sum(1 for t in os.listdir(os.path.join(A, "tasks")) if os.path.exists(os.path.join(A, "tasks", t, "done.json")))
        elif o == "--cpu": CPU = float(args.pop(0))
        elif o == "--no-record": record = False
        else: raise SystemExit(f"unknown option {o}")
    for arg in args:
        d, _, s = arg.rpartition("=")
        extra.append((float(s), float(nodes_in(d)), os.path.basename(d.rstrip("/"))))
    anchors = (RECORD + MEASURED if record else []) + extra
    if len(anchors) < 2:
        raise SystemExit("the fit needs two measured subtrees or more")
    a, b, s2 = fit(anchors)
    S = shortfalls()
    models = ([("record anchors", interp(RECORD))] if record else []) + [("per-level 47^(s/130)", lambda s: 47 ** (max(s, 0) / 130)),
              ("fit", lambda s: 1.0 if s <= 0 else max(1.0, math.exp(a + b * s))),
              ("fit, mean", lambda s: 1.0 if s <= 0 else max(1.0, math.exp(a + b * s + s2 / 2)))]
    print(f"frontier {FRONT:,} nodes, {TOP:,} above it; stage B: {len(S)} sampled, {sum(x < 0 for x in S)} fathomed in the first solve, the farthest {max(S):.1f} short")
    print("measured subtrees (shortfall, DNN nodes): " + "; ".join(f"{x:.1f} {n:,.0f} {l}" for x, n, l in sorted(anchors)))
    print(f"fit: log(size) = {a:.3f} + {b:.5f} s (x{math.exp(130 * b):.1f} per 130 units against 47), residual sigma {math.sqrt(s2):.2f}"
          f" (a factor {math.exp(math.sqrt(s2)):.1f}), mean factor {math.exp(s2 / 2):.2f}")
    print(f"| size model | mean subtree | total DNN nodes (95%) | CPU s at {CPU} s per node | core-hours |")
    print("|---|---|---|---|---|")
    rng = random.Random(1)
    for name, f in models:
        z = [f(s) for s in S]; m = sum(z) / len(z)
        boots = sorted(sum(z[rng.randrange(len(z))] for _ in z) / len(z) for _ in range(4000))
        tot = TOP + FRONT * m
        print(f"| {name} | {m:,.1f} | {tot:,.0f} ({TOP + FRONT * boots[100]:,.0f}-{TOP + FRONT * boots[3899]:,.0f}) | {tot * CPU:.2e} | {tot * CPU / 3600:,.0f} |")

if __name__ == "__main__":
    main()
