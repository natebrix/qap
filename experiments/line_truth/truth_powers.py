#!/usr/bin/env python3
"""truth_powers.py: the six truth nodes (tho30_solve/line_truth, tho40_line_truth: 8 lines each, every child's
subtree solved) scored from the children's DNN first-solve shortfalls with the exponents between the plain sum
and power 8 that branching_lab/ucheap/truthla.py did not try, and the size law's own form, sum of exp(k s) - 1."""
import glob, json, math, os
import numpy as np
E = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
def rank(a):
    a = np.asarray(a, float); r = np.empty(len(a)); o = np.argsort(a, kind='stable'); r[o] = np.arange(len(a))
    for v in set(a):
        m = a == v; r[m] = r[m].mean()
    return r
sp = lambda a, b: float(np.corrcoef(rank(a), rank(b))[0, 1])
FORMS = [(f"power {p:g}", (lambda p: lambda s, m: ((s / m) ** p).sum())(p)) for p in (1, 1.5, 2, 3, 4, 8)]
FORMS += [(f"exp(k s) - 1, k smax = {L:g}", (lambda L: lambda s, m: np.expm1(L * s / m).sum())(L)) for L in (2, 4, 8)]
FORMS += [("largest", lambda s, m: s.max())]
agg = {nm: [[], [], 0] for nm, _ in FORMS}
for inst, pat in (("tho30", "tho30_solve/line_truth/results/*.json"), ("tho40", "tho40_line_truth/results/*.json")):
    for path in sorted(glob.glob(os.path.join(E, pat))):
        rec = json.load(open(path)); cut = None
        for ln in open(os.path.join(E, "branching_lab/ucheap/dumps", inst + "_" + rec["node"].replace(":", "-").replace(",", "-") + ".err")):
            f = ln.split()
            if f and f[0] == "ls-node:": cut = float(f[f.index("inc") + 1]) - float(f[f.index("delta") + 1])
        lines, truth, sd = [], [], {}
        for line, L in rec["lines"].items():
            kids = L["children"]
            if not all("truth" in k and "nodes" in k["truth"] for k in kids.values()): continue
            sd[line] = np.array([0.0 if k["score"].get("fathomed") else max(0.0, cut - k["score"]["lb"]) for k in kids.values()])
            lines.append(line); truth.append(1 + sum(k["truth"]["nodes"] for k in kids.values()))
        best = min(truth); smax = max(v.max() for v in sd.values()); row = []
        for nm, fn in FORMS:
            v = [fn(sd[l], smax) for l in lines]; top = int(np.argmin(v))
            a = agg[nm]; a[0].append(sp(v, truth)); a[1].append(truth[top] / best); a[2] += truth[top] == best
            row.append(f"{lines[top]} {truth[top]/best:.2f}x")
        print(f"{inst} {rec['node']}: best {lines[truth.index(best)]} {best:,} nodes, worst {max(truth)/best:.2f}x, largest shortfall {smax:.0f} | first by form: " + ", ".join(row))
print("\n| score | mean Spearman with the truth | best line first | its first line / the best, mean (worst) |\n|---|---|---|---|")
for nm, (rho, loss, first) in agg.items():
    print(f"| {nm} | {np.mean(rho):+.2f} | {first}/{len(rho)} | {np.mean(loss):.2f}x ({max(loss):.2f}x) |")
