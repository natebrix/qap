#!/usr/bin/env python3
"""price_check.py: does the pricing of lines_price.py hold where the truth is known? At the six truth nodes
(tho30_solve/line_truth, tho40_line_truth: 8 lines each, every child's subtree solved) fit the size law on one
line's children alone, ln(nodes under a child) = a + k shortfall, price every other line with it (sum over its
open children of exp(a + k s)), and compare price / the fitted line's price with the true ratio of the work.
The fitted line plays the part of the run's line in wil50. Rows (r) and columns (c) are both among the lines."""
import glob, json, math, os, sys
import numpy as np
E = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
allr = []; byk = {}
for inst, pat in (("tho30", "tho30_solve/line_truth/results/*.json"), ("tho40", "tho40_line_truth/results/*.json")):
    for path in sorted(glob.glob(os.path.join(E, pat))):
        rec = json.load(open(path)); cut = None
        for ln in open(os.path.join(E, "branching_lab/ucheap/dumps", inst + "_" + rec["node"].replace(":", "-").replace(",", "-") + ".err")):
            f = ln.split()
            if f and f[0] == "ls-node:": cut = float(f[f.index("inc") + 1]) - float(f[f.index("delta") + 1])
        L = {}
        for line, d in rec["lines"].items():
            ch = d["children"]
            if not all("truth" in c and "nodes" in c["truth"] for c in ch.values()): continue
            s = np.array([0.0 if c["score"].get("fathomed") else max(0.0, cut - c["score"]["lb"]) for c in ch.values()])
            w = np.array([float(c["truth"]["nodes"]) for c in ch.values()])
            L[line] = (s, w)
        tot = {l: w.sum() for l, (s, w) in L.items()}
        print(f"\n{inst} {rec['node']}: {len(L)} lines, true work {min(tot.values()):,.0f} to {max(tot.values()):,.0f} nodes")
        print("  fitted on  children  k        rms factor | predicted / true ratio of each other line's work to the fitted line's: median, worst | lines ordered right vs the fitted line")
        for base, (s, w) in L.items():
            m = (s > 0) & (w > 1)
            if m.sum() < 5: continue
            k, a = np.polyfit(s[m], np.log(w[m]), 1); res = np.log(w[m]) - (a + k * s[m])
            price = lambda v: float(np.exp(a + k * v[v > 0]).sum())
            pb = price(s); errs = []; right = 0; n = 0
            for l, (s2, w2) in L.items():
                if l == base: continue
                pred = price(s2) / pb; true = tot[l] / tot[base]; errs.append(pred / true); n += 1
                right += (pred > 1) == (true > 1)
                allr.append((inst, rec['node'], base, l, pred, true))
            e = np.array(errs); le = np.abs(np.log(e))
            print(f"  {base:9s} {int(m.sum()):6d}   {k:.4f}   {math.exp(res.std()):6.2f}     | {np.median(e):5.2f}, {e[np.argmax(le)]:5.2f} | {right}/{n}")
A = np.array([(p, t) for *_, p, t in allr]); lr = np.log(A[:, 0] / A[:, 1])
print(f"\nall {len(A)} (fitted line, other line) pairs: predicted / true ratio: median {math.exp(np.median(lr)):.2f}, rms factor {math.exp(math.sqrt((lr**2).mean())):.2f}, "
      f"within 1.25x {np.mean(np.abs(lr) < math.log(1.25)):.0%}, within 1.5x {np.mean(np.abs(lr) < math.log(1.5)):.0%}")
for lo, hi, nm in ((0, 0.7, "predicted under 0.7 of the fitted line"), (0.7, 0.9, "predicted 0.7 to 0.9"), (0.9, 1.1, "predicted 0.9 to 1.1"), (1.1, 99, "predicted over 1.1")):
    m = (A[:, 0] >= lo) & (A[:, 0] < hi)
    if m.sum(): print(f"  {nm}: {m.sum()} pairs, true ratio median {np.median(A[m, 1]):.2f}, range {A[m, 1].min():.2f} to {A[m, 1].max():.2f}, truly cheaper in {np.mean(A[m, 1] < 1):.0%}")
