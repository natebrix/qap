#!/usr/bin/env python3
"""ucheap/truthla.py: the size score on the children's own QPB bounds, at the six line-truth nodes (README.md).

  python3 truthla.py

dumps/tho30_*.err, tho40_*.err (not *_root*): ls-LAcell row col bound for every free cell, the QPB bound rule 4's
lookahead gives the child. Each of the study's 8 lines scored from its children's shortfalls (sum, power, exp,
largest), against the truth (every child's subtree solved to the end) and against the same forms on the
children's cold first DNN solve.
"""
import glob, json, math, os, sys
sys.path.insert(0, "/home/user/wt-port/experiments/line_signals")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import signals as S, ucheap as UC

FORMS = [("sum (P = 1)", lambda s, smax, G: sum(s)),
         ("power 8", lambda s, smax, G: sum((x / smax) ** 8 for x in s if x > 0)),
         ("power 32", lambda s, smax, G: sum((x / smax) ** 32 for x in s if x > 0)),
         ("exp, L = 16", lambda s, smax, G: sum(math.exp(16 * (x - smax) / G) for x in s if x > 0)),
         ("exp, L = 32", lambda s, smax, G: sum(math.exp(32 * (x - smax) / G) for x in s if x > 0)),
         ("largest shortfall", lambda s, smax, G: max(s))]
agg = {(src, nm): [[], [], 0] for src in ("QPB", "DNN") for nm, _ in FORMS}
cells = []
for inst, pat in (("tho30", "tho30_solve/line_truth/results/*.json"), ("tho40", "tho40_line_truth/results/*.json")):
    for path in sorted(glob.glob(os.path.join(UC.EXP, pat))):
        rec = json.load(open(path))
        q, bound, cut = {}, None, None
        for ln in open(os.path.join(HERE, "dumps", inst + "_" + rec["node"].replace(":", "-").replace(",", "-") + ".err")):
            f = ln.split()
            if f[0] == "ls-LAcell": q[(int(f[1]), int(f[2]))] = float(f[3])
            elif f[0] == "ls-node:":
                bound = float(f[f.index("bound") + 1]); cut = float(f[f.index("inc") + 1]) - float(f[f.index("delta") + 1])
        G = cut - bound
        lines, truth, sq, sd = [], [], {}, {}
        for line, L in rec["lines"].items():
            kids = L["children"]
            if not all("truth" in k and "nodes" in k["truth"] for k in kids.values()): continue
            a, b = [], []
            for cell, k in kids.items():
                r, c = (int(x) for x in cell.split(":"))
                a.append(max(0.0, cut - q[(r, c)]))
                b.append(0.0 if k["score"].get("fathomed") else max(0.0, cut - k["score"]["lb"]))
                cells.append((a[-1], b[-1], k["truth"]["nodes"]))
            lines.append(line); truth.append(1 + sum(k["truth"]["nodes"] for k in kids.values())); sq[line], sd[line] = a, b
        best = min(truth)
        out = []
        for src, sv in (("QPB", sq), ("DNN", sd)):
            smax = max(max(v) for v in sv.values())
            Gs = G if src == "QPB" else smax * 4      # the DNN's exp on its own scale
            for nm, fn in FORMS:
                v = [-fn(sv[l], smax, Gs) for l in lines]
                top = max(range(len(lines)), key=lambda i: v[i])
                rho = S.spearman(v, [-t for t in truth])
                a = agg[(src, nm)]; a[0].append(rho); a[1].append(truth[top] / best); a[2] += truth[top] == best
                if src == "QPB" and nm in ("sum (P = 1)", "exp, L = 32"): out.append(f"{nm}: {lines[top]} {truth[top] / best:.2f}x ({rho:+.2f})")
        print(f"- {inst} {rec['node']}: best line {lines[truth.index(best)]} {best:,} nodes, worst {max(truth) / best:.2f}x, rule 9's "
              f"{rec['rule9']} {truth[lines.index(rec['rule9'])] / best:.2f}x; QPB children " + "; ".join(out))
print("\n| children bounded by | score | mean Spearman with the truth | best line first | its first line / the best, mean (worst) |\n|---|---|---|---|---|")
for (src, nm), (rho, loss, first) in agg.items():
    n = len(rho)
    print(f"| {src} | {nm} | {sum(rho) / n:+.2f} | {first}/{n} | {sum(loss) / n:.2f}x ({max(loss):.2f}x) |")
print(f"\nChild by child ({len(cells)} children): Spearman of the subtree's nodes with the child's QPB shortfall "
      f"{S.spearman([c[0] for c in cells], [c[2] for c in cells]):+.2f}, with its DNN first solve's "
      f"{S.spearman([c[1] for c in cells], [c[2] for c in cells]):+.2f}; QPB's with the DNN's {S.spearman([c[0] for c in cells], [c[1] for c in cells]):+.2f}.")
