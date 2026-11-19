#!/usr/bin/env python3
"""Step 2: can a child's first-solve shortfall be predicted out of sample, well enough to rank a parent's
children? Trained only on what a run would already have bounded:

  time   train on the children traced before the first X% of the trace's lines, test on the parents
         traced after it (all their children come later still)
  depth  train on the children at depth <= d, test on those below (d: about a third of the rows above)
  run    train on one run of an instance, test on another (tai50b-261116 <-> tai50b-r9: rule 2 below
         depth 6 against rule 9 everywhere; tho40ab-0 <-> tho40ab-32: the same subtree, two DNN state
         precisions; tho40's Nov 8 DNN trace -> tho40ab-32: another part of the tree)

  OMP_NUM_THREADS=2 python3 predict.py > predict_out.md

Scores (larger: more likely open): random; qpb = the child's QPB shortfall (cutoff - QPB bound, what rule 9
reads on its two lines); pc = a pseudocost table (row, column and pair effects with ridge shrinkage);
lin = ridge on the context and structural features; gbt = gradient-boosted trees on them, with the row
and column as categories; gbt-w = the same fitted to the target less its parent's mean (within-parent); "+q" = with the child's QPB bound as a feature. Metrics within each test parent
with 3 or more children (common.within_parent): AUC open vs fathomed, Spearman with the shortfall,
R-precision, and recall / precision of a shortlist of each parent's top quarter."""
import sys
import numpy as np, pandas as pd
from common import *

MODELS = [
    ("random", lambda tr: m_random(tr)),
    ("qpb", m_qpb),
    ("pc", lambda tr: m_pseudocost(tr)),
    ("lin", lambda tr: m_linear(tr, False)),
    ("lin+q", lambda tr: m_linear(tr, True)),
    ("gbt", lambda tr: m_gbt(tr, False)),
    ("gbt+q", lambda tr: m_gbt(tr, True)),
    ("gbt-w", lambda tr: m_gbt(tr, False, within=True)),
    ("gbt-w+q", lambda tr: m_gbt(tr, True, within=True)),
]


def run_split(title, tr, te, out):
    has_q = te.cqpb_s.notna().all() and tr.cqpb_s.notna().any()
    res = []
    for name, fit in MODELS:
        if ("q" in name.split("+")[-1] or name == "qpb") and not has_q:
            continue
        s = fit(tr)(te)
        r = within_parent(te, s)
        res.append((name, r))
    r0 = res[0][1]
    print(f"\n### {title}\n", file=out)
    print(f"train {len(tr):,} children; test {len(te):,} children of {r0['parents']:,} parents "
          f"({r0['scored']:,} with an open child)\n", file=out)
    print("| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |", file=out)
    print("|---|---|---|---|---|---|", file=out)
    for name, r in res:
        print(f"| {name} | {fmt(r['auc'])} | {fmt(r['spearman'])} | {fmt(r['rprec'])} | {fmt(r['recall25'])} | "
              f"{fmt(r['prec25'])} |", file=out)
    out.flush()
    return res


def main():
    out = sys.stdout
    summary = []
    print("# Step 2: out-of-sample prediction of a child's first solve (predict.py)", file=out)
    for run in ["tai60b-261117", "tai50b-261116", "tai50b-r9", "tho40ab-32", "tho40ab-0", "tho40-dnn",
                "tho30-dnn", "nug30-trace-a"]:
        d = load(run)
        print(f"\n## {run}", file=out)
        for frac in (0.25, 0.5):
            cut = d.cidx.quantile(frac)
            tr = d[d.cidx < cut]; te = d[d.pidx >= cut]
            if len(te) < 50 or len(tr) < 200: continue
            summary.append((run, f"time {frac:.0%}", run_split(f"time: train on the first {frac:.0%} of the trace", tr, te, out)))
        # depth: the smallest d leaving at least 30% of the children at or above it
        cum = d.depth.value_counts().sort_index().cumsum() / len(d)
        dd = int(cum[cum >= 0.3].index[0])
        tr = d[d.depth <= dd]; te = d[d.depth > dd]
        if len(te) >= 50:
            summary.append((run, f"depth <= {dd}", run_split(f"depth: train on depths <= {dd}, test below", tr, te, out)))
    print("\n## Across runs of one instance", file=out)
    for a, b in [("tai50b-261116", "tai50b-r9"), ("tai50b-r9", "tai50b-261116"), ("tho40ab-0", "tho40ab-32"),
                 ("tho40-dnn", "tho40ab-32")]:
        tr, te = load(a), load(b)
        if a == "tho40-dnn":
            te = te.copy(); te["cqpb_s"] = np.nan; te["pdnn_minus_qpb"] = np.nan
        summary.append((f"{a} -> {b}", "run", run_split(f"train {a}, test {b}", tr, te, out)))
    print("\n## Summary: AUC within parent (open vs fathomed), each score\n", file=out)
    names = [m for m, _ in MODELS]
    print("| test | split | " + " | ".join(names) + " |", file=out)
    print("|---|---|" + "---|" * len(names), file=out)
    for run, split, res in summary:
        r = dict(res)
        print(f"| {run} | {split} | " + " | ".join(fmt(r[m]["auc"]) if m in r else "-" for m in names) + " |", file=out)
    print("\n## Summary: Spearman within parent with the shortfall\n", file=out)
    print("| test | split | " + " | ".join(names) + " |", file=out)
    print("|---|---|" + "---|" * len(names), file=out)
    for run, split, res in summary:
        r = dict(res)
        print(f"| {run} | {split} | " + " | ".join(fmt(r[m]["spearman"]) if m in r else "-" for m in names) + " |", file=out)


if __name__ == "__main__":
    main()
