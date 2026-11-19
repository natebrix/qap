#!/usr/bin/env python3
"""Step 1: describe the child tables and how much of the variance each feature explains, out of fold
(5 folds grouped by parent, so a parent's children are never split between fit and test).

  OMP_NUM_THREADS=2 python3 describe.py > describe_out.md

Targets: ylog = log1p(1e4 x shortfall / |cutoff|) (0 when the first solve fathoms), open (the first solve
did not fathom), gainrel = (min(child DNN, cutoff) - parent bound) / |cutoff|. R2 total: over all children;
R2 within: of the target less its parent's mean, against the prediction less its parent's mean (what
ranking a parent's children needs: a parent-level feature explains nothing within)."""
import sys
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from common import *

RUNS = ["tai60b-261117", "tai50b-261116", "tai50b-r9", "tho40ab-32", "tho40-dnn", "tho30-dnn", "nug30-trace-a"]


def r2s(d, y, pred):
    yt = d[y].values
    tot = 1 - ((yt - pred) ** 2).sum() / ((yt - yt.mean()) ** 2).sum()
    g = d.pkey.values
    ym = pd.Series(yt).groupby(g).transform("mean").values
    pm = pd.Series(pred).groupby(g).transform("mean").values
    sw = ((yt - ym) ** 2).sum()
    win = 1 - (((yt - ym) - (pred - pm)) ** 2).sum() / sw if sw > 0 else np.nan
    return tot, win


def oof(d, fit, y):
    pred = np.zeros(len(d))
    for tr, te in GroupKFold(5).split(d, groups=d.pkey):
        pred[te] = fit(d.iloc[tr], y)(d.iloc[te])
    return pred


def main():
    print("# Step 1: the child tables (describe.py)\n")
    print("| run | children | parents | open (first solve missed) | of the open: C / B | depths | lines: rows |")
    print("|---|---|---|---|---|---|---|")
    ds = {}
    for run in RUNS:
        d = load(run); ds[run] = d
        o = d[d.open == 1]
        print(f"| {run} | {len(d):,} | {d.pkey.nunique():,} | {d.open.mean():.3f} | "
              f"{(o.outcome == 'C').mean():.2f} / {(o.outcome != 'C').mean():.2f} | {d.depth.min()}-{d.depth.max()} | "
              f"{fmt(d.line_is_row.mean())} |")
    print("\nShare of the variance between parents (a parent fixed effect), the ceiling of any parent-level feature:\n")
    print("| run | ylog | open | gainrel |")
    print("|---|---|---|---|")
    for run, d in ds.items():
        row = []
        for y in ("ylog", "open", "gainrel"):
            m = d.groupby("pkey")[y].transform("mean")
            row.append(1 - ((d[y] - m) ** 2).sum() / ((d[y] - d[y].mean()) ** 2).sum())
        print(f"| {run} | " + " | ".join(f"{x:.2f}" for x in row) + " |")

    sets = [
        ("row", lambda tr, y: m_pseudocost(tr, ("row",), target=y)),
        ("column", lambda tr, y: m_pseudocost(tr, ("col",), target=y)),
        ("row + column", lambda tr, y: m_pseudocost(tr, ("row", "col"), target=y)),
        ("row + column + pair (pseudocost)", lambda tr, y: m_pseudocost(tr, ("row", "col", "pair"), target=y)),
        ("depth", lambda tr, y: m_gbt(tr, target=y, fs=["depth"])),
        ("parent's shortfall", lambda tr, y: m_gbt(tr, target=y, fs=["pshort"])),
        ("child's QPB shortfall", lambda tr, y: m_gbt(tr, target=y, fs=["cqpb_s"])),
        ("structure: sums, nnz, clusters, GL, fixed interaction", lambda tr, y: m_gbt(tr, target=y, fs=["rsum", "rnnz", "csum", "cnnz", "rclu", "cclu", "gl", "dfix"])),
        ("agrees with the best known", lambda tr, y: m_gbt(tr, target=y, fs=["agree", "par_agree"])),
        ("all but the child's QPB (GBT)", lambda tr, y: m_gbt(tr, False, target=y)),
        ("all, with the child's QPB (GBT)", lambda tr, y: m_gbt(tr, True, target=y)),
    ]
    for y in ("ylog", "open", "gainrel"):
        print(f"\n## Out-of-fold R2 of {y}: total / within parent\n")
        print("| predictor | " + " | ".join(RUNS) + " |")
        print("|---|" + "---|" * len(RUNS))
        for name, fit in sets:
            cells = []
            for run, d in ds.items():
                need = {"child's QPB shortfall": "cqpb_s", "agrees with the best known": "agree"}.get(name)
                if need and not d[need].notna().any():
                    cells.append("-"); continue
                if name.startswith("all, with") and not d.cqpb_s.notna().any():
                    cells.append("-"); continue
                dd = d if len(d) < 60000 else d[d.pkey.isin(pd.Series(d.pkey.unique()).sample(frac=60000 / len(d), random_state=1))]  # tai60b: a sample of parents
                t, w = r2s(dd, y, oof(dd, fit, y))
                cells.append(f"{t:.2f} / {w:.2f}")
            print(f"| {name} | " + " | ".join(cells) + " |")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
