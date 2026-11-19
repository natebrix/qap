#!/usr/bin/env python3
"""The cross-run tests of predict.py, split by whether the test child's subproblem (its set of fixes) was
also bounded in the training run: there a model can recall the node instead of predicting it.
  OMP_NUM_THREADS=1 python3 overlap.py > overlap_out.md"""
import numpy as np, pandas as pd
from common import *

def key(d):
    return [frozenset(([tuple(map(int, p.split(":"))) for p in pk.split(",")] if pk else []) + [(r, c)])
            for pk, r, c in zip(d.pkey, d.row, d.col)]

print("| train -> test | test parents | of which every child's subproblem is new | score | AUC, all | AUC, new parents | Spearman, new parents |")
print("|---|---|---|---|---|---|---|")
for a, b in [("tai50b-261116", "tai50b-r9"), ("tai50b-r9", "tai50b-261116"), ("tho40ab-0", "tho40ab-32")]:
    tr, te = load(a), load(b)
    seen = set(key(tr)); te = te.assign(seen=[k in seen for k in key(te)])
    newp = te.groupby("pkey").seen.transform("max") == False
    for name, fit in [("qpb", m_qpb), ("pc", m_pseudocost), ("gbt", lambda t: m_gbt(t, False)), ("gbt+q", lambda t: m_gbt(t, True))]:
        s = fit(tr)(te)
        ra = within_parent(te, s); rn = within_parent(te[newp], s[newp.values])
        print(f"| {a} -> {b} | {ra['parents']} | {rn['parents']} | {name} | {fmt(ra['auc'])} | {fmt(rn['auc'])} | {fmt(rn['spearman'])} |", flush=True)
