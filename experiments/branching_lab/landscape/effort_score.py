#!/usr/bin/env python3
"""effort_score.py RUN: the survive rule on effort.py's sample (all B/C children and a random sample of D children;
D children are reweighted to their full count for precision), excess by search effort: ex_k (descents from p0
and the parent's completion), ex_ilsK (an ILS of K kicks from the p0-adjusted descent), ex_rilsK (the same from
a random start: no knowledge of p0). Slack: pathd and par; threshold ex <= slack. Writes out/effort_RUN.md."""
import sys, numpy as np, pandas as pd, lsq

run = sys.argv[1]; TAG = sys.argv[2] if len(sys.argv) > 2 else ''
k = pd.read_pickle(f'{lsq.H}/out/kid_{run}.pkl')
e = pd.read_csv(f'{lsq.H}/out/effort{TAG}_{run}.tsv', sep='\t', dtype={'fixes': str})
d = k.merge(e, on='fixes')
nD = int((k.out == 'D').sum()); sD = int((d.out == 'D').sum()); wD = nD / max(sD, 1)
nBC = int(k.out.isin(['B', 'C', 'E']).sum()); sBC = int(d.out.isin(['B', 'C', 'E']).sum())
d['w'] = np.where(d.out == 'D', wD, nBC / sBC)
cols = ['ex_k'] + [c for c in e.columns if c.startswith('ex_')]
out = [f'## {run}: effort sample {len(d)} children ({sD} of {nD} D children, weight {wD:.1f}); survivors {int(d.surv.sum())}\n',
       '| excess | median excess of the survivors / pathd | slack pathd: precision / recall | slack par: P / R | rows where it is lower than ex_k |',
       '|---|---|---|---|---|']
for c in cols:
    cells = []
    for s in ['s_pathd', 's_par']:
        ok = d[s].notna(); pred = (d[c] <= d[s]) & ok; tr = d.surv & ok
        tp = (d.w * (pred & tr)).sum(); fp = (d.w * (pred & ~tr)).sum(); fn = (d.w * (~pred & tr)).sum()
        cells.append(f'{tp / (tp + fp):.2f} / {tp / (tp + fn):.2f}')
    lower = int((d[c] < d.ex_k - 0.5).sum())
    med = (d[c] / d.s_pathd)[d.surv & ~d.onpath].median()
    out.append(f'| {c} | {med:.2f} | ' + ' | '.join(cells) + f' | {lower} |')
txt = '\n'.join(out); print(txt); open(f'{lsq.H}/out/effort{TAG}_{run}.md', 'w').write(txt + '\n')
