#!/usr/bin/env python3
"""split.py: the survive rule (ex_k <= s_par) split by where the parent is: a chain node (its own best completion
within 3,000 of the best known: it holds the best known solution or tai60b's mirror) or off the chains (inside a
subtree whose best completion is above the cutoff). Also "gap": ex_k <= s_par + ex_k(parent). Writes out/split.md."""
import numpy as np, pandas as pd, lsq
out = ['| run | parents | children | branched | rule (s_par): precision / recall | gap rule: P / R |', '|---|---|---|---|---|---|']
for run in ['tai35b_base', 'tai40b_base', 'tai50b-261116', 'tai50b-r9', 'tai60b-261117']:
    k = pd.read_pickle(f'{lsq.H}/out/kid_{run}.pkl')
    P = pd.read_csv(f'{lsq.H}/out/{run}.tsv', sep='\t', usecols=['fixes', 'ex_k'], dtype={'fixes': str}, keep_default_na=False)
    exP = dict(zip(P.fixes, pd.to_numeric(P.ex_k.replace('-', np.nan)))); k['exP'] = k.parent.map(exP)
    for lab, m in [('chain', k.exP <= 3000), ('off-chain', k.exP > 3000)]:
        g = k[m]
        if not len(g): continue
        def pr(pred):
            tp = (pred & g.surv).sum(); fp = (pred & ~g.surv).sum(); fn = (~pred & g.surv).sum()
            return f'{tp / max(tp + fp, 1):.2f} / {tp / max(tp + fn, 1):.2f}'
        out.append(f'| {run} {lab} | {g.parent.nunique()} | {len(g)} | {int(g.surv.sum())} | {pr(g.ex_k <= g.s_par)} | {pr(g.ex_k <= g.s_par + g.exP)} |')
txt = '\n'.join(out); print(txt); open(lsq.H + '/out/split.md', 'w').write(txt + '\n')
