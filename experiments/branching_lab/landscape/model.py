#!/usr/bin/env python3
"""model.py RUN [RUN...]: fit the size model of a child from x = ex_k / s (s = the parent's shortfall, s_par, or the
on-path shortfall at the child's depth, s_pathd; SLACK env, default par) on the production traces:
  P(branch | x) = logistic(p0 + p1 log max(x, 0.01))   (all children of branched nodes)
  log subtree(child) = a_parent + c log max(x, 0.05)   (branched off-path children; a_parent a fixed effect per
  parent: demeaned; b, the same in x rather than log x, is printed too but fits worse)
so a line's predicted cost at a node is proportional to sum over its children of P(x) max(x, 0.05)^c (the fathomed
children's cost, one node each, is left out: it is the same m nodes for every line). Writes out/model.json and
prints the fit, the subtree sizes by x bin, and how much of each parent's survivors' subtree the model explains."""
import sys, os, json, numpy as np, pandas as pd, lsq
from scipy.optimize import minimize

SL = os.environ.get('SLACK', 'par')
def sizes(run):
    d = pd.read_csv(f'{lsq.H}/out/{run}.tsv', sep='\t', usecols=['depth', 'fixes'], dtype={'fixes': str})
    size = dict.fromkeys(d.fixes, 1)
    for f, dp in sorted(zip(d.fixes, d.depth), key=lambda x: -x[1]):
        if dp >= 1:
            p = ','.join(f.split(',')[:-1]) if dp > 1 else '-'
            if p in size: size[p] += size[f]
    return size
K = []
for run in sys.argv[1:]:
    k = pd.read_pickle(f'{lsq.H}/out/kid_{run}.pkl'); sz = sizes(run)
    k['size'] = k.fixes.map(sz); k['x'] = k.ex_k / k[f's_{SL}']; k['run'] = run
    K.append(k[['run', 'depth', 'fixes', 'parent', 'surv', 'onpath', 'x', 'size', 'ex_k', 's_par', 's_pathd']])
K = pd.concat(K); K = K[np.isfinite(K.x)]
x, y = np.log(K.x.clip(lower=0.01).values), K.surv.values.astype(float)   # the logistic is in u = log max(x, 0.01)
nll = lambda w: np.sum(np.logaddexp(0, -(w[0] + w[1] * x)) * y + np.logaddexp(0, w[0] + w[1] * x) * (1 - y))
w = minimize(nll, [0.0, -1.0], method='Nelder-Mead').x
S = K[K.surv & ~K.onpath].copy()
S['ly'] = np.log(S['size']); S['lyd'] = S.ly - S.groupby(['run', 'parent']).ly.transform('mean')
S['xd'] = S.x - S.groupby(['run', 'parent']).x.transform('mean')
b = float((S.xd * S.lyd).sum() / (S.xd ** 2).sum())
S['lx'] = np.log(S.x.clip(lower=0.01)); S['lxd'] = S.lx - S.groupby(['run', 'parent']).lx.transform('mean')
c = float((S.lxd * S.lyd).sum() / (S.lxd ** 2).sum())     # log size = a_parent + c log x: the model used
print(f'within-parent log size on log x: slope c = {c:.2f}'); print(f'slack {SL}: {len(K)} children, {int(y.sum())} branched; P(branch|x) = logistic({w[0]:.2f} + {w[1]:.2f} log x); '
      f'log size of the off-path branched = a_parent + {b:.2f} x ({len(S)} children, {S.groupby(["run", "parent"]).ngroups} parents)')
bins = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10, np.inf]
K['bin'] = pd.cut(K.x, bins, right=False)
t = K.groupby('bin', observed=True).agg(children=('surv', 'size'), branched=('surv', 'sum'), mean_size=('size', 'mean'),
                                        model_p=('x', lambda v: float(np.mean(1 / (1 + np.exp(-(w[0] + w[1] * np.log(np.clip(v, 0.01, None)))))))))
print(t.to_string())
json.dump({'slack': SL, 'p0': float(w[0]), 'p1': float(w[1]), 'b': b, 'c': c, 'xmin': 0.05, 'scale': 1e9, 'runs': sys.argv[1:]},
          open(f'{lsq.H}/out/model{"" if SL == "par" else "_" + SL}.json', 'w'), indent=1)
