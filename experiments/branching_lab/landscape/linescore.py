#!/usr/bin/env python3
"""linescore.py: the landscape signal against the line truth sets (truth/tho30.tsv, tho40.tsv, tai60b.tsv from
truth.py), in the form of ../../../../wt-cognac/experiments/line_signals/signals_out.md: per node, Spearman of the
signal with the truth (positive: agrees) and the truth of the signal's first line over the best line's (ties
in the signal: the mean over the tied lines). Child by child: precision / recall of "branches iff ex <= slack".
Slack s at a node: cutoff - the node's DNN bound (tho: a cold dnntest solve, give-up 300, truth/*_nodes_dnntest.txt;
tai60b 43:58: the line study's arms, all from one parent state, 1.76M).
Signals (smaller = better line):
  cnt1, cnt2  children with exK <= a x s (a = 1, 2), the predicted survivors; cnt1_d0 the same with ex0 (descents only)
  cntg, cntg2  children with exK - ex_node <= a x s: the excess over the node's own best completion (the min over
          its children), i.e. the slack taken as the node's gap s + ex_node
  soft    sum over children of exp(-exK / s)
  model   sum over children of the fitted size model of model.py (out/model.json), if present
  (a)     reference: the children's cold first DNN solve, shortfall sum (truth-side, costly; tho only)
Writes out/linescore.md."""
import json, os, re, numpy as np, pandas as pd, lsq
from scipy.stats import spearmanr

H = lsq.H
def node_lb(f):
    d = {}
    for ln in open(f):
        m = re.match(r'(\S+) RESULT .* lb=(\S+)', ln)
        if m: d[m[1]] = float(m[2])
    return d
LB = {**node_lb(H + '/truth/tho30_nodes_dnntest.txt'), **node_lb(H + '/truth/tho40_nodes_dnntest.txt')}
CUT = {'tho30': 149936 - 1.95, 'tho40': 240516 - 1.95, 'tai60b': 608215055 - 0.95}
TAI_TRUTH = {'r1': 19600, 'r45': 22900, 'c26': 38800, 'r8': 67100}       # ../tai60b_lines rough subtree
TAI_S = 1.76e6
MODEL = json.load(open(H + '/out/model.json')) if os.path.exists(H + '/out/model.json') else None

def model_cost(x):
    """expected subtree of a child with x = excess / slack, relative (model.py)"""
    if MODEL is None: return np.nan
    p = 1 / (1 + np.exp(-(MODEL['p0'] + MODEL['p1'] * np.log(max(x, 0.01)))))
    return p * np.exp(MODEL['b'] * x) + (1 - p) * 1.0 / MODEL['scale']

def md(R):
    f = lambda v: (f'{v:.3g}' if isinstance(v, float) else str(v))
    return '\n'.join(['| ' + ' | '.join(R.columns) + ' |', '|' + '---|' * len(R.columns)] +
                     ['| ' + ' | '.join(f(v) for v in r) + ' |' for r in R.itertuples(index=False)])

def top_over_best(sig, truth):
    sig, truth = np.asarray(sig, float), np.asarray(truth, float)
    t = truth[np.isclose(sig, sig.min())].mean()
    return t / truth.min()

SIGS = ['cnt1', 'cnt1_d0', 'cnt2', 'cntg', 'cntg2', 'soft', 'model', 'a']
out = ['# The landscape signal against the line truth sets\n']
summ = []
for SET in ['tho30', 'tho40', 'tai60b']:
    f = f'{H}/truth/{SET}.tsv'
    if not os.path.exists(f): continue
    d = pd.read_csv(f, sep='\t', dtype={'node': str, 'line': str, 'child': str})
    exK = [c for c in d.columns if c.startswith('ex')][-1]
    for node, g in d.groupby('node', sort=False):
        s = CUT[SET] - (LB[node] if SET != 'tai60b' else (CUT[SET] - TAI_S))
        g = g.copy(); exn = g[exK].min()      # the node's own best completion: the best over any line's children
        g['x'] = g[exK] / s; g['xg'] = (g[exK] - exn) / s
        # child by child
        br = g.branched.astype(int) == 1
        rows = []
        for line, L in g.groupby('line', sort=False):
            truth = (1 + L.truth.sum()) if SET != 'tai60b' else TAI_TRUTH[line]
            sh = (np.clip(CUT[SET] - L.lb.astype(float), 0, None)).sum() if SET != 'tai60b' else np.nan
            rows.append(dict(line=line, truth=truth, branched=int(L.branched.astype(int).sum()),
                             cnt1=int((L[exK] <= s).sum()), cnt2=int((L[exK] <= 2 * s).sum()),
                             cnt1_d0=int((L.ex0 <= s).sum()), cntg=int((L.xg <= 1).sum()), cntg2=int((L.xg <= 2).sum()), soft=float(np.exp(-L.x).sum()),
                             model=float(sum(model_cost(x) for x in L.x)), a=sh))
        R = pd.DataFrame(rows).sort_values('truth')
        out.append(f'## {SET} {node}: slack s = {s:,.0f}, the node\'s own excess {exn:,.0f}, {len(g)} children in {len(R)} lines\n')
        tp = int(((g[exK] <= s) & br).sum()); fp = int(((g[exK] <= s) & ~br).sum()); fn = int(((g[exK] > s) & br).sum())
        out.append(f'Child by child, branches iff {exK} <= s: tp {tp}, fp {fp}, fn {fn} (precision {tp / max(tp + fp, 1):.2f}, '
                   f'recall {tp / max(tp + fn, 1):.2f}); children that branch {int(br.sum())} of {len(g)}; '
                   f'median x of the branched {g.x[br].median():.2f}, of the rest {g.x[~br].median():.2f}; with the node\'s excess '
                   f'subtracted (xg = (ex - ex_node) / s <= 1): tp {int(((g.xg <= 1) & br).sum())}, fp {int(((g.xg <= 1) & ~br).sum())}, '
                   f'fn {int(((g.xg > 1) & br).sum())}\n')
        out.append(md(R))
        line = {'set': SET, 'node': node}
        for sig in SIGS:
            if R[sig].isna().all(): continue
            rho = spearmanr(R[sig], R.truth).correlation
            line[sig] = (rho, top_over_best(R[sig], R.truth))
        summ.append(line)
        out.append('\nSpearman with the truth / top line over best: ' + ', '.join(
            f'{k} {v[0]:+.2f} / {v[1]:.2f}x' for k, v in line.items() if isinstance(v, tuple)) + '\n')
out.append('## Summary\n')
sigs = SIGS
out.append('| node | ' + ' | '.join(sigs) + ' |'); out.append('|---|' + '---|' * len(sigs))
for l in summ:
    out.append(f"| {l['set']} {l['node']} | " + ' | '.join(
        (f'{l[s][0]:+.2f} / {l[s][1]:.2f}x' if s in l and l[s][0] == l[s][0] else (f'nan / {l[s][1]:.2f}x' if s in l else '-')) for s in sigs) + ' |')
for grp, sets in [('tho', ('tho30', 'tho40')), ('tai', ('tai60b',))]:
    L = [l for l in summ if l['set'] in sets]
    if not L: continue
    cells = []
    for s in sigs:
        v = [l[s] for l in L if s in l]
        cells.append('-' if not v else f'{np.nanmean([x[0] for x in v]):+.2f} / {np.mean([x[1] for x in v]):.2f}x (worst {max(x[1] for x in v):.2f}x)')
    out.append(f'| mean, {grp} | ' + ' | '.join(cells) + ' |')
txt = '\n'.join(out); print(txt); open(H + '/out/linescore.md', 'w').write(txt + '\n')
