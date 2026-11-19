#!/usr/bin/env python3
"""chainscore.py RUN [RUN...]: chain/RUN.tsv (chain.py) at each chain node: every free line scored by the landscape
signal, s = cutoff - the node's own bound (after its cut phase: what the node knows when it branches):
  cnt   children with exK <= s (predicted survivors)
  cnt0  the same with ex0 (descents only)
  mod   sum over children of P(x) exp(b x), x = exK / s (out/model.json): a predicted cost, relative
For the line production took: predicted survivors against realised (the trace's B children), its rank among all
free lines by cnt and by mod; the line the signal would take (least mod, ties by cnt), its cnt, and the model's
ratio mod(taken) / mod(signal's line). Writes out/chain_RUNS.md."""
import sys, json, numpy as np, pandas as pd, lsq

M = json.load(open(lsq.H + '/out/model.json'))
def P(x): return 1 / (1 + np.exp(-(M['p0'] + M['p1'] * np.log(np.clip(x, 0.01, None)))))
out = []
for run in sys.argv[1:]:
    d = pd.read_csv(f'{lsq.H}/chain/{run}.tsv', sep='\t', dtype={'node': str, 'child': str, 'out': str, 'taken': str, 'line': str})
    out.append(f'## {run}\n')
    out.append('| depth | node (last fix) | s | taken | taken: predicted / realised survivors | taken: rank by cnt, by mod (of lines) | '
               'signal\'s line: cnt | mod(taken) / mod(signal) | lines with cnt <= taken\'s |')
    out.append('|---|---|---|---|---|---|---|---|---|')
    allk = []
    for node, g in d.groupby('node', sort=False):
        s = g.cutoff.iloc[0] - g.lb.iloc[0]
        g = g.assign(x=g.exK / s)
        L = g.groupby('line').agg(cnt=('exK', lambda v: int((v <= s).sum())), cnt0=('ex0', lambda v: int((v <= s).sum())),
                                  mod=('x', lambda v: float((P(v.values) * np.exp(M['b'] * v.values)).sum())),
                                  m=('x', 'size'))
        tk = g.taken.iloc[0]
        tkk = g[g.line == tk]; real = int(tkk.out.isin(['B', 'E']).sum())
        allk.append(tkk.assign(s=s))
        rc = int((L.cnt < L.loc[tk, 'cnt']).sum()) + 1; rm = int((L['mod'] < L.loc[tk, 'mod']).sum()) + 1
        best = L.sort_values(['mod', 'cnt']).index[0]
        out.append(f"| {g.depth.iloc[0]} | {node.split(',')[-1]} | {s:,.0f} | {tk} | {L.loc[tk, 'cnt']} / {real} | {rc}, {rm} (of {len(L)}) | "
                   f"{best}: {L.loc[best, 'cnt']} | {L.loc[tk, 'mod'] / L.loc[best, 'mod']:.2f} | {int((L.cnt <= L.loc[tk, 'cnt']).sum())} |")
    K = pd.concat(allk); br = K.out.isin(['B', 'E']); pr = K.exK <= K.s
    out.append(f'\nThe taken lines\' children, child by child: tp {int((pr & br).sum())}, fp {int((pr & ~br).sum())}, '
               f'fn {int((~pr & br).sum())}, of {len(K)} children ({int(br.sum())} branched)\n')
txt = '\n'.join(out); print(txt); open(f'{lsq.H}/out/chain_{"_".join(sys.argv[1:])}.md', 'w').write(txt + '\n')
