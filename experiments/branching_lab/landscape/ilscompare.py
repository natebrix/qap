#!/usr/bin/env python3
"""ilscompare.py: the signal at tai60b's chain nodes with descents started from the best known solution and its
mirror (chain/tai60b-261117.tsv) against descents started only from truth/tai60b_ils.perm, a solution 13,483
above the best known that an ILS from random starts found in 25 s without knowing it (noknow.py); 30 kicks
either way. Per node: Spearman of the two line counts, whether the least-count line is the same, and the taken
line's children child by child. Writes out/ilscompare.md."""
import numpy as np, pandas as pd, lsq
from scipy.stats import spearmanr
rd = lambda f: pd.read_csv(f, sep='\t', dtype={'node': str, 'line': str, 'taken': str, 'out': str, 'child': str})
a = rd(lsq.H + '/chain/tai60b-261117.tsv'); b = rd(lsq.H + '/chain/tai60b-261117_ils.tsv')
b = b[b.node.isin(set(a.node))]
out = ['| depth | node | Spearman of the line counts | least-count line: p0 start / ILS start | taken line: predicted (p0) / (ILS) / realised |',
       '|---|---|---|---|---|']
tp = fp = fn = 0
for node, g in a.groupby('node', sort=False):
    h = b[b.node == node]
    if not len(h): continue
    s = g.cutoff.iloc[0] - g.lb.iloc[0]
    ca = g.groupby('line').exK.apply(lambda v: int((v <= s).sum())); cb = h.groupby('line').exK.apply(lambda v: int((v <= s).sum()))
    tk = g.taken.iloc[0]; real = int(g[g.line == tk].out.isin(['B', 'E']).sum())
    hb = h[h.line == tk].merge(g[g.line == tk][['child', 'out']], on='child', suffixes=('_x', ''))
    br = hb.out.isin(['B', 'E']); pr = hb.exK <= s
    tp += int((pr & br).sum()); fp += int((pr & ~br).sum()); fn += int((~pr & br).sum())
    out.append(f'| {g.depth.iloc[0]} | {node.split(",")[-1]} | {spearmanr(ca, cb.reindex(ca.index)).correlation:+.2f} | '
               f'{ca.idxmin()} ({ca.min()}) / {cb.idxmin()} ({cb.min()}) | {ca[tk]} / {cb[tk]} / {real} |')
out.append(f'\nThe taken lines\' children with the ILS start: tp {tp}, fp {fp}, fn {fn}\n')
txt = '\n'.join(out); print(txt); open(lsq.H + '/out/ilscompare.md', 'w').write(txt + '\n')
