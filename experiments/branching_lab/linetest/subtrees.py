#!/usr/bin/env python3
"""subtrees.py: production's (tai60b-261117) subtree under each branching child of 10:39 (nodes, CPU s), the
child's shortfall there, and a log-linear fit log(nodes) = a + b * shortfall over the 9 off-path children; then
that fit applied to the strong arms' branching children (sc54, sc47), a rough [e] size per line. Writes
out/subtrees.txt."""
import gzip, math, collections, numpy as np
NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39,'
cnt = collections.Counter(); cpu = collections.Counter()
for ln in gzip.open('/home/user/wt-cognac/infra/phase1/runs/tai60b-261117/nodes.txt.gz', 'rt'):
    if NODE not in ln: continue
    f = ln.split()
    if not f[2].startswith(NODE): continue
    k = f[2][len(NODE):].split(',')[0]; cnt[k] += 1; cpu[k] += float(f[11])
tot_n = sum(cnt.values()); tot_c = sum(cpu.values())
print(f'production, under 10:39: {tot_n:,} nodes, {tot_c/3600:,.0f} CPU h (children and below)')
rows = []
for f in (l.split() for l in open('out/score_c54.tsv').read().split('\n')[1:] if l): pass
sc = {}
for arm in ('sc54', 'sc47'):
    sc[arm] = [l.split('\t') for l in open(f'out/score_{arm}.tsv').read().split('\n')[1:] if l]
prodsf = {}
for l in open('out/score_sc54.tsv').read().split('\n')[1:]:
    if not l: continue
    t = l.split('\t'); po = t[10]
    if po.startswith('B'): prodsf[t[0]] = float(po.split()[1].replace(',', ''))
print('child   production shortfall   nodes   CPU h')
for k in sorted(prodsf, key=lambda c: -cnt[c]):
    print(f'{k:6s} {prodsf[k]:>12,.0f} {cnt[k]:>9,} {cpu[k]/3600:8.1f}')
path = '14:54'
xs = [prodsf[k] for k in prodsf if k != path]; ys = [math.log(cnt[k]) for k in prodsf if k != path]
b, a = np.polyfit(xs, ys, 1)
print(f'fit, off-path children: log nodes = {a:.2f} + {b*1e5:.3f} per 100k of shortfall (e-fold every {1/b:,.0f}); corr {np.corrcoef(xs, ys)[0,1]:.2f}')
for arm in ('sc54', 'sc47'):
    br = [(t[0], float(t[3])) for t in sc[arm] if t[1] == 'B']
    off = [(c, s) for c, s in br if c not in ('14:54', '47:47')]
    est = sum(math.exp(a + b * s) for c, s in off)
    print(f'{arm}: {len(br)} branching; off-path {len(off)}, fit sum {est:,.0f} nodes [e]; path child', [ (c, f"{s:,.0f}") for c, s in br if c in ("14:54", "47:47")])
