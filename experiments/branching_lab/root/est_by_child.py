#!/usr/bin/env python3
"""est_by_child.py DIVES.jsonl.gz CHILDREN.txt ROW [DEPTH]: an estimate's rest (coord/dives.py: per open node its
exact part plus the mean of its dives, over pi) totalled under each depth-1 node, beside the root lookahead's
shortfall for the same child of ROW (CHILDREN: QAP_LA_PROGRESS lines, "child depth 0 cap C cell R:C lb .. shortfall ..").
Prints the table and the rank correlations. The estimate's "done" part (the solve's own nodes) is not in the
dives file and is left out: on tho40-est2-261119 it is 18,767 nodes of 2.07 million."""
import sys, gzip, json, statistics, collections
dv, ch, row = sys.argv[1], sys.argv[2], sys.argv[3]
strata, dives = {}, collections.defaultdict(list)
for l in gzip.open(dv, 'rt'):
    d = json.loads(l)
    if d['type'] == 'stratum': strata[d['key']] = d
    elif d['type'] == 'dive': dives[d['key']].append(d)
N, C, K = collections.Counter(), collections.Counter(), collections.Counter()
for k, s in strata.items():
    pi = float(s.get('pi', 1.0)) or 1.0
    ds = dives.get(k, [])
    n = s['exact_nodes'] + (statistics.fmean(x['nodes'] for x in ds) if ds else 0.0)
    c = s['exact_cpu'] + (statistics.fmean(x['cpu'] for x in ds) if ds else 0.0)
    top = k.split(',')[0]
    N[top] += n / pi; C[top] += c / pi; K[top] += 1
sh = {}
for l in open(ch):
    t = l.split()
    if t[:1] == ['child'] and t[2] == '0' and t[6].split(':')[0] == row: sh[t[6]] = (float(t[10]), float(t[8]))
def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]: j += 1
        for q in range(i, j + 1): r[s[q]] = (i + j) / 2 + 1
        i = j + 1
    return r
def corr(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    c = sum((x - ma) * (y - mb) for x, y in zip(a, b)); va = sum((x - ma) ** 2 for x in a); vb = sum((y - mb) ** 2 for y in b)
    return c / (va * vb) ** 0.5
tn, tc = sum(N.values()), sum(C.values())
print(f'estimate rest from the dives file: {tn:,.0f} nodes, {tc / 3600:,.1f} CPU h, {len(strata)} open nodes')
print('| depth-1 node | lookahead shortfall (150 it) | lb | open nodes sampled | estimated nodes | share | estimated CPU h | share |')
print('|---|---|---|---|---|---|---|---|')
keys = sorted(N, key=lambda k: -C[k])
for k in keys:
    s = sh.get(k, (float('nan'), float('nan')))
    print(f'| {k} | {s[0]:,.1f} | {s[1]:,.1f} | {K[k]} | {N[k]:,.0f} | {100 * N[k] / tn:.1f}% | {C[k] / 3600:,.1f} | {100 * C[k] / tc:.1f}% |')
both = [k for k in keys if k in sh]
if len(both) > 2:
    a = [sh[k][0] for k in both]
    import math
    print(f'\n{len(both)} children. Spearman(shortfall, estimated nodes) = {corr(rank(a), rank([N[k] for k in both])):+.2f}; '
          f'(shortfall, estimated CPU) = {corr(rank(a), rank([C[k] for k in both])):+.2f}; '
          f'Pearson(shortfall, ln nodes) = {corr(a, [math.log(N[k]) for k in both]):+.2f}')
