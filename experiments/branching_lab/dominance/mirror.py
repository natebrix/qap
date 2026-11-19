"""tai60b: the (44,43) exchange along the mirror chain 44:58 (and the path chain 44:45, its twin): min Delta over
the node's completions for the pair (44, 43) itself, by depth; and the mirror's best completion's Delta.
Also the pair test's time by depth (from the replay). usage: mirror.py <replay.tsv.gz>"""
import sys, gzip, collections, numpy as np
from qapdata import read_dat, read_sln, cost
from dom import pair_min_general
A, B = read_dat('tai60b'); opt, p = read_sln('tai60b')
rows = {}
with gzip.open(sys.argv[1], 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        r = dict(zip(hdr, ln.rstrip('\n').split('\t'))); rows[r['key']] = r
kids = collections.defaultdict(list)
for k in rows:
    if k != '-': kids[k.rsplit(',', 1)[0] if ',' in k else '-'].append(k)
sub = {}
for k in sorted(rows, key=lambda k: -int(rows[k]['depth'])): sub[k] = 1 + sum(sub[c] for c in kids[k])
q = p.copy(); q[[43, 44]] = q[[44, 43]]
print('mirror of the best known solution (rows 43, 44 exchanged): cost - opt =', cost(A, B, q) - opt)
for root in ['44:58', '44:45']:
    k = root; out = []
    while True:
        F = {int(a): int(b) for a, b in (s.split(':') for s in k.split(','))}
        if 43 in F and 44 in F: out.append((int(rows[k]['depth']), pair_min_general(A, B, F, 44, 43)))
        if not kids[k]: break
        k = max(kids[k], key=lambda c: sub[c])
    print(f'chain {root}: min Delta of the (44,43) exchange by depth: ' + ' '.join(f'd{d}:{v}' for d, v in out[::3] + [out[-1]]))
t = collections.defaultdict(list)
for r in rows.values(): t[int(r['depth'])].append((float(r['t_pair']), float(r['t_one']), float(r['cpu'])))
print('depth  nodes  pair ms  one-fixed ms  solver s/node (means)')
for d in sorted(t)[::5]:
    a = np.array(t[d]); print(f'{d:5d} {len(a):6d} {1e3*a[:,0].mean():8.3f} {1e3*a[:,1].mean():10.2f} {a[:,2].mean():10.2f}')
