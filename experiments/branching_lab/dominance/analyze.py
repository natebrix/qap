"""Summarise a replay: what the pair / one-fixed tests would have saved.
usage: analyze.py <label> <instance> <replay.tsv.gz> [chain-root-key ...]
A node is flagged by a test if its value >= 1. Savings: the nodes and CPU of the subtrees under the highest
flagged node on each branch (the flagged node itself included: its own bound would not be computed)."""
import sys, gzip, collections, numpy as np, itertools
from qapdata import read_dat, read_sln, cost
label, name, path = sys.argv[1:4]
chains = sys.argv[4:]
rows = []
with gzip.open(path, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        rows.append(dict(zip(hdr, ln.rstrip('\n').split('\t'))))
node = {r['key']: r for r in rows}
kids = collections.defaultdict(list)
for r in rows:
    k = r['key']
    if k == '-': continue
    par = k.rsplit(',', 1)[0] if ',' in k else '-'
    kids[par].append(k)
# subtree totals (keys sorted by depth descending)
order = sorted(node, key=lambda k: -int(node[k]['depth']))
sub_n, sub_c = {}, {}
for k in order:
    sub_n[k] = 1 + sum(sub_n[c] for c in kids[k]); sub_c[k] = float(node[k]['cpu']) + sum(sub_c[c] for c in kids[k])
roots = [k for k in node if k == '-' or (k.rsplit(',', 1)[0] if ',' in k else '-') not in node]
TN = sum(sub_n[k] for k in roots); TC = sum(sub_c[k] for k in roots)
val = lambda r, col: (int(r[col]) if r[col] != '' else None)
def flagged(r, test):
    v1 = val(r, 'pairmin'); v2 = val(r, 'onefix')
    if test == 'pair': return v1 is not None and v1 >= 1
    if test == 'one': return v2 is not None and v2 >= 1
    return (v1 is not None and v1 >= 1) or (v2 is not None and v2 >= 1)
print(f'## {label} ({name}): {TN} nodes, {TC/3600:.1f} CPU h in the trace')
# best known solution path (and, for symmetric instances, the images of it that appear)
from autos import optimal_orbit
opt, p0 = read_sln(name); A, B = read_dat(name)
orbit, _, _ = optimal_orbit(A, B, p0)
onpath = [k for k in node if k == '-' or any(all(p[int(a)] == int(b) for a, b in (s.split(':') for s in k.split(','))) for p in orbit)]
print(f'known optima (the QAPSOLN solution and its images under the matrices automorphisms): {len(orbit)}; trace nodes containing one: {len(onpath)}')
for test in ['pair', 'one', 'either']:
    top = []  # highest flagged nodes
    stack = list(roots)
    while stack:
        k = stack.pop()
        if flagged(node[k], test): top.append(k); continue
        stack.extend(kids[k])
    sn = sum(sub_n[k] for k in top); sc = sum(sub_c[k] for k in top)
    allflag = [k for k in node if flagged(node[k], test)]
    br = [k for k in allflag if node[k]['outcome'] == 'B']
    nB = sum(1 for k in node if node[k]['outcome'] == 'B')
    print(f'{test:6s}: flagged nodes {len(allflag)}, highest flagged {len(top)}; saved {sn} nodes ({100*sn/TN:.3f}%), '
          f'{sc/3600:.2f} CPU h ({100*sc/TC:.3f}%); branched nodes flagged {len(br)} of {nB}; '
          f'best-known path nodes flagged {sum(flagged(node[k], test) for k in onpath)} of {len(onpath)}')
    byd = collections.defaultdict(lambda: [0, 0, 0.0])
    for k in top:
        d = int(node[k]['depth']); byd[d][0] += 1; byd[d][1] += sub_n[k]; byd[d][2] += sub_c[k]
    if top:
        print('   by depth of the highest flagged node: ' + '; '.join(f'd{d}: {v[0]} tops, {v[1]} nodes, {v[2]/3600:.2f} h' for d, v in sorted(byd.items())))
        oc = collections.Counter(node[k]['outcome'] for k in top)
        print('   outcomes of the highest flagged nodes:', dict(oc))
# distribution of the best pair value among branched nodes, by depth
print('branched nodes: best pair min Delta / best one-fixed value (median, max) by depth')
byd = collections.defaultdict(list)
for k, r in node.items():
    if r['outcome'] == 'B': byd[int(r['depth'])].append((val(r, 'pairmin'), val(r, 'onefix')))
for d in sorted(byd):
    pv = [a for a, b in byd[d] if a is not None]; ov = [b for a, b in byd[d] if b is not None]
    print(f'  d{d:2d} n={len(byd[d]):5d}  pair med {np.median(pv) if pv else float("nan"):14.0f} max {max(pv) if pv else float("nan"):14.0f}   '
          f'one med {np.median(ov) if ov else float("nan"):14.0f} max {max(ov) if ov else float("nan"):14.0f}')
tp = np.array([float(r['t_pair']) for r in rows]); to = np.array([float(r['t_one']) for r in rows])
cp = np.array([float(r['cpu']) for r in rows])
print(f'test time per node: pair mean {1e3*tp.mean():.3f} ms max {1e3*tp.max():.2f} ms; one-fixed mean {1e3*to.mean():.2f} ms max {1e3*to.max():.1f} ms;'
      f' solver CPU per node mean {cp.mean():.2f} s median {np.median(cp):.2f} s')
for ch in chains:
    print(f'chain under {ch}: the child with the largest subtree at each depth')
    k = ch
    while True:
        r = node[k]
        print(f"  d{r['depth']:>2} {k.split(',')[-1]:>6} {r['outcome']} sub {sub_n[k]:7d} pair {r['pairmin']:>12} ({r['pair']}) one {r['onefix']:>12} ({r['onefix_ki']})")
        if not kids[k]: break
        k = max(kids[k], key=lambda c: sub_n[c])
