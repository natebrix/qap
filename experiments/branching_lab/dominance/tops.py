"""tai60b: where the pair test's savings are. For each highest flagged node that the solver branched on: its chain
(root child), depth, subtree, the pair, its min Delta, and the node's twin (the fixes with the pair's columns exchanged):
is the twin in the trace, with what outcome and subtree.  usage: tops.py <replay.tsv.gz>"""
import sys, gzip, collections
rows = {}
with gzip.open(sys.argv[1], 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        r = dict(zip(hdr, ln.rstrip('\n').split('\t'))); rows[r['key']] = r
kids = collections.defaultdict(list)
for k in rows:
    if k != '-': kids[k.rsplit(',', 1)[0] if ',' in k else '-'].append(k)
sub = {}; cpu = {}
for k in sorted(rows, key=lambda k: -int(rows[k]['depth'])):
    sub[k] = 1 + sum(sub[c] for c in kids[k]); cpu[k] = float(rows[k]['cpu']) + sum(cpu[c] for c in kids[k])
fl = lambda k: rows[k]['pairmin'] != '' and int(rows[k]['pairmin']) >= 1
tops = []; st = ['-']
while st:
    k = st.pop()
    if fl(k): tops.append(k); continue
    st.extend(kids[k])
agg = collections.defaultdict(lambda: [0, 0, 0.0])
for k in tops:
    fx = k.split(','); key = (fx[0], fx[1] if len(fx) > 1 else '', rows[k]['pair'])
    a = agg[key]; a[0] += 1; a[1] += sub[k]; a[2] += cpu[k]
print('saved by (root child, depth-2 fix, flagged pair):  tops  nodes  CPU h')
for key, a in sorted(agg.items(), key=lambda x: -x[1][1])[:15]:
    print(f'  {key}  {a[0]:6d} {a[1]:7d} {a[2]/3600:8.1f}')
print('\nbranched highest flagged nodes with subtree >= 300: depth, subtree, CPU h, pair, minDelta, twin outcome/subtree')
for k in sorted([k for k in tops if rows[k]['outcome'] == 'B'], key=lambda k: -sub[k]):
    if sub[k] < 300: break
    i, kk = map(int, rows[k]['pair'].split(','))
    F = [tuple(map(int, s.split(':'))) for s in k.split(',')]
    d = dict(F); a, b = d[i], d[kk]
    tw = ','.join(f'{r}:{(b if r == i else a if r == kk else c)}' for r, c in F)
    # the twin's highest ancestor present in the trace
    t = tw
    while t not in rows and ',' in t: t = t.rsplit(',', 1)[0]
    tinfo = f"{rows[t]['outcome']} at depth {rows[t]['depth']} (sub {sub[t]})" if t in rows else 'absent'
    print(f"  d{rows[k]['depth']} sub {sub[k]:6d} {cpu[k]/3600:6.1f} h pair ({i},{kk}) min {rows[k]['pairmin']:>8} | {k[:60]}... | twin's deepest trace ancestor: {tinfo}{' (the twin itself)' if t == tw else ''}")
