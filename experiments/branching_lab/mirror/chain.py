"""The heaviest chain under a node of a trace (at each step the child with the largest subtree), with each chain
node's subtree size, CPU and shortfall (cutoff - its best DNN/cut value). usage: chain.py TRACE.gz ROOTKEY > out"""
import sys, gzip, collections
tr, root = sys.argv[1], sys.argv[2]
cnt = collections.Counter(); cpu = collections.Counter(); kids = collections.defaultdict(set); info = {}
with gzip.open(tr, 'rt') as f:
    for ln in f:
        t = ln.split()
        if not t or t[0] != 'node' or t[2] == '-': continue
        parts = t[2].split(','); c = float(t[11]) if t[11] != '-' else 0.0
        for d in range(1, len(parts) + 1):
            k = ','.join(parts[:d]); cnt[k] += 1; cpu[k] += c
            if d > 1: kids[','.join(parts[:d - 1])].add(k)
        vals = [float(x) for x in t[4:6] if x != '-']
        info[t[2]] = (t[7], float(t[6]) - max(vals) if vals else float('nan'))
k = root
while True:
    print(len(k.split(',')), k, cnt[k], f'{cpu[k]/3600:.1f}', info.get(k, ('?', float('nan')))[0], f'{info.get(k, ("?", float("nan")))[1]:.0f}')
    if not kids[k]: break
    k = max(kids[k], key=lambda c: cnt[c])
