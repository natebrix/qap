"""Subtree sizes and CPU under named nodes of a node trace (fixes in fixing order: the parent is the prefix).
usage: subtrees.py TRACE.gz KEY [KEY ...]   (KEY = fixes as in the trace, e.g. 44:58,43:45)"""
import sys, gzip, collections
tr, keys = sys.argv[1], sys.argv[2:]
cnt = collections.Counter(); cpu = collections.Counter(); line = {}
with gzip.open(tr, 'rt') as f:
    for ln in f:
        t = ln.split()
        if not t or t[0] != 'node': continue
        fx = t[2]; c = float(t[11]) if t[11] != '-' else 0.0
        parts = [] if fx == '-' else fx.split(',')
        for d in range(len(parts) + 1):
            k = ','.join(parts[:d]) if d else '-'
            cnt[k] += 1; cpu[k] += c
        line[fx] = ln.strip()
tot = cnt['-']
for k in keys:
    print(f'{k}: {cnt[k]} nodes ({100*cnt[k]/tot:.1f}% of {tot}), {cpu[k]/3600:.1f} CPU h ({100*cpu[k]/cpu["-"]:.1f}%)')
    print('   ', line.get(k, '(not in trace)'))
