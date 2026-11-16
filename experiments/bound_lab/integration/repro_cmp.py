#!/usr/bin/env python3
"""repro_cmp.py MINE.dnn.txt AUTHOR.dnn.txt [--depth D] : two QAP_DNN_TRACE files node by node (keyed by depth and fixes):
DNN nodes, fathomed, first-solve iterations, cut phases, phase iterations, and how many nodes agree in
(first status, first iterations, phase iterations, phase best to 4 decimals). --depth D: only nodes of depth <= D."""
import re, sys
def load(p, dmax):
    out = {}
    for l in open(p):
        m = re.search(r'depth (\d+) fixes (.*?) \| first lb (\S+) primal \S+ status (\w+) iters (\d+)', l)
        if not m or int(m.group(1)) > dmax: continue
        c = re.search(r'cut 1 best (\S+) rounds \d+ iters (\d+) cpu \S+ status (\w+)(?: cp \d+ active \d+ stalled (\d+)(?: gaveup (\d+))?)?', l)
        out[(int(m.group(1)), m.group(2))] = dict(lb=m.group(3), st=m.group(4), it=int(m.group(5)),
            cut=(c.group(1), int(c.group(2)), c.group(3), c.group(5) or '0') if c else None, fath=' fathomed 1' in l)
    return out
dmax = int(sys.argv[sys.argv.index('--depth') + 1]) if '--depth' in sys.argv else 10**9
a, b = load(sys.argv[1], dmax), load(sys.argv[2], dmax)
def tot(x):
    return 'DNN nodes %d, fathomed %d, first-solve iterations %d, cut phases %d (gave up %d), phase iterations %d, all iterations %d' % (
        len(x), sum(v['fath'] for v in x.values()), sum(v['it'] for v in x.values()), sum(v['cut'] is not None for v in x.values()),
        sum(v['cut'] is not None and v['cut'][3] == '1' for v in x.values()),
        sum(v['cut'][1] for v in x.values() if v['cut']), sum(v['it'] + (v['cut'][1] if v['cut'] else 0) for v in x.values()))
print('here:   ' + tot(a)); print('author: ' + tot(b))
both = set(a) & set(b)
same = sum(a[k] == b[k] for k in both)
print('nodes in both %d (only here %d, only author %d); identical in every compared field %d' % (len(both), len(set(a) - set(b)), len(set(b) - set(a)), same))
for k in sorted(both):
    if a[k] != b[k]: print('  first difference: depth %d fixes %s\n    here   %s\n    author %s' % (k[0], k[1], a[k], b[k])); break
