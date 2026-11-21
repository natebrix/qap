#!/usr/bin/env python3
"""lines.py ERR [DEPTH] [LINE ...]: per free line of a node's DNN marginals x-bar (the solver's QAP_LINE_SCORES
dump in its stderr: ls-rows / ls-cols / ls-X k:), the entries above 0.001, 0.0001 and 0.01, the sum of squares,
the largest entry and its cell. Sorted as the rule "tail" orders them (fewest entries above 0.001, then the larger
sum of squares, rows before columns). LINE (r23, c23, ...) marks lines to show wherever they rank; without any,
the first 12 lines and rule 9's two candidates (the row and the column with the largest sum of squares)."""
import sys
err = sys.argv[1]
depth = int(sys.argv[2]) if len(sys.argv) > 2 else 0
want = set(sys.argv[3:])
blocks, cur = [], None
for l in open(err):
    if l.startswith('ls-node:'):
        cur = {'hdr': l.strip(), 'X': {}}
    elif cur is None:
        continue
    elif l.startswith('ls-rows:'): cur['rows'] = [int(x) for x in l.split()[1:]]
    elif l.startswith('ls-cols:'): cur['cols'] = [int(x) for x in l.split()[1:]]
    elif l.startswith('ls-X '):
        k, v = l[5:].split(':', 1); cur['X'][int(k)] = [float(x) for x in v.split()]
    elif l.startswith('ls-end'):
        blocks.append(cur); cur = None
for b in blocks:
    if int(b['hdr'].split()[2]) != depth: continue
    m = len(b['X']); X = [b['X'][k] for k in range(m)]; rows, cols = b['rows'], b['cols']
    L = []
    for side in (0, 1):
        for k in range(m):
            v = X[k] if side == 0 else [X[a][k] for a in range(m)]
            j = max(range(m), key=lambda i: v[i])
            name = ('r%d' % rows[k]) if side == 0 else ('c%d' % cols[k])
            cell = ('%d:%d' % (rows[k], cols[j])) if side == 0 else ('%d:%d' % (rows[j], cols[k]))
            L.append(dict(name=name, side=side, k=k, c3=sum(x > 1e-3 for x in v), c4=sum(x > 1e-4 for x in v),
                          c2=sum(x > 1e-2 for x in v), c1=sum(x > 0.1 for x in v), ss=sum(x * x for x in v), mx=v[j], cell=cell))
    L.sort(key=lambda d: (d['c3'], -d['ss'], d['side'], d['k']))
    for i, d in enumerate(L): d['rank'] = i + 1
    r9 = {max((d for d in L if d['side'] == s), key=lambda d: d['ss'])['name'] for s in (0, 1)}
    print(b['hdr']); print('m %d, %d lines; rule 9 candidates: %s' % (m, len(L), ' '.join(sorted(r9))))
    print('%4s %-5s %7s %7s %7s %5s %7s %7s  %s' % ('rank', 'line', '>0.001', '>0.0001', '>0.01', '>0.1', 'sumsq', 'max', 'cell'))
    for d in L:
        if (not want and d['rank'] <= 12) or d['name'] in want or d['name'] in r9:
            print('%4d %-5s %7d %7d %7d %5d %7.4f %7.4f  %s' % (d['rank'], d['name'], d['c3'], d['c4'], d['c2'], d['c1'], d['ss'], d['mx'], d['cell']))
    cs = sorted(d['c3'] for d in L)
    print('entries above 0.001 over the lines: min %d, quartiles %d %d %d, max %d' % (cs[0], cs[len(cs) // 4], cs[len(cs) // 2], cs[3 * len(cs) // 4], cs[-1]))
