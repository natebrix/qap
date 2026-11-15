#!/usr/bin/env python3
"""kids.py its log...: per child of a dnntest -K log, the certified running best at the given iterations
and how the child ended (status@iterations). The parent's final value comes first."""
import sys, re
its = [int(x) for x in sys.argv[1].split(',')]
for f in sys.argv[2:]:
    name = f.split('/')[-1].replace('.log', '')
    cur = 'parent'; d = {}; out = []
    def flush(res):
        if not d: return
        row = ['%.0f' % d[max(j for j in d if j <= i)] if any(j <= i for j in d) and i <= max(d) else '-' for i in its]
        out.append((cur, row, res))
    for l in open(f):
        m = re.match(r'it\s+(\d+) lb ([-\d.e+]+)', l)
        if m: d[int(m.group(1))] = float(m.group(2)); continue
        if l.startswith('CHILD-SPEC'): nxt = l.split()[1]
        if l.startswith('RESULT'):
            st = re.search(r'status=(\w+) iters=(\d+)', l)
            flush(st.group(1) + '@' + st.group(2)); d = {}
            cur = None
        if l.startswith('CHILDNODE'): cur = l.split()[1]
    for c, row, res in out:
        print(name, c, *row, res, sep=' | ')
