#!/usr/bin/env python3
"""tho.py log...: per dnntest -K log, the children's status and iterations (RESULT lines after the parent),
the histogram of iterations to fathom, their sum, the extra certificates, and the children's CPU."""
import sys, re, collections
for f in sys.argv[1:]:
    its = []; st = collections.Counter(); extra = 0; cpu = 0.0; first = True; par = ''
    for l in open(f):
        if l.startswith('admm cross:'):
            if not first: extra += int(l.split()[2])
        if l.startswith('RESULT'):
            m = re.search(r'lb=([-\d.e+]+) .*status=(\w+) iters=(\d+) cpu=([\d.]+)', l)
            if first: first = False; par = 'parent lb %.1f %s@%s' % (float(m.group(1)), m.group(2), m.group(3)); continue
            st[m.group(2)] += 1; its.append(int(m.group(3))); cpu += float(m.group(4))
    h = collections.Counter(its)
    print(f.split('/')[-1], '|', par, '| children', len(its), dict(st), '| sum its', sum(its), '| median', sorted(its)[len(its)//2] if its else '-',
          '| extra certs', extra, '| cpu %.1f' % cpu)
    print('   iterations: count ->', ' '.join('%d:%d' % (k, h[k]) for k in sorted(h)))
