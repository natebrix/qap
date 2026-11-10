#!/usr/bin/env python3
"""ucmp.py TRUTHPREFIX GROUPPREFIX: each union node against its members bounded alone."""
import re, sys, glob, os
U = os.path.dirname(os.path.abspath(__file__))
CUT = 48814.05
kid = {}; par = None
for fn in sorted(glob.glob(os.path.join(U, sys.argv[1] + '*.log'))):
    name = None
    for ln in open(fn):
        if ln.startswith('CHILDNODE'): name = ln.split()[1]
        m = re.match(r'RESULT node=(\d+) m=\d+ lb=([\d.\-]+) .*status=(\w+) iters=(\d+) cpu=([\d.]+)', ln)
        if m and int(m.group(1)) == 0: par = (float(m.group(2)), m.group(3), int(m.group(4)), float(m.group(5)))
        if m and name and int(m.group(1)) > 0:
            kid[name] = (m.group(3), int(m.group(4)), float(m.group(2)), float(m.group(5))); name = None
print('parent: %.1f short after %d it (%s); children %d, fathomed %d' % (CUT - par[0], par[2], par[1], len(kid), sum(k[0] == 'fathomed' for k in kid.values())))
rows = []
for fn in sorted(glob.glob(os.path.join(U, sys.argv[2] + '*.log'))):
    mem = None; hist = {}
    for ln in open(fn):
        m = re.match(r'GROUP-PASS (\d+) size=(\d+) line=\w+ members=(\S+)', ln)
        if m: mem = m.group(3).split(','); hist = {}; continue
        mi = re.match(r'\s*it\s+(\d+) lb ([\d.\-]+)', ln)
        if mi and mem is not None: hist[int(mi.group(1))] = float(mi.group(2)); continue
        m = re.match(r'RESULT group=(\d+) m=\d+ lb=([\d.\-]+) .*status=(\w+) iters=(\d+) cpu=([\d.]+)', ln)
        if m and mem:
            rows.append((len(mem), mem, m.group(3), int(m.group(4)), float(m.group(2)), float(m.group(5)), dict(hist))); mem = None
rows.sort(key=lambda r: (r[0], r[3]))
print('size | union: status, it, CPU, end (vs cutoff) | members alone: it (sum), CPU (sum), slowest it, weakest end | union/members CPU | union at it 10/20/30 vs parent')
for n, mem, st, it, lb, cpu, hist in rows:
    ks = [kid[x] for x in mem if x in kid]
    si = sum(k[1] for k in ks); sc = sum(k[3] for k in ks)
    h = ' '.join(('%+.0f' % (hist[i] - par[0])) if i in hist else '.' for i in (10, 20, 30))
    print('%2d | %-9s %4d it %6.1f s %+7.1f | %5d it %6.1f s  slowest %3d  weakest %+6.1f | %.2f | %s | %s' % (
        n, st, it, cpu, lb - CUT, si, sc, max(k[1] for k in ks), min(k[2] for k in ks) - CUT, cpu / sc if sc else 0, h, ','.join(mem[:4]) + (',...' if len(mem) > 4 else '')))
