#!/usr/bin/env python3
"""Read the truth logs (every child of one line bounded alone) and print the children
by speed, plus four -G strings of union groups for the second stage."""
import re, sys, glob, os
U = os.path.dirname(os.path.abspath(__file__))
CUT = 48814.05
kids = []
for fn in sorted(glob.glob(os.path.join(U, sys.argv[1] + '*.log'))):
    name = None
    for ln in open(fn):
        if ln.startswith('CHILDNODE'): name = ln.split()[1]
        m = re.match(r'RESULT node=(\d+) m=\d+ lb=([\d.\-]+) .*status=(\w+) iters=(\d+) cpu=([\d.]+)', ln)
        if m and name and int(m.group(1)) > 0:
            kids.append((name, m.group(3), int(m.group(4)), float(m.group(2)), float(m.group(5)))); name = None
kids.sort(key=lambda k: (k[1] != 'fathomed', k[2], -k[3]))
f = [k for k in kids if k[1] == 'fathomed']
print('children %d, fathomed %d; iterations (fathomed) %d, CPU %.0f s; misses %d, %d it' % (len(kids), len(f), sum(k[2] for k in f), sum(k[4] for k in f), len(kids) - len(f), sum(k[2] for k in kids if k[1] != 'fathomed')))
for k in kids: print('  %-7s %-9s %4d it  %+7.1f  %5.1f s' % (k[0], k[1], k[2], k[3] - CUT, k[4]))
n = [k[0] for k in f]
def grp(a, b): return ','.join(n[a:b])
G = []
if len(n) >= 2: G += [grp(0, 2), grp(2, 4), grp(len(n)//2, len(n)//2 + 2), grp(len(n) - 2, len(n))]
if len(n) >= 4: G += [grp(0, 4), grp(len(n)//2, len(n)//2 + 4), ','.join([n[0], n[1], n[-2], n[-1]])]
if len(n) >= 8: G += [grp(0, 8), grp(len(n) - 8, len(n))]
if len(n) >= 16: G += [grp(0, 16)]
if len(n) >= 3: G += [grp(0, len(n))]
G = [g for g in G if g.count(',') >= 1]
order = sorted(range(len(G)), key=lambda i: -G[i].count(','))
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 4
bins = [[] for _ in range(NB)]
for j, i in enumerate(order): bins[j % NB].append(G[i])
for b in bins: print('G ' + ';'.join(b))
