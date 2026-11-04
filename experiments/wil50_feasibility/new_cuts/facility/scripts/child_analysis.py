#!/usr/bin/env python3
"""Certificate child bounds (dnntest -C 2: 'CHILD i j lb', original indices)
against the node's bound, overall and on chosen lines; with --lines-parent,
also against the lines study's 300-iteration child bounds.
  python3 child_analysis.py LOG [--cutoff 48814.05] [--row 49 --col 45]"""
import argparse, re, json, glob, os
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('log'); ap.add_argument('--cutoff', type=float, default=48814.05)
ap.add_argument('--row', type=int, action='append', default=[]); ap.add_argument('--col', type=int, action='append', default=[])
a = ap.parse_args()
txt = open(a.log).read()
lb = float(re.search(r'RESULT .*?lb=(\S+)', txt).group(1))
ch = {(int(i), int(j)): float(v) for i, j, v in re.findall(r'^CHILD (\d+) (\d+) (\S+)$', txt, re.M)}
v = np.array(list(ch.values())); ex = v - lb
print('node lb %.2f; %d certificate child bounds: fathomed %d; excess over lb: median %.3f, 90%% %.3f, 99%% %.3f, max %.3f; within 1 of lb: %d' % (
    lb, len(v), (v > a.cutoff).sum(), np.median(ex), np.quantile(ex, .9), np.quantile(ex, .99), ex.max(), (ex < 1).sum()))
for r in a.row:
    e = sorted(ch[(r, j)] - lb for (i, j) in ch if i == r)
    print('  row %d: %d children, excess min %.2f median %.2f max %.2f' % (r, len(e), e[0], e[len(e) // 2], e[-1]))
for c in a.col:
    e = sorted(ch[(i, c)] - lb for (i, j) in ch if j == c)
    print('  col %d: %d children, excess min %.2f median %.2f max %.2f' % (c, len(e), e[0], e[len(e) // 2], e[-1]))
