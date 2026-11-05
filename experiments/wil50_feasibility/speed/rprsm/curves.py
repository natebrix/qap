#!/usr/bin/env python3
"""Tabulate nug20 root curves: certified lb at fixed iterations, and the first
check reaching a few levels, for logs/nug20_<name>.log.
   usage: curves.py name [name ...]"""
import os, re, sys
D = os.path.dirname(os.path.abspath(__file__))
PTS = [10, 30, 50, 70, 90, 130, 170, 210, 250, 290, 330, 370, 400]
LEVELS = [2500.0, 2504.0, 2505.0, 2505.5, 2505.9]


def load(name):
    tr, ms = {}, None
    for line in open(os.path.join(D, "logs", "nug20_%s.log" % name)):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m: tr[int(m.group(1))] = float(m.group(2))
        m = re.search(r"ms_per_iter=(\S+)", line)
        if m: ms = float(m.group(1))
    return tr, ms


print("# it\t" + "\t".join(sys.argv[1:]))
data = {n: load(n) for n in sys.argv[1:]}
for p in PTS:
    print("%d\t" % p + "\t".join("%.4f" % data[n][0].get(p, float("nan")) for n in sys.argv[1:]))
print("# first check at or above the level")
for lv in LEVELS:
    row = []
    for n in sys.argv[1:]:
        its = [i for i, v in sorted(data[n][0].items()) if v >= lv]
        row.append(str(its[0]) if its else "-")
    print("%.1f\t" % lv + "\t".join(row))
print("# ms/iteration\t" + "\t".join("%.2f" % data[n][1] for n in sys.argv[1:]))
