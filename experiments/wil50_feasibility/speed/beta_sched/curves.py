#!/usr/bin/env python3
"""Certified lb at given iterations, one column per log: curves.py DIR name ... [-i 10,50,...]"""
import os, re, sys

args = sys.argv[1:]
its = [10, 30, 50, 70, 90, 110, 130, 150, 170, 200, 250, 300, 350, 400]
if "-i" in args:
    k = args.index("-i"); its = [int(x) for x in args[k + 1].split(",")]; del args[k:k + 2]
d, names = args[0], args[1:]
tr = {}
for n in names:
    tr[n] = {}
    for line in open(os.path.join(d, n + ".log")):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m: tr[n][int(m.group(1))] = float(m.group(2))
print("it\t" + "\t".join(names))
for i in its:
    print("%d\t" % i + "\t".join("%.4f" % tr[n][i] if i in tr[n] else "-" for n in names))
