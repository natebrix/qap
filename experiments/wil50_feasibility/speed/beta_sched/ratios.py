#!/usr/bin/env python3
"""Residual ratios from the bsched lines of DNN_BETA_SCHED runs (any setting
prints them at every check with -v): ratios.py lo hi log ...
Prints, over the checks with lo <= it <= hi, the median of pres/(beta dY)
and of rhat/shat (the controller's normalized ratio), and the beta used."""
import re, statistics, sys

lo, hi = int(sys.argv[1]), int(sys.argv[2])
allr, alln = [], []
for f in sys.argv[3:]:
    r, n, b = [], [], set()
    for line in open(f):
        m = re.match(r"bsched it (\d+) beta (\S+) pres/\(beta dY\) (\S+) rhat \S+ shat \S+ rhat/shat (\S+)", line)
        if m and lo <= int(m.group(1)) <= hi:
            b.add(float(m.group(2))); r.append(float(m.group(3))); n.append(float(m.group(4)))
    if r:
        print("%-40s beta %-18s pres/(beta dY) med %.4f  rhat/shat med %.3f [%.3f, %.3f]" % (
            f, ",".join("%.4g" % x for x in sorted(b)), statistics.median(r), statistics.median(n), min(n), max(n)))
        allr += r; alln += n
if alln:
    print("all: pres/(beta dY) med %.4f  rhat/shat med %.3f (%d checks)" % (
        statistics.median(allr), statistics.median(alln), len(alln)))
