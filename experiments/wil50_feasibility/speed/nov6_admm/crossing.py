#!/usr/bin/env python3
"""From dnntest -v logs: for every fathomed solve, the fathom check F, the certified
lb at F-10 and F, the cutoff, and the linearly interpolated crossing iteration.
Prints per-solve rows and the mean/median of the iterations spent past the crossing."""
import re, sys, statistics
rows = []
for f in sys.argv[1:]:
    cutoff = None; trace = {}
    for line in open(f, errors="replace"):
        m = re.search(r"cutoff=(\S+)", line)
        if m and line.startswith("node"):
            cutoff = float(m.group(1)); trace = {}
        m = re.match(r"it\s+(\d+) lb (\S+) \(cert ([0-9.eE+-]+)", line)
        if m:
            trace[int(m.group(1))] = float(m.group(3))
        if line.startswith("RESULT") and "fathomed=1" in line or (line.startswith("RESULT") and "status=fathomed" in line):
            it = int(re.search(r"iters=(\d+)", line).group(1))
            if cutoff is None or it not in trace or (it - 10) not in trace:
                continue
            a, b = trace[it - 10], trace[it]
            if b <= cutoff: continue
            frac = (cutoff - a) / (b - a) if b > a else 1.0
            cross = it - 10 + 10 * frac
            waste = it - cross
            rows.append((f.split("/")[-1], it, a, b, cutoff, cross, waste))
for r in rows:
    print("%-28s F %4d lb(F-10) %.2f lb(F) %.2f cutoff %.2f cross %.1f waste %.1f" % r)
w = [r[6] for r in rows]
if w:
    print("n %d  mean waste %.2f  median %.2f  mean F %.1f  waste/F %.3f" % (
        len(w), statistics.mean(w), statistics.median(w), statistics.mean(r[1] for r in rows),
        sum(w) / sum(r[1] for r in rows)))
