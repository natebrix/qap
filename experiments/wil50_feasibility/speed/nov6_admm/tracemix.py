#!/usr/bin/env python3
"""Node mix of a pipeline_ab trace.txt: first-solve iterations by outcome and warm status,
phase iterations, and the fathom-iteration distribution of first-solve fathoms."""
import re, sys, statistics
tot = {"first_fathom": [0, 0], "first_miss": [0, 0], "phase": [0, 0]}
F = []; miss_it = []; miss_short = []; warm_miss = 0
for line in open(sys.argv[1]):
    m = re.search(r"first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+) warm (\d) pool (\d)", line)
    if not m: continue
    lb, st, it, cpu, warm, pool = float(m.group(1)), m.group(2), int(m.group(3)), float(m.group(4)), int(m.group(5)), int(m.group(6))
    cut = float(re.search(r"cutoff (\S+)", line).group(1))
    if st == "fathomed":
        tot["first_fathom"][0] += it; tot["first_fathom"][1] += 1; F.append(it)
    else:
        tot["first_miss"][0] += it; tot["first_miss"][1] += 1; miss_it.append(it); miss_short.append(cut - lb)
        warm_miss += warm
    mc = re.search(r"cut 1 best \S+ rounds (\d+) iters (\d+)", line)
    if mc:
        tot["phase"][0] += int(mc.group(2)); tot["phase"][1] += 1
allit = sum(v[0] for v in tot.values())
for k, v in tot.items():
    print("%-13s nodes %4d iterations %6d (%.1f%%)" % (k, v[1], v[0], 100.0 * v[0] / allit))
F.sort()
print("first-solve fathoms: n %d mean %.1f median %d quartiles %d %d min %d max %d" % (
    len(F), statistics.mean(F), statistics.median(F), F[len(F) // 4], F[3 * len(F) // 4], F[0], F[-1]))
print("  fathoms at <= 50: %d, 51-100: %d, 101-200: %d, > 200: %d" % (
    sum(1 for x in F if x <= 50), sum(1 for x in F if 50 < x <= 100), sum(1 for x in F if 100 < x <= 200), sum(1 for x in F if x > 200)))
if miss_it:
    print("first-solve misses: n %d (warm %d) iterations %s short at give-up %s" % (
        len(miss_it), warm_miss, sorted(miss_it), [round(x) for x in sorted(miss_short)]))
