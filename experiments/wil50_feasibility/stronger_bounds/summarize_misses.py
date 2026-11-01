#!/usr/bin/env python3
"""Table of the step-2 regional-cut runs on the 14 misses (misses/polyak_<id>.log):
uncut converged bound (round 0), its shortfall against the cutoff, the best
certified bound with cuts, the gain, the round it was reached in, whether the
node fathomed, and the CPU of the whole run.

  python3 summarize_misses.py [misses] [cutoff]
"""
import glob, os, re, sys

d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "misses")
cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 48814.05
rows = []
for f in sorted(glob.glob(os.path.join(d, "polyak_*.log")), key=lambda s: int(re.search(r"_(\d+)\.log", s).group(1))):
    k = int(re.search(r"_(\d+)\.log", f).group(1))
    rounds = []
    for line in open(f):
        m = re.match(r"round\s+(\d+)\s+lb (\S+) .* status (\S+)\s+iters (\d+)\s+cpu (\S+)", line)
        if m:
            rounds.append((int(m.group(1)), float(m.group(2)), m.group(3), int(m.group(4)), float(m.group(5))))
    if not rounds:
        continue
    first = rounds[0][1]
    best_r, best = max(((r, lb) for r, lb, *_ in rounds), key=lambda t: t[1])
    fath = any(st == "fathomed" for _, _, st, _, _ in rounds)
    rows.append((k, first, cutoff - first, best, best - first, best_r, fath, sum(c for *_, c in rounds), sum(i for *_, i, _ in rounds)))

print("node  uncut(1500)  short   best_lb    gain  round  fathomed  rounds  iterations   cpu_s")
for k, first, short, best, gain, r, fath, cpu, its in rows:
    print("%4d  %11.1f  %5.1f  %9.1f  %+6.1f  %5d  %-8s  %6d  %10d  %6.0f" % (
        k, first, short, best, gain, r, "yes" if fath else "no", r if fath else 8, its, cpu))
near = [x for x in rows if x[2] < 130]
clear = [x for x in rows if x[2] >= 130]
g = sorted(x[4] for x in clear)
print("\nnear misses (< 130 short): %d, fathomed %d" % (len(near), sum(1 for x in near if x[6])))
print("clear misses: %d, gain median %.1f, range %.1f to %.1f, fathomed %d" % (
    len(clear), g[len(g) // 2], g[0], g[-1], sum(1 for x in clear if x[6])))
print("CPU per node: fathomed mean %.0f s, not fathomed mean %.0f s" % (
    sum(x[7] for x in rows if x[6]) / max(1, sum(1 for x in rows if x[6])),
    sum(x[7] for x in rows if not x[6]) / max(1, sum(1 for x in rows if not x[6]))))
