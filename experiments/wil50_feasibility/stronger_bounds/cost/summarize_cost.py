#!/usr/bin/env python3
"""Table of the cost experiment (run_cost.sh): per protocol and node, whether
the node fathomed, the round, total iterations and CPU, the best certified
bound; then per protocol the fathom count and the mean iterations and CPU of
the fathoms, next to step 2's 700-iteration protocol.

  python3 summarize_cost.py [dir]
"""
import glob, os, re, sys

d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
step2 = {6: (5000, 3011), 9: (6400, 3736), 15: (5000, 2850), 16: (6400, 3681)}
runs = {}
for f in sorted(glob.glob(os.path.join(d, "*_*.log"))):
    m = re.match(r"(.+)_(\d+)\.log$", os.path.basename(f))
    if not m:
        continue
    proto, k = m.group(1), int(m.group(2))
    rounds = []
    for line in open(f):
        r = re.match(r"round\s+(\d+)\s+lb (\S+) .* status (\S+)\s+iters (\d+)\s+cpu (\S+)", line)
        if r:
            rounds.append((int(r.group(1)), float(r.group(2)), r.group(3), int(r.group(4)), float(r.group(5))))
    if rounds:
        runs[(proto, k)] = rounds
print("proto     node  fathomed  round  iterations   cpu_s   best_lb    (step 2: iterations, cpu)")
protos = sorted(set(p for p, _ in runs))
for p in protos:
    for k in sorted(set(k for q, k in runs if q == p)):
        rounds = runs[(p, k)]
        fath = any(st == "fathomed" for _, _, st, _, _ in rounds)
        its = sum(i for *_, i, _ in rounds); cpu = sum(c for *_, c in rounds)
        best = max(lb for _, lb, *_ in rounds)
        print("%-8s  %4d  %-8s  %5d  %10d  %6.0f  %9.2f    (%d, %d)" % (
            p, k, "yes" if fath else "no", rounds[-1][0], its, cpu, best, *step2.get(k, (0, 0))))
print()
for p in protos:
    rs = [(k, runs[(p, k)]) for q, k in runs if q == p]
    fs = [(k, r) for k, r in rs if any(st == "fathomed" for _, _, st, _, _ in r)]
    if fs:
        its = sum(sum(i for *_, i, _ in r) for _, r in fs) / len(fs)
        cpu = sum(sum(c for *_, c in r) for _, r in fs) / len(fs)
        print("%-8s: %d of %d fathomed; fathoms mean %.0f iterations, %.0f CPU-s" % (p, len(fs), len(rs), its, cpu))
    else:
        print("%-8s: 0 of %d fathomed" % (p, len(rs)))
print("step 2  : 4 of 4 fathomed; mean %.0f iterations, %.0f CPU-s" % (
    sum(v[0] for v in step2.values()) / 4, sum(v[1] for v in step2.values()) / 4))
