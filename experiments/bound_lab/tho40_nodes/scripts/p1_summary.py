#!/usr/bin/env python3
"""Part 1 table: every certified value of logs/p1_*.log against the best known value."""
import glob, os, re, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")
BEST = 240516.0
rows = []
for f in sorted(glob.glob(os.path.join(D, "p1_*.log"))):
    t = open(f).read()
    tag = os.path.basename(f)[:-4]
    m = re.search(r"-x (\S+)", t); spec = m.group(1)
    r0 = re.search(r"^round  0  lb (\S+)", t, re.M)
    certs = [(int(i), float(c)) for i, c in re.findall(r"^cp it +(\d+) cert (\S+) ", t, re.M)]
    lbs = [float(x) for x in re.findall(r"^it +\d+ lb (\S+)", t, re.M)]
    fin = re.search(r"^cp solve  lb (\S+) .*status (\S+)  iters (\d+)", t, re.M)
    done = os.path.exists(f[:-4] + ".done")
    if not certs: continue
    mx = max(c for _, c in certs); allmax = max([mx] + lbs + ([float(fin.group(1))] if fin else []))
    last = certs[-1]
    # the certified value 200 iterations before the end: is it still rising?
    prev = [c for i, c in certs if i <= last[0] - 200]
    rows.append((tag, spec.count(":"), 40 - spec.count(":"), float(r0.group(1)) if r0 else float("nan"), last[0], last[1],
                 last[1] - prev[-1] if prev else float("nan"), allmax, len(certs), "done" if done else "running",
                 "PASS" if allmax <= BEST else "FAIL"))
print("| run | depth | m | plain@300 | CP iters | CP certified (last) | gain over last 200 it | max certified, all checks | 240,516 - max | checks | |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (r[0].split("_")[1], r[1])):
    print("| %s | %d | %d | %.1f | %d | %.4f | %+.1f | %.4f | %.4f | %d | %s%s |" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], BEST - r[7], r[8], r[10], "" if r[9] == "done" else " (running)"))
print("\nworst: max certified - 240516 = %.6f over %d certified values in %d runs" % (max(r[7] for r in rows) - BEST, sum(r[8] for r in rows), len(rows)))
