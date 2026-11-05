#!/usr/bin/env python3
"""Per-instance node ratios and root-bound differences of bench results
against a reference results.csv (default bench/reference/results.csv).

  python3 ratio.py RESULTS_DIR [REFERENCE_DIR]
"""
import csv
import math
import os
import sys

res = sys.argv[1]
ref = sys.argv[2] if len(sys.argv) > 2 else "bench/reference"
load = lambda d: {r["instance"]: r for r in csv.DictReader(open(os.path.join(d, "results.csv")))}
a, b = load(ref), load(res)
logs, worst_lo, worst_hi, tlogs = [], (1e9, ""), (0, ""), []
print("%-8s %9s %9s %7s %12s %12s %9s" % ("instance", "nodes ref", "nodes", "ratio", "root ref", "root", "rel diff"))
for k in a:
    if k not in b or a[k]["status"] != "ok" or b[k]["status"] != "ok":
        continue
    na, nb = int(a[k]["nodes"]), int(b[k]["nodes"])
    ra, rb = float(a[k]["root_bound"]), float(b[k]["root_bound"])
    r = nb / na
    logs.append(math.log(r))
    worst_lo = min(worst_lo, (r, k)); worst_hi = max(worst_hi, (r, k))
    ta, tb = float(a[k]["solver_time"]), float(b[k]["solver_time"])
    if ta > 1.0:
        tlogs.append(math.log(tb / ta))
    print("%-8s %9d %9d %7.3f %12.6g %12.6g %9.2e" % (k, na, nb, r, ra, rb, (rb - ra) / abs(ra)))
print("geometric mean node ratio %.3f over %d instances (min %.3f %s, max %.3f %s)" %
      (math.exp(sum(logs) / len(logs)), len(logs), worst_lo[0], worst_lo[1], worst_hi[0], worst_hi[1]))
if tlogs:
    print("geometric mean solver-time ratio %.3f over %d instances with reference time > 1 s "
          "(only meaningful when both ran on the same machine)" % (math.exp(sum(tlogs) / len(tlogs)), len(tlogs)))
