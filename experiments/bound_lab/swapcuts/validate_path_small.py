#!/usr/bin/env python3
"""All cut families (pairwise, -fp, -rot) on small on-path nodes: a node that contains a global
optimum must never get a certified value above the optimum, whatever the multipliers.
  python3 validate_path_small.py > logs/validate_path_small.out      (uses ./swapdrv2 = the current swapdrv.C)
"""
import os, re, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(here, "../../.."))
worst_all = -1e300; total = 0
for inst in ("nug12", "scr12", "had12", "rou12", "tai12a", "tai12b", "chr12a"):
    rows = []; opt = None
    for line in open(os.path.join(here, "logs", "stage2_%s.out" % inst)):
        if line.startswith("#"):
            if "optimum" in line:
                opt = float(line.split("optimum")[1].split()[0])
            continue
        t = line.split()
        if len(t) == 11 and t[10] == "1" and int(t[1]) >= 2:
            rows.append((t[3] + "," if t[3] != "-" else "") + "%s:%s" % (t[4], t[5]))
    worst = -1e300; n = 0; above_plain = 0
    for spec in sorted(set(rows)):
        for frac in (0.0005, 0.01):
            cmd = [os.path.join(here, "swapdrv2"), "-x", spec, "-c", "1e12", "-I", "400", "-er", "6", "-ei", "150", "-ek", "400",
                   "-edelta", "%g" % (frac * opt), "-rot", "-fp", "-T", "1", os.path.join(ROOT, "qaplib", inst + ".dat")]
            out = subprocess.run(cmd, capture_output=True, text=True).stdout
            lbs = [float(x) for x in re.findall(r"^(?:round +\d+|quiet round) +lb (-?[\d.]+)", out, re.M)]
            if not lbs or "max |err| 0" not in out:
                print("PROBLEM", inst, spec, out[-300:]); continue
            n += len(lbs); worst = max(worst, max(lbs) - opt)
    total += n; worst_all = max(worst_all, worst)
    print("%-7s optimum %.0f: %d on-path nodes (depths 2-5), %d certified values with pairwise + fixed-pair + rotation cuts: max(certified - optimum) = %.6f" %
          (inst, opt, len(set(rows)), n, worst))
    sys.stdout.flush()
print("%d certified values, max over all %.6f: %s" % (total, worst_all, "PASS" if worst_all <= 1e-6 else "FAIL"))
