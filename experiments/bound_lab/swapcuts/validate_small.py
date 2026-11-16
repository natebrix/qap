#!/usr/bin/env python3
"""Stage 3 validity on small instances: every certified value of the exchange-cut solve
must stay at or below the node's optimum over exchange-optimal completions (v_ex1 for the
pairwise cuts, v_ex2 with -fp), which exch.C enumerated (logs/stage2_*.out).
  python3 validate_small.py [per_instance=20] > logs/validate_small.out
Also reports how much of the exact rise the relaxation captures, and how many nodes end
ABOVE the node's own optimum v_node (impossible for any inequality valid for all permutations).
"""
import os, re, subprocess, sys, random, statistics as st
here = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(here, "../../.."))
per = int(sys.argv[1]) if len(sys.argv) > 1 else 20
random.seed(5)
INF = float("inf")
worst_all = -INF
for inst in ("nug12", "scr12", "had12", "rou12", "tai12a", "tai12b", "chr12a"):
    rows = []; opt = None
    for line in open(os.path.join(here, "logs", "stage2_%s.out" % inst)):
        if line.startswith("#"):
            if "optimum" in line:
                opt = float(line.split("optimum")[1].split()[0])
            continue
        t = line.split()
        if len(t) != 11 or int(t[1]) not in (2, 3, 4):
            continue
        v = [INF if x == "inf" else float(x) for x in t[6:10]]
        spec = (t[3] + "," if t[3] != "-" else "") + "%s:%s" % (t[4], t[5])
        rows.append((spec, int(t[1]), v, int(t[10])))
    rows = list({r[0]: r for r in rows}.values())
    random.shuffle(rows)
    rows = rows[:per]
    for fp in (0, 1):
        worst = -INF; nchk = 0; above = 0; cap = []; ncert = 0; aff = 0
        for spec, d, v, onp in rows:
            vex = v[2] if fp else v[1]
            cmd = [os.path.join(here, "swapdrv"), "-x", spec, "-c", "1e12", "-I", "300", "-er", "8", "-ei", "150", "-ek", "300",
                   "-edelta", "%g" % (0.01 * opt), "-eq", "300", "-T", "1"] + (["-fp"] if fp else []) + [os.path.join(ROOT, "qaplib", inst + ".dat")]
            out = subprocess.run(cmd, capture_output=True, text=True).stdout
            lbs = [float(x) for x in re.findall(r"^(?:round +\d+|quiet round) +lb (-?[\d.]+)", out, re.M)]
            if not lbs:
                print("NO OUTPUT", inst, spec); continue
            ncert += len(lbs)
            plain, best = lbs[0], max(lbs)
            if vex < INF:
                nchk += 1
                worst = max(worst, best - vex)
                if vex > v[0]:
                    aff += 1
                    if vex - plain > 1e-6:
                        cap.append((best - plain) / (vex - plain))
            if best > v[0] + 1e-6 * abs(v[0]) + 1e-3:
                above += 1
        worst_all = max(worst_all, worst)
        print("%-7s %-9s nodes %2d  certified values %4d  checked against v_ex %2d: max(certified - v_ex) = %10.4f   "
              "ends above v_node: %2d   capture of the exact rise (median over %d affected finite nodes): %s" %
              (inst, "pair+fp" if fp else "pairwise", len(rows), ncert, nchk, worst, above, len(cap), "%.2f" % st.median(cap) if cap else "-"))
        sys.stdout.flush()
print("max over all: %.6f  %s" % (worst_all, "PASS" if worst_all <= 1e-6 else "FAIL"))
