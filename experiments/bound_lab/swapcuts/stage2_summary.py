#!/usr/bin/env python3
"""Stage 2 table from logs/stage2_*.out (exch.C).  python3 stage2_summary.py > logs/stage2_summary.txt"""
import glob, os, statistics as st
here = os.path.dirname(os.path.abspath(__file__))
tot = {}
print("affected = v_ex > v_node (or no exchange-optimal completion below v_node + 5% of the optimum = inf);")
print("ex1 = fixed-free exchanges only (the plain pairwise cuts), ex2 = + fixed-fixed, 2opt = every exchange (reference).")
print("rise = (v_ex - v_node) / optimum, median over affected nodes with a finite value; off-path nodes only.\n")
print("%-8s %-4s %5s | %5s %5s %5s | %5s %5s %5s | %8s %8s" % ("inst", "kind", "depth", "nodes", "aff1", "inf1", "aff2", "inf2", "aff2o", "rise1%", "rise2%"))
for f in sorted(glob.glob(os.path.join(here, "logs", "stage2_*.out"))):
    inst = os.path.basename(f)[7:-4]
    if inst == "summary":
        continue
    rows = []
    opt = None
    for line in open(f):
        if line.startswith("#"):
            if "optimum" in line:
                opt = float(line.split("optimum")[1].split()[0])
            continue
        t = line.split()
        if len(t) != 11:
            continue
        v = [float("inf") if x == "inf" else float(x) for x in t[6:10]]
        rows.append((t[0], int(t[1]), v, int(t[10])))
    for kind in "GW":
        for d in range(1, 6):
            sel = [r for r in rows if r[0] == kind and r[1] == d and not r[3]]
            if not sel:
                continue
            n = len(sel)
            a1 = sum(r[2][1] > r[2][0] for r in sel); i1 = sum(r[2][1] == float("inf") for r in sel)
            a2 = sum(r[2][2] > r[2][0] for r in sel); i2 = sum(r[2][2] == float("inf") for r in sel)
            a3 = sum(r[2][3] > r[2][0] for r in sel)
            r1 = [100 * (r[2][1] - r[2][0]) / opt for r in sel if r[2][0] < r[2][1] < float("inf")]
            r2 = [100 * (r[2][2] - r[2][0]) / opt for r in sel if r[2][0] < r[2][2] < float("inf")]
            print("%-8s %-4s %5d | %5d %5d %5d | %5d %5d %5d | %8s %8s" % (inst, kind, d, n, a1, i1, a2, i2, a3,
                  "%.2f" % st.median(r1) if r1 else "-", "%.2f" % st.median(r2) if r2 else "-"))
            k = (kind, d)
            e = tot.setdefault(k, [0, 0, 0, 0, 0, 0])
            for x, y in enumerate((n, a1, i1, a2, i2, a3)):
                e[x] += y
print("\nall instances")
print("%-4s %5s | %5s %7s %7s | %7s %7s | %7s" % ("kind", "depth", "nodes", "aff1", "inf1", "aff2", "inf2", "aff2opt"))
for (kind, d), e in sorted(tot.items()):
    n = e[0]
    print("%-4s %5d | %5d %6.0f%% %6.0f%% | %6.0f%% %6.0f%% | %6.0f%%" % (kind, d, n, 100 * e[1] / n, 100 * e[2] / n, 100 * e[3] / n, 100 * e[4] / n, 100 * e[5] / n))
