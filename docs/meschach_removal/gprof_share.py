#!/usr/bin/env python3
"""Sum a gprof flat profile's self time by origin: Meschach (functions
defined in src/MESCHACH/meschach.a), the solver's eigen wrapper and the
Meschach-shaped helpers in linalg.C, and the rest.

  python3 gprof_share.py FLAT.txt MESCHACH.a

gprof samples only the executable's own text, so time in libc (malloc,
free, memcpy) and in OpenBLAS is not in the profile; the script prints the
sampled total so the caller can compare it with the run's CPU time.
"""
import re
import subprocess
import sys

flat, archive = sys.argv[1], sys.argv[2]
nm = subprocess.run(["nm", "-g", "--defined-only", archive], capture_output=True, text=True).stdout
mes = {l.split()[2] for l in nm.splitlines() if len(l.split()) == 3 and l.split()[1] in "TtWw"}

row = re.compile(r"^\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(?:(\d+)\s+[\d.]+\s+[\d.]+\s+)?(.+?)\s*$")
total = 0.0
groups = {"meschach": 0.0, "dlap (LAP solver)": 0.0, "fw (Frank-Wolfe)": 0.0, "other solver": 0.0}
mes_rows = []
for line in open(flat):
    m = row.match(line)
    if not m:
        continue
    self_s = float(m.group(3))
    name = m.group(5)
    base = name.split("(")[0]
    total += self_s
    if base in mes:
        groups["meschach"] += self_s
        mes_rows.append((self_s, base))
    elif base == "dlap":
        groups["dlap (LAP solver)"] += self_s
    elif base == "fw":
        groups["fw (Frank-Wolfe)"] += self_s
    else:
        groups["other solver"] += self_s
print("sampled total %.2f s" % total)
for g, s in groups.items():
    print("  %-20s %8.2f s  %5.1f%%" % (g, s, 100 * s / total))
print("Meschach functions by self time:")
for s, n in sorted(mes_rows, reverse=True):
    if s > 0:
        print("  %8.2f  %s" % (s, n))
