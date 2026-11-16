#!/usr/bin/env python3
"""Every dive of every finished batch, one line each -> results/dives.csv (python3 scripts/dives_table.py).
Columns: arm, run, stratum, dive, last depth, nodes visited (unweighted), process CPU s (unweighted), then the dive's
estimates for its stratum's subtree below the depth-1 node (weighted by 1/probability): nodes, first-solve iterations,
phase iterations, CPU s here."""
import glob, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import parse, dive_tot, q, HERE
os.makedirs(HERE + "/results", exist_ok=True)
with open(HERE + "/results/dives.csv", "w") as f:
    f.write("arm,run,stratum,dive,last_depth,nodes_visited,process_cpu_s,est_nodes,est_first_solve_it,est_phase_it,est_cpu_s\n")
    for d in sorted(glob.glob(HERE + "/runs/*_s*_*")):
        m = re.match(r"(cut06|hand|comb|base)_s(\d+)_(\w+)$", os.path.basename(d))
        if not m or "exit 0" not in open(d + "/cmd.txt").read(): continue
        top, dives = parse(d + "/trace.txt")
        for i, dv in enumerate(dives):
            f.write(f"{m[1]},{os.path.basename(d)},26:{m[2]},{i},{max(n['depth'] for n in dv)},{len(dv)},{sum(q(n)['cpu'] for n in dv):.1f},"
                    f"{dive_tot(dv, 'nodes'):.0f},{dive_tot(dv, 'it1'):.0f},{dive_tot(dv, 'itc'):.0f},{dive_tot(dv, 'cpu'):.0f}\n")
print(open(HERE + "/results/dives.csv").read().count("\n") - 1, "dives")
