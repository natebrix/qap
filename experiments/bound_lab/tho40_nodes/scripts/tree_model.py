#!/usr/bin/env python3
"""[e] What the measured CP fathoms would do to the production subtree 26:10,8:0 (trace
infra/phase1/runs/tho40ab-32/nodes.txt.gz, 22,650 nodes, 78,800 CPU s), by a threshold model:
a branched node whose give-up shortfall is at most T is fathomed by the CP phase instead, at
the CP iterations interpolated from part 2's sample (part2.json: production shortfall -> CP
iterations to pass the cutoff), each costing CPF plain iterations of that depth (the trace's own
seconds per DNN iteration by depth, results.txt). Its subtree disappears. Nodes above the
threshold are left as they were (no CP attempt, no inherited pool): conservative on both counts.
  python3 tree_model.py [CPF=1.5]
"""
import gzip, json, os, sys, collections, bisect
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ROOT = os.path.abspath(os.path.join(D, "../../.."))
CUT = 240514.05
CPF = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
rows = [l.split() for l in gzip.open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/nodes.txt.gz"), "rt") if l.startswith("node ")]
# seconds per DNN iteration by depth, from the run's "DNN bound by level" table (level = depth)
sec_it = {}
txt = open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/results.txt")).read()
for l in txt.split("DNN bound by level:")[1].splitlines()[2:]:
    t = l.split()
    if len(t) >= 7 and t[0].isdigit(): sec_it[int(t[0])] = float(t[6]) / max(1, int(t[4]))
# the sample: (production shortfall, CP iterations to fathom) for nodes that fathomed
pts = sorted((float(r["pshort"]), r["cpfath"]) for r in json.load(open(os.path.join(D, "part2.json")))
             if r.get("cpfath") and r["pshort"] != "-" and r["d"] >= 3)
xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
# monotone envelope (running maximum), linear interpolation
for i in range(1, len(ys)): ys[i] = max(ys[i], ys[i - 1])
def cp_iters(s):
    k = bisect.bisect_left(xs, s)
    if k == 0: return ys[0]
    if k >= len(xs): return None
    return ys[k - 1] + (ys[k] - ys[k - 1]) * (s - xs[k - 1]) / (xs[k] - xs[k - 1])
print("sample (production shortfall -> CP iterations, running maximum):", " ".join("%.0f:%d" % (x, y) for x, y in zip(xs, ys)))
print("seconds per DNN iteration by depth (production):", " ".join("%d:%.3f" % kv for kv in sorted(sec_it.items())))
tot_n = len(rows); tot_cpu = sum(float(r[11]) for r in rows)
print("production: %d nodes, %.0f CPU s; CP iteration = %.1f plain iterations" % (tot_n, tot_cpu, CPF))
print("| threshold T | give-up nodes fathomed by CP | nodes left | node ratio | CPU left | CP CPU added | CPU ratio |")
print("|---|---|---|---|---|---|---|")
for T in (250, 500, 1000, 1500, max(xs)):
    cut_at = set(); cp_cpu = 0.0
    for r in sorted(rows, key=lambda r: int(r[1])):
        if r[7] != "B": continue
        f = r[2].split(",")
        if any(",".join(f[:k]) in cut_at for k in range(2, len(f))): continue
        s = CUT - float(r[4])
        if s <= T and cp_iters(s) is not None:
            cut_at.add(r[2]); cp_cpu += cp_iters(s) * CPF * sec_it.get(int(r[1]), 0.15)
    n = 0; cpu = 0.0
    for r in rows:
        f = r[2].split(",")
        if any(",".join(f[:k]) in cut_at for k in range(2, len(f))): continue
        n += 1; cpu += float(r[11])
    print("| %.0f | %d | %d | %.3f | %.0f | %.0f | %.3f |" % (T, len(cut_at), n, n / tot_n, cpu, cp_cpu, (cpu + cp_cpu) / tot_cpu))
