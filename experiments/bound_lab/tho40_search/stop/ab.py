"""The stop-rule A/B on warm-rooted subtrees: every arm run on each node, iterations against z030, zall and
production (nodes x 25.6 [e], as ../wtable.py).  python3 ab.py"""
import os, sys, glob
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, H); sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
from summarize import stats
import ntrace
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
idx = {",".join(fx): i for i, fx in enumerate(t.fixes)}
ORDER = ["z030", "zall", "zallh", "zallhg", "zallhge", "zallhw", "zallhwg", "zallhws"]
print("| node | arm | nodes | fathomed: first solve / phase | phases (fathomed) | iterations: first + phases = total | root: first + phase (gap left) | CPU s | vs z030 | vs zall | vs zallhge | vs production [e] |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for w in ["w18", "w38", "w30"]:   # (w18r_zallhg, the attribution run, is in the text)
    rows = {}
    for d in glob.glob(os.path.join(H, "runs", w + "_*")):
        arm = os.path.basename(d).split("_", 1)[1]
        s, nodes, dnn = stats(d)
        if nodes: rows[arm] = (s, nodes, dnn)
    if not rows: continue
    any_nodes = next(iter(rows.values()))[1]
    i = idx[any_nodes[0]["fixes"]]
    pn = t.size[i]; pit = pn * 25.6
    print("| `%s` | production | %s | | | %s [e] | %.0f short | %s (c7a) | | | | |" % (any_nodes[0]["fixes"].split(",")[-1], format(pn, ","), format(round(pit), ","), t.cutoff[i] - t.dnn[i], format(round(t.sub_cpu[i]), ",")))
    for arm in ORDER:
        if arm not in rows: continue
        s, nodes, dnn = rows[arm]
        r = dnn[0]
        left = "fathomed" if r["cstatus"] == "fathomed" else "%.0f left" % (r["cutoff"] - r["cbest"]) if r["cut"] else "-"
        def ratio(a):
            return "%.2fx" % (rows[a][0]["it_total"] / s["it_total"]) if a in rows and a != arm and s["done"] and rows[a][0]["done"] else ""
        print("| | %s%s | %d | %d / %d | %d (%d) | %s + %s = %s | %d + %d (%s, from %.0f) | %.0f | %s | %s | %s | %s |" % (
            arm, "" if s["done"] else " (part-way)", s["nodes"], s["D"], s["C"], s["phases"], s["phase_fath"], format(s["it_first"], ","), format(s["it_cut"], ","),
            format(s["it_total"], ","), r["iters"], r["cit"], left, r["cutoff"] - r["lb"], s["cpu"], ratio("z030"), ratio("zall"), ratio("zallhge"), "%.2fx" % (pit / s["it_total"]) if s["done"] else ""))
