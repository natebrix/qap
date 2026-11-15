"""The README's table: subtree x arm, markdown.  python3 table.py"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import stats, H
PROD = {"gate": (36, 406, "d4"), "t1": (106, 796, "d4"), "t2": (211, 1931, "d4"), "t3": (287, 2427, "d3"), "t4": (688, 2675, "d4")}
# t4 base was not rerun here: production 688 nodes (warm root); iterations [e] = 688 x 25.6 (wtable.py) + the cold root's 200
PRODBASE = {"t4": (688, 688 * 25.6)}
order = ["base", "z003", "z010", "z030", "z030g", "z030c", "zall"]
print("| subtree | arm | nodes | branched | fathomed: first solve / phase | phases (fathomed) | iterations: first solves + phases = total | of which the cold root | CPU s | iterations vs base | nodes vs base |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for t in ["gate", "t1", "t2", "t3", "t4"]:
    base = None
    for arm in order:
        d = os.path.join(H, "runs", "%s_%s" % (t, arm))
        if not os.path.exists(os.path.join(d, "done")): continue
        s, _, _ = stats(d)
        if arm == "base": base = s
        r_it = "%.2fx" % (base["it_total"] / s["it_total"]) if base and arm != "base" else ""
        r_n = "%.1fx" % (base["nodes"] / s["nodes"]) if base and arm != "base" else ""
        if not base and arm != "base" and t in PRODBASE:   # no base run here: the production trace's subtree [e]
            pn, pit = PRODBASE[t]
            r_it = "%.2fx [e]" % (pit / s["it_total"]); r_n = "%.1fx" % (pn / s["nodes"])
        print("| %s | %s | %d | %d | %d / %d | %d (%d) | %s + %s = %s | %s | %.0f | %s | %s |" % (
            t, arm, s["nodes"], s["B"], s["D"], s["C"], s["phases"], s["phase_fath"], format(s["it_first"], ","), format(s["it_cut"], ","),
            format(s["it_total"], ","), format(s["root_it"], ","), s["cpu"], r_it, r_n))
