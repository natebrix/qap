"""Warm-rooted complete subtrees under a cut arm against the production trace's own subtree (same node, same
warm first solve), markdown.  python3 wtable.py
Production iterations are not in its trace: estimated [e] as nodes x the base arm's iterations per node here
(IT_PER_NODE, from t1-t3 base). Production CPU is c7a's; x 0.70 puts it on this box (calib.py)."""
import os, sys, glob, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
from summarize import stats, H
import ntrace
IT_PER_NODE = 25.6   # t1-t3 base, cold roots excluded: (2626+6187+7281-720) / (106+211+287-3)
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
idx = {",".join(fx): i for i, fx in enumerate(t.fixes)}
print("| node (warm root) | production: shortfall, nodes, CPU s (x 0.70 here) | arm | nodes | fathomed: first solve / phase | phases (fathomed) | iterations: first + phases = total | CPU s here | nodes vs production | iterations vs production [e] | CPU vs production x 0.70 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for d in sorted(glob.glob(os.path.join(H, "runs", "w[0-9]*_*"))):
    s, nodes, dnn = stats(d)
    if not nodes: continue
    i = idx[nodes[0]["fixes"]]
    pn, pc, ps = t.size[i], t.sub_cpu[i], t.cutoff[i] - t.dnn[i]
    done = s["done"]
    print("| `%s` | %.0f, %s, %s (%s) | %s | %d%s | %d / %d | %d (%d) | %s + %s = %s | %.0f | %.1fx | %.1fx | %.1fx |" % (
        nodes[0]["fixes"], ps, format(pn, ","), format(round(pc), ","), format(round(pc * 0.70), ","), os.path.basename(d).split("_")[1],
        s["nodes"], "" if done else " (RUNNING)", s["D"], s["C"], s["phases"], s["phase_fath"], format(s["it_first"], ","), format(s["it_cut"], ","),
        format(s["it_total"], ","), s["cpu"], pn / max(s["nodes"], 1), pn * IT_PER_NODE / max(s["it_total"], 1), pc * 0.70 / max(s["cpu"], 1e-9)))
