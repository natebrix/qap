"""A warm-rooted run's root and its children, node by node, against the production trace's same nodes:
production shortfall, subtree nodes and CPU | here: outcome, first solve (pooled) and phase.
  python3 wchildren.py RUN [RUN2]      e.g. w18_z030 [w18_zall]: RUN2's columns beside RUN's
Production iterations [e] = subtree nodes x 25.6 (wtable.py)."""
import sys, os
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ntrace
from summarize import load, H
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
idx = {",".join(fx): i for i, fx in enumerate(t.fixes)}
runs = sys.argv[1:]
data = []
for r in runs:
    n, d = load(os.path.join(H, "runs", r))
    data.append(({x["fixes"]: (x, e) for x, e in zip(n, d)}, n))
root = data[0][1][0]["fixes"]; ri = idx[root]
def cell(m, fx):
    if fx not in m: return "-"
    x, e = m[fx]
    s = "%s %d it" % (x["out"], e["iters"])
    if x["out"] != "D": s += " (%.0f short)" % (e["cutoff"] - e["lb"])
    if e["cut"]: s += " + %d it %s" % (e["cit"], e["cstatus"]) + ("" if e["cstatus"] == "fathomed" else " %.0f short" % (e["cutoff"] - e["cbest"]))
    return s
def its(m, fx):
    if fx not in m: return 0
    return m[fx][1]["iters"] + m[fx][1]["cit"]
print("| node | production: shortfall, subtree nodes, iterations [e] | " + " | ".join(runs) + " |")
print("|---|---|" + "---|" * len(runs))
print("| %s (root) | %.0f, %d, - | %s |" % (root, t.cutoff[ri] - t.dnn[ri], t.size[ri], " | ".join(cell(m, root) for m, _ in data)))
tot = [0] * len(runs); ptot = 0; pn = 0; seen = [0] * len(runs)
order = [x["fixes"] for x in data[0][1][1:]]
kids = sorted(t.children[ri], key=lambda c: order.index(",".join(t.fixes[c])) if ",".join(t.fixes[c]) in order else 999)
below = {}
for m, n in data:
    for x in n: pass
for c in kids:
    fx = ",".join(t.fixes[c]); s = t.cutoff[c] - t.dnn[c]
    if t.outcome[c] != "B" and all(fx not in m or m[fx][0]["out"] == "D" for m, _ in data): continue   # fathomed at once everywhere: summed below
    print("| %s | %.0f, %d, %s | %s |" % (fx[len(root) + 1:], s, t.size[c], format(round(t.size[c] * 25.6), ","), " | ".join(cell(m, fx) for m, _ in data)))
for k, (m, n) in enumerate(data):
    tot[k] = sum(e["iters"] + e["cit"] for x, e in m.values()); seen[k] = len(n)
print("| all nodes traced so far | %d nodes, %s it [e] | %s |" % (t.size[ri], format(round(t.size[ri] * 25.6), ","), " | ".join("%d nodes, %s it" % (seen[k], format(tot[k], ",")) for k in range(len(runs)))))
