"""The A/B subtree's CPU under a cut arm, by the class of its depth-3 nodes (production CPU shares,
../prodstats.py), and the whole tree if it responds alike.  python3 tree.py

A class's factor is production CPU x 0.70 (this box per matched node, ../calib.py) over the CPU measured
here on the warm-rooted subtree named, or an assumption [e]. The A/B root's own solve and phase
(runs/wroot_ARM) is added: under zallh it runs 3,030 iterations.
"""
import os, sys
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, H); sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
from summarize import stats
import ntrace
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
idx = {",".join(fx): i for i, fx in enumerate(t.fixes)}
TOT = sum(t.cpu)
def measured(run):
    d = os.path.join(H, "runs", run)
    if not os.path.exists(os.path.join(d, "done")): return None
    s, nodes, dnn = stats(d)
    i = idx[nodes[0]["fixes"]]
    return t.sub_cpu[i] * 0.70 / s["cpu"], t.size[i] * 25.6 / s["it_total"], t.size[i] / s["nodes"]
def classes(arms):
    """(label, share of CPU, (low, mid, high) factor, source)"""
    def m(name):
        for a in arms:
            r = measured(name + "_" + a)
            if r: return r, name + "_" + a
        return None, None
    w38, n38 = m("w38"); w30, n30 = m("w30"); w18, n18 = m("w18")
    out = [("fathomed at once, or <= 1,500", 0.012, (1.5, 2.0, 4.0), "[e]")]
    out.append(("1,500-2,500", 0.074, (0.6 * w38[0], w38[0], w38[0]), "[m] %s (the farthest of them) %.1fx; low: 0.6 of it for the small ones" % (n38, w38[0])))
    if w30: out.append(("2,500-3,500", 0.123, (0.8 * w30[0], w30[0], 1.2 * w30[0]), "[m] %s %.1fx, +-20%%" % (n30, w30[0])))
    else: out.append(("2,500-3,500", 0.123, (min(w38[0], w18[0]), (w38[0] + w18[0]) / 2, max(w38[0], w18[0])), "[e] between the two measured"))
    out.append(("3,500-5,000", 0.537, (0.8 * w18[0], w18[0], 1.2 * w18[0]), "[m] %s %.1fx, +-20%% for the 8 other nodes" % (n18, w18[0])))
    out.append(("above 5,000 (2:32)", 0.253, (0.5 * w18[0], w18[0], 1.2 * w18[0]), "[e] as w18; low: half of it"))
    return out
def total(arms, root_run):
    cl = classes(arms)
    rs, _, _ = stats(os.path.join(H, "runs", root_run))
    root = rs["cpu"] / 0.70 / TOT          # the A/B root's solve and phase, as a share of production's CPU
    res = []
    for k in range(3):
        res.append(1.0 / (root + sum(sh / f[k] for _, sh, f, _ in cl)))
    return cl, root, res
if __name__ == "__main__":
    for label, arms, root_run in (("first study: z030 (its own class table gave 2.4x, 1.9-2.9)", None, None),
                                  ("zall + hand-over rule + gu_zone 0 + giveup_warm 100 (zallhge; zallhg where hge was not run)", ["zallhge", "zallhg"], "wroot_zallh")):
        if arms is None: continue
        cl, root, res = total(arms, root_run)
        print(label)
        print("| depth-3 class (short) | share of the subtree's CPU | CPU factor: low, taken, high | from |")
        print("|---|---|---|---|")
        for lab, sh, f, src in cl:
            print("| %s | %.1f%% | %.1fx, **%.1fx**, %.1fx | %s |" % (lab, 100 * sh, f[0], f[1], f[2], src))
        print("| the A/B root's own solve and phase | +%.1f%% of production's CPU | | [m] %s |" % (100 * root, root_run))
        print("\nA/B subtree: %.1fx less CPU (%.1f-%.1fx); whole tree if alike: %.0f c7a CPU-days (%.0f-%.0f) from 1,300" % (
            res[1], res[0], res[2], 1300 / res[1], 1300 / res[2], 1300 / res[0]))
