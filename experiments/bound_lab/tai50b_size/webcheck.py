"""The web model's rule against the nodes actually bounded (every DNN trace from the long root, depths 1-8):
for each node, delta = (best solution found inside it by web.py's descent) - best known value, against g(depth),
the spine's shortfall at that depth.  python3 webcheck.py  -> logs/webcheck.txt"""
import glob, os, sys, collections, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, web as W
H = os.path.dirname(os.path.abspath(__file__))
G = [None, 3222232, 2652428, 2185798, 1029797, 626571, 199906, 112268, 107349, 107349]
fs = [f for f in glob.glob(H + "/runs/**/dnn.txt", recursive=True) if not any(x in f for x in ("/t30_", "/s1_", "/reh/", "/root_"))]
seen = {}; rows = []
for d in N.load(fs):
    key = tuple(d["fixes"])
    if key in seen or not 1 <= d["depth"] <= 9: continue
    seen[key] = 1
    fx = {int(x.split(":")[0]): int(x.split(":")[1]) for x in d["fixes"]}
    p = W.BK.copy()
    for r, l in fx.items():
        j = int(np.where(p == l)[0][0]); p[r], p[j] = p[j], p[r]
    assert all(p[r] == l for r, l in fx.items())
    q, v = W.descent(p, [x for x in range(50) if x not in fx])
    rows.append((d["depth"], ",".join(d["fixes"]), v - W.V0, G[d["depth"]], d["out"], d["fit"] + d["cit"], d["cutoff"] - d["final"], d["cstatus"]))
out = ["%d distinct nodes bounded at depths 1-9 under the long root" % len(rows), "",
       "| delta / g(depth) | nodes | branched | of those: at the phase limit or stalled within 1.5M (near misses) | fathomed | mean iterations of the fathomed |", "|---|---|---|---|---|---|"]
bins = [(0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 1e9)]
for lo, hi in bins:
    R = [r for r in rows if lo <= r[2] / r[3] < hi]; Bn = [r for r in R if r[4] == "B"]; F = [r for r in R if r[4] != "B"]
    out.append("| %s to %s | %d | %d | %d | %d | %s |" % (lo, hi if hi < 1e9 else "inf", len(R), len(Bn), sum(1 for r in Bn if r[6] < 1.5e6), len(F), "%.0f" % (sum(r[5] for r in F) / len(F)) if F else "-"))
out += ["", "branched nodes: depth, fixes, delta, g(depth), iterations, final shortfall, phase end"]
out += ["  %d %s delta %.0f g %d it %d short %.0f %s" % (r[0], r[1], r[2], r[3], r[5], r[6], r[7]) for r in sorted(rows) if r[4] == "B"]
out += ["", "fathomed nodes with delta < 1.5 g:"] + ["  %d %s delta %.0f g %d it %d" % (r[0], r[1], r[2], r[3], r[5]) for r in sorted(rows) if r[4] != "B" and r[2] < 1.5 * r[3]]
open(H + "/logs/webcheck.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
