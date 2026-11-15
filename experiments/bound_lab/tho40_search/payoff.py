"""Each fathoming phase of a cut arm against what the base arm spent under the same node (matched by fixes):
shortfall, phase iterations, base nodes and DNN iterations below the node.   python3 payoff.py"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, H
print("run        depth  short  phase it | base below: nodes  iterations | iterations saved per phase iteration")
for d in sorted(glob.glob(os.path.join(H, "runs", "t*_z*")) + glob.glob(os.path.join(H, "runs", "gate_z*"))):
    name = os.path.basename(d); t = name.split("_")[0]
    b = os.path.join(H, "runs", t + "_base")
    if not os.path.exists(os.path.join(b, "done")): continue
    bn, bd = load(b)
    below = {}
    for x, e in zip(bn, bd):
        f = x["fixes"].split(",")
        for k in range(1, len(f)):
            p = ",".join(f[:k]); v = below.setdefault(p, [0, 0]); v[0] += 1; v[1] += e["iters"]
    cn, cd = load(d)
    for x, e in zip(cn, cd):
        if e["cut"] and e["cstatus"] == "fathomed" and x["fixes"] in below:
            v = below[x["fixes"]]
            print("%-10s %3d %7.1f %6d   | %5d %8d | %.1f%s" % (name, x["depth"], e["cutoff"] - e["lb"], e["cit"], v[0], v[1], v[1] / e["cit"], "  (cold root)" if x is cn[0] else ""))
