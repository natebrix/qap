"""Matched nodes (same fixes, in order) between the production trace and the base runs here:
CPU ratio (this box, shared, / production c7a), node counts, and iterations per node by depth and outcome."""
import os, sys, glob, collections, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, H
prod = {}
for ln in gzip.open("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz", "rt"):
    f = ln.split()
    if len(f) >= 12 and f[0] == "node": prod[f[2]] = (f[7], float(f[11]), f[4], int(f[1]))
tot = collections.defaultdict(lambda: [0, 0.0, 0.0, 0, 0, 0])   # (depth,outcome) -> n, cpu here, cpu prod, iters, same outcome, n matched
for d in sorted(glob.glob(os.path.join(H, "runs", "*_base"))):
    nodes, dnn = load(d)
    n_in = sum(1 for x in nodes[1:] if x["fixes"] in prod)
    same = sum(1 for x in nodes[1:] if x["fixes"] in prod and prod[x["fixes"]][0] == x["out"])
    ch = sum(x["cpu"] for x in nodes[1:] if x["fixes"] in prod); cp = sum(prod[x["fixes"]][1] for x in nodes[1:] if x["fixes"] in prod)
    print("%-10s nodes here %4d, of which in production (root aside) %4d, same outcome %4d; CPU here %7.1f, production %7.1f, ratio %.2f" % (os.path.basename(d), len(nodes), n_in, same, ch, cp, ch / max(cp, 1e-9)))
    for x, e in zip(nodes, dnn):   # traces are in the same order (one DNN line per node here: QPB fathoms nothing)
        if x is nodes[0]: continue
        t = tot[(x["depth"], x["out"])]
        t[0] += 1; t[1] += x["cpu"]; t[3] += e["iters"]
        if x["fixes"] in prod: t[2] += prod[x["fixes"]][1]; t[5] += 1
print("depth outcome  n  it/node  cpu/node here  cpu/node production (matched)  s/iteration here")
for k in sorted(tot):
    t = tot[k]
    print("%4d %s %5d %7.1f %8.2f %8.2f  %.3f" % (k[0], k[1], t[0], t[3] / t[0], t[1] / t[0], t[2] / max(t[5], 1), t[1] / max(t[3], 1)))
