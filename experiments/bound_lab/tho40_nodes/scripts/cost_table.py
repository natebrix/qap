#!/usr/bin/env python3
"""Per sampled node: what the production run spent below it (trace CPU of its subtree less the
node itself, converted to DNN iterations of the children's depth by the run's own seconds per
iteration) against the CP iterations that fathom the node itself (x 1.2: a CP iteration costs
1.19 plain ones in dnncuts, logs/p2cp_*.log)."""
import json, os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ROOT = os.path.abspath(os.path.join(D, "../../.."))
import gzip
rows = [l.split() for l in gzip.open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/nodes.txt.gz"), "rt") if l.startswith("node ")]
own = {r[2]: float(r[11]) for r in rows}
sec_it = {}
txt = open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/results.txt")).read()
for l in txt.split("DNN bound by level:")[1].splitlines()[2:]:
    t = l.split()
    if len(t) >= 7 and t[0].isdigit(): sec_it[int(t[0])] = float(t[6]) / max(1, int(t[4]))
nodes = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(os.path.join(D, "nodes.tsv")).read().splitlines()[1:]}
print("| node | depth | production shortfall | subtree nodes | CPU below the node (s) | = DNN iterations at depth+1 | CP iterations to fathom | x 1.2 | iterations saved | ratio |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in json.load(open(os.path.join(D, "part2.json"))):
    if r["d"] < 3 or not r.get("cpfath"): continue
    t = nodes[r["tag"]]
    below = float(t[8]) - own[r["spec"]]
    its = below / sec_it[r["d"] + 1]
    cp = 1.2 * r["cpfath"]
    print("| %s | %d | %s | %s | %.0f | %.0f | %d | %.0f | %.0f | %.1f |" % (r["tag"], r["d"], t[5], t[7], below, its, r["cpfath"], cp, its - cp, its / cp))
