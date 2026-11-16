#!/usr/bin/env python3
"""Part 2 sample: branched nodes of the production trace infra/phase1/runs/tho40ab-32/nodes.txt.gz
(subtree 26:10,8:0 of tho40t.dat) spread over the give-up shortfall range at depths 3-6, plus the
depth-2 node itself, two depth-1 nodes and the root. Writes nodes.tsv:
tag, depth, m, spec, production DNN bound, production shortfall, kids, subtree nodes, subtree cpu."""
import gzip, os, collections
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../.."))
CUT = 240514.05
rows = []
for l in gzip.open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/nodes.txt.gz"), "rt"):
    t = l.split()
    if t and t[0] == "node": rows.append(t)
# subtree size and cpu: every traced node whose fixes extend the node's
size = collections.Counter(); cpu = collections.Counter()
for t in rows:
    f = t[2].split(",")
    for k in range(2, len(f) + 1):
        size[",".join(f[:k])] += 1; cpu[",".join(f[:k])] += float(t[11])
B = collections.defaultdict(list)
for t in rows:
    if t[7] == "B": B[int(t[1])].append((CUT - float(t[4]), t))
pick = {3: [0.0, 0.25, 0.5, 0.75, 1.0], 4: [0.05, 0.5, 0.75, 0.9, 1.0], 5: [0.1, 0.5, 0.8, 1.0], 6: [0.2, 0.5, 0.8, 1.0]}
out = [("root", 0, "-", None), ("d1_path", 1, "26:2", None), ("d1_sub", 1, "26:10", None)]
def add(tag, t):
    out.append((tag, int(t[1]), t[2], t))
add("d2_sub", [t for t in rows if t[2] == "26:10,8:0"][0])
for d in sorted(pick):
    s = sorted(B[d], key=lambda x: x[0])
    for k, qn in enumerate(pick[d]):
        add("d%d_%s" % (d, "abcde"[k]), s[min(len(s) - 1, int(qn * (len(s) - 1) + 0.5))][1])
    if d == 4:
        add("d4_gate", [t for t in rows if t[2] == "26:10,8:0,39:32,38:7"][0])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "nodes.tsv"), "w") as fh:
    fh.write("tag\tdepth\tm\tspec\tprod_dnn\tprod_short\tkids\tsubtree_nodes\tsubtree_cpu\n")
    for tag, d, spec, t in out:
        if t is None: fh.write("%s\t%d\t%d\t%s\t-\t-\t-\t-\t-\n" % (tag, d, 40 - d, spec))
        else: fh.write("%s\t%d\t%d\t%s\t%.2f\t%.1f\t%s\t%d\t%.0f\n" % (tag, d, 40 - d, spec, float(t[4]), CUT - float(t[4]), t[10], size[spec], cpu[spec]))
