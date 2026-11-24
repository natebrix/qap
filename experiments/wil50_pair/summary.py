#!/usr/bin/env python3
"""experiments/wil50_pair/summary.py: every node with a finished subtree in both arms, rule 10 against rule 9:
the line each arm branched the depth-2 node on, nodes, CPU, and rule 10's share of rule 9's."""
import os, re, sys
D = os.path.dirname(os.path.abspath(__file__))


def arm(d):
    res = open(os.path.join(d, "results.txt")).read()
    cpu = float(re.search(r"Total Time:\s*([\d.]+)", res).group(1))
    nodes = int(re.search(r"Nodes:\s*(\d+)", res).group(1))
    kids = set()                      # the depth-3 nodes' fix that the depth-2 node lacks
    fixes = [set(m.group(1).split()) for m in re.finditer(r"depth \d+ fixes (.*?) \|", open(os.path.join(d, "trace.txt")).read())]
    root = min(fixes, key=len)
    for f in fixes:
        if len(f) == len(root) + 1 and root < f:
            kids.add(next(iter(f - root)))
    fac = {k.split(":")[0] for k in kids}; loc = {k.split(":")[1] for k in kids}
    line = "r" + fac.pop() if len(fac) == 1 and len(kids) > 1 else "c" + loc.pop() if len(loc) == 1 and len(kids) > 1 else "?"
    return line, nodes, cpu


rows = []; tot = [0, 0.0, 0, 0.0]
names = sorted(n for n in os.listdir(os.path.join(D, "r10", "subtrees")) if os.path.isdir(os.path.join(D, "r10", "subtrees", n))) \
    if os.path.isdir(os.path.join(D, "r10", "subtrees")) else []
print("| node | rule 10: line, nodes, CPU s | rule 9: line, nodes, CPU s | nodes, rule 10 / rule 9 | CPU |\n|---|---|---|---|---|")
for n in names:
    a, b = (os.path.join(D, x, "subtrees", n) for x in ("r10", "r9"))
    if not all(os.path.isfile(os.path.join(x, f)) for x in (a, b) for f in ("results.txt", "trace.txt")):
        continue
    (la, na, ca), (lb, nb, cb) = arm(a), arm(b)
    if not na or not nb:
        print(f"| {n} | {na:,} nodes | {nb:,} nodes | an arm ran nothing: not a result | |")
        continue
    node = n.replace("_", ",").replace("-", ":")
    print(f"| {node} | {la}, {na:,}, {ca:,.0f} | {lb}, {nb:,}, {cb:,.0f} | {na / nb:.2f} | {ca / cb:.2f} |" + (" same line" if la == lb != "?" else ""))
    tot[0] += na; tot[1] += ca; tot[2] += nb; tot[3] += cb
if tot[2]:
    print(f"| all | {tot[0]:,}, {tot[1]:,.0f} | {tot[2]:,}, {tot[3]:,.0f} | {tot[0] / tot[2]:.2f} | {tot[1] / tot[3]:.2f} |")
else:
    print("no node is finished in both arms yet")
