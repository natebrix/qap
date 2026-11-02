#!/usr/bin/env python3
"""Tree model of ../wil50_feasibility/NOTES.md section 5, applied to a depth-3 sample.

A depth-3 node that converges s short of the cutoff fathoms k = ceil(s / g) levels lower
(g = median gain per level; nodes with s <= 0 fathom at depth 3). Its subtree holds
prod_{j=1..k} (n - 3 - j + 1) ... DNN nodes down to that level, every one of which runs the DNN
bound (misses at the last level are what fathom their children). Extra bound strength X per
node (a cut phase at every DNN node) is modelled as s -> s - X at every level. Node counts at
depth 3 are n(n-1)(n-2)/|G| for the symmetry group G. Estimates, not measurements.

usage: model.py n sym g cost_per_node [cut_cost_per_node] -- s1 s2 ... (shortfalls, negative = fathomed)
"""
import math, sys
a = sys.argv[1:]; i = a.index("--")
n, sym, g, cost = int(a[0]), int(a[1]), float(a[2]), float(a[3])
cutcost = float(a[4]) if i > 4 else cost
S = [float(x) for x in a[i+1:]]
d3 = n * (n - 1) * (n - 2) / sym
def subtree(s, X):
    """DNN nodes below one depth-3 node (the node itself excluded), and the depth it fathoms at"""
    k = 0 if s - X <= 0 else math.ceil((s - X) / g)
    nodes, width = 0, 1
    for j in range(1, k + 1):
        width *= (n - 3 - j + 1); nodes += width
    return nodes, 3 + k
print(f"n={n} |G|={sym} depth-3 nodes {d3:,.0f}; g={g} per level; {cost:.0f} CPU-s/node ({cutcost:.0f} with a cut phase)")
print(f"| extra bound per node | DNN nodes | CPU-s | CPU-years | 1,000 cores | fathoming depths (sample) |")
print("|---|---|---|---|---|---|")
for X, c in ((0, cost), (130, cutcost), (260, cutcost)):
    tot = d3 * (1 + sum(subtree(s, X)[0] for s in S) / len(S))
    depths = sorted(subtree(s, X)[1] for s in S)
    cpu = tot * c
    print(f"| +{X:.0f} | {tot:.1e} | {cpu:.1e} | {cpu / 3.156e7:,.0f} | {cpu / 1000 / 86400:,.0f} days | {' '.join(map(str, depths))} |")
