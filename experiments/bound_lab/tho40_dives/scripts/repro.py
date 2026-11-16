#!/usr/bin/env python3
"""The reproduction run (runs/repro_s0: stratum 26:0 exactly as ../../tho40b_analysis/knuth_tho40/s0, the laptop's:
tho40-b's recipe, -u 240518, start level 3, q = 10, children by their predicted bound) against the laptop's."""
import os, re, statistics as st, sys
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAP = HERE + "/../../tho40b_analysis/knuth_tho40/s0"
def top(path):
    s = open(path + "/est.out", errors="replace").read()
    g = lambda p: (re.search(p, s) or [None, "?"])[1]
    return dict(sum=g(r"\nsum = (\S+)"), stack=g(r"\nstack = (\d+)"), fw=g(r"fw iterations: (\d+)"),
                d1=g(r"DNN bound by level:\n.*\n 1\s+\d+\s+\S+\s+\S+\s+(\d+)"), d2=g(r"DNN bound by level:\n.*\n.*\n 2\s+\d+\s+\S+\s+\S+\s+(\d+)"),
                t=g(r"Total Time: (\S+) seconds"))
def dives(path):
    return [tuple(map(float, l.split())) for l in open(path + "/dives.txt") if l.strip()]
a, b = top(LAP), top(HERE + "/runs/repro_s0")
print("| the top (the depth-1 node and its 39 children), deterministic | laptop | here |\n|---|---|---|")
for k, name in (("stack", "depth-3 nodes the dives start from"), ("sum", "sum of their scores"), ("fw", "FW iterations"), ("d1", "DNN iterations, depth 1"),
                ("d2", "DNN iterations, depth 2"), ("t", "CPU s")):
    print(f"| {name} | {a[k]} | {b[k]} |")
ta, tb = float(a["t"]), float(b["t"])
f = tb / ta
print(f"\nCPU here / laptop on the top's 40 nodes: {f:.2f}\n")
print("| dives | laptop | here |\n|---|---|---|")
A, B = dives(LAP), dives(HERE + "/runs/repro_s0")
row = lambda name, fa, fb: print(f"| {name} | {fa} | {fb} |")
row("dives", len(A), len(B))
for nm, i in (("node estimate", 1), ("CPU estimate, s (laptop s; here s / %.2f)" % f, 2)):
    sc = (1, 1) if i == 1 else (1, 1 / f)
    for stat, fn in (("median", st.median), ("mean", st.fmean), ("lower quartile", lambda v: st.quantiles(v, n=4)[0]), ("upper quartile", lambda v: st.quantiles(v, n=4)[2]), ("largest", max)):
        row(f"{nm}: {stat}", f"{fn([x[i] for x in A]) * sc[0]:,.0f}", f"{fn([x[i] for x in B]) * sc[1]:,.0f}")
dep = lambda V: " ".join(f"{d}: {sum(1 for x in V if x[0] == d) / len(V):.0%}" for d in range(3, 8))
row("last depth of a dive", dep(A), dep(B))
