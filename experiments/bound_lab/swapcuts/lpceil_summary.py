#!/usr/bin/env python3
"""Table from logs/lpceil.out (also while lpceil.py is still running).  python3 lpceil_summary.py > logs/lpceil_summary.txt"""
import os, statistics as st
here = os.path.dirname(os.path.abspath(__file__))
OPT = {}
S = {}
INF = float("inf")
for line in open(os.path.join(here, "logs", "lpceil.out")):
    t = line.split()
    if line.startswith("#") or len(t) != 9:
        continue
    inst, d, m = t[0], int(t[1]), int(t[2])
    vnode, vlp = float(t[4]), float(t[5]); vex = INF if t[6] == "inf" else float(t[6])
    if vlp != vlp:
        vlp = INF      # the LP is infeasible: the linear cuts alone exclude the whole node
    nok, N = map(int, t[7].split("/"))
    onp = int(t[8])
    if onp:
        OPT[inst] = vnode
    S.setdefault((inst[-3:-1] if inst[-1] in "ab" else inst[-2:], d, m), []).append((inst, vnode, vlp, vex, nok / N, onp))
print("off-path nodes. logical = v_ex1 > v_node (stage 2's ceiling); linear = v_lp > v_node (the ceiling of the cuts as")
print("inequalities on Y over the exact lifted hull); kept = median of (v_lp - v_node) / (v_ex1 - v_node) where v_ex1 is finite")
print("and above v_node; x-opt = median share of the node's completions that are exchange-optimal against the fixed facilities.")
print("infeasible = the LP has no solution: no mixture of completions satisfies the cuts (the node is excluded outright).")
print("rise = median of (v_lp - v_node) / v_node in % over the nodes with a finite rise.")
print("%-4s %5s %3s %5s %8s %8s %10s %8s %6s %8s" % ("n", "depth", "m", "nodes", "logical", "linear", "infeasible", "rise%", "kept", "x-opt%"))
for (n, d, m), e in sorted(S.items()):
    e = [v for v in e if not v[5]]
    if not e:
        continue
    k = len(e)
    alog = sum(v[3] > v[1] + 1e-6 for v in e)
    alin = [v for v in e if v[2] > v[1] + 1e-6 * abs(v[1]) + 1e-6]
    # the optimum per instance: the smallest v_node seen (the on-path node's, when it was in the sample)
    rise = [100 * (v[2] - v[1]) / v[1] for v in alin if v[2] < INF]
    ninf = sum(v[2] == INF for v in e)
    kept = [(v[2] - v[1]) / (v[3] - v[1]) for v in e if v[1] + 1e-6 < v[3] < INF and v[2] < INF]
    print("%-4s %5d %3d %5d %7.0f%% %7.0f%% %9.0f%% %8s %6s %8.2f" % (n, d, m, k, 100 * alog / k, 100 * len(alin) / k, 100 * ninf / k,
          "%.2f" % st.median(rise) if rise else "-", "%.2f" % st.median(kept) if kept else "-", 100 * st.median(v[4] for v in e)))
