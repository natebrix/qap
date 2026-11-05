#!/usr/bin/env python3
"""Compare the tail-probe arms (runs/ARM/summary.txt): estimated DNN iterations
and CPU per level and in total, and the ratios against base and baseship.
  python3 compare.py [RUNS]"""
import os
import re
import sys

runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
arms = [a for a in ("base", "baseship", "full", "fullwide") if os.path.exists(os.path.join(runs, a, "summary.txt"))]
data = {}
for a in arms:
    txt = open(os.path.join(runs, a, "summary.txt")).read()
    lv = {}
    for m in re.finditer(r"^\| (\d+) \| ([\d,]+) \| (\d+) \| (\d+) \| (\d+) \| ([\d,]+) .*? \| ([\d.,]+) \| ([\d,]+) \| ([\d,]+) \|", txt, re.M):
        n = lambda x: float(x.replace(",", ""))
        lv[int(m.group(1))] = dict(n=n(m.group(2)), bounded=int(m.group(3)), fathomed=int(m.group(4)),
                                   branched=int(m.group(5)), it=n(m.group(8)), cpu=n(m.group(9)))
    below = re.search(r"Open below level (\d+) \(estimated\): ([\d,]+) nodes", txt)
    data[a] = dict(levels=lv, below=(int(below.group(1)) + 1, float(below.group(2).replace(",", ""))) if below else None)
print("| arm | level: est. nodes, fathomed/bounded, est. iterations | total est. iterations | total est. CPU s | open below |")
print("|---|---|---|---|---|")
for a in arms:
    L = data[a]["levels"]
    cells = "; ".join(f"L{k}: {v['n']:,.0f}, {v['fathomed']}/{v['bounded']}, {v['it']:,.0f}" for k, v in sorted(L.items()))
    ti = sum(v["it"] for v in L.values()); tc = sum(v["cpu"] for v in L.values())
    b = data[a]["below"]
    print(f"| {a} | {cells} | {ti:,.0f} | {tc:,.0f} | {f'{b[1]:,.0f} at L{b[0]}' if b else '-'} |")
for ref in ("base", "baseship"):
    if ref not in data:
        continue
    R = data[ref]
    ti_r = sum(v["it"] for v in R["levels"].values())
    for a in arms:
        if a == ref:
            continue
        ti = sum(v["it"] for v in data[a]["levels"].values())
        br = data[a]["below"][1] if data[a]["below"] else 0
        br_r = R["below"][1] if R["below"] else 0
        print(f"- {a} against {ref}: {ti_r / ti:.2f}x fewer estimated iterations through the probed levels; "
              f"open below: {br:,.0f} against {br_r:,.0f}" + (f" ({br_r / br:.2f}x)" if br else ""))
