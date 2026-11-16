#!/usr/bin/env python3
"""Sanity of the dive runs (python3 scripts/sanity.py):
1. no run printed a new objective value (the incumbent is the best known 240,516, or 240,518 in path_s2/repro_s0);
2. every node seen whose fixes agree with the best known solution (the image with 26 at 2, ../tho40_nodes/scripts/path.py)
   has a bound of at most 240,516 (its optimum), in every arm;
3. the path run (runs/path_s2: the cut arm, incumbent 240,518, the dive made to follow the best known solution):
   every node on it must branch, with its first solve and its phase at most 240,516."""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import parse, HERE
P = [27, 38, 32, 14, 5, 18, 34, 10, 15, 26, 25, 17, 4, 30, 28, 20, 24, 23, 37, 29, 11, 13, 36, 35, 8, 3, 2, 22, 9, 19, 0, 12, 21, 33, 1, 16, 6, 7, 31, 39]
OPT = 240516.0
ok = True
runs = sorted(glob.glob(HERE + "/runs/*"))
new = [r for r in runs if "New objective" in open(r + "/est.out", errors="replace").read()]
print(f"1. runs with a 'New objective value' line: {len(new)} of {len(runs)} {[os.path.basename(r) for r in new]}")
ok &= not new
onpath = lambda nd: all(P[int(a)] == int(b) for a, b in (x.split(":") for x in nd["fixes"].split()))
seen = bad = 0
for r in runs:
    if not os.path.exists(r + "/trace.txt"): continue
    top, dives = parse(r + "/trace.txt")
    for nd in top + [n for d in dives for n in d]:
        if onpath(nd):
            seen += 1
            vals = [v for v in (nd["bound"], nd["lb1"], nd["cbest"]) if v is not None]
            if max(vals) > OPT + 1e-6:
                bad += 1; print("   ABOVE THE OPTIMUM:", os.path.basename(r), nd["depth"], nd["fixes"], vals)
print(f"2. nodes on the path to the best known solution seen in any run: {seen}; with a bound above 240,516: {bad}")
ok &= bad == 0
r = HERE + "/runs/path_s2"
if os.path.exists(r + "/trace.txt"):
    top, dives = parse(r + "/trace.txt")
    print("3. the path run (cut arm, incumbent 240,518):")
    for nd in top + [n for d in dives for n in d]:
        good = onpath(nd) and nd["fate"] == "B" and max(v for v in (nd["bound"], nd["lb1"], nd["cbest"]) if v is not None) <= OPT + 1e-6
        ok &= good
        ph = f"phase best {nd['cbest']:.4f} after {nd['itc']} it ({nd['cstat']})" if nd["cut"] else "no phase"
        print(f"   depth {nd['depth']:2d} {'on path' if onpath(nd) else 'OFF PATH'} fate {nd['fate']} kids {nd['kids']:2d} first solve "
              f"{nd['lb1'] if nd['lb1'] is not None else float('nan'):.4f} ({nd['it1']} it, pool {nd['pool']}), {ph}, {nd['cpu']:.0f} s -> {'ok' if good else 'FAIL'}")
    miss = "DIVE_PATH: no child" in open(r + "/est.out", errors="replace").read()
    if miss: print("   the path child was missing at some node: FAIL"); ok = False
print("PASS" if ok else "FAIL")
