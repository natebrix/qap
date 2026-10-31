#!/usr/bin/env python3
"""A gap-ratio give-up rule, evaluated on -v traces: from iteration H on, give
up when the gap to the cutoff shrank by less than a factor rho over the last
50 iterations (g(t) / g(t-50) > rho, gap of the running best)."""
import re, sys
cut = 48814.05
def load(f):
    pts = [(int(m.group(1)), float(m.group(2))) for m in re.finditer(r"^it\s+(\d+) lb (\S+)", open(f).read(), re.M)]
    out = {}; b = -1e300
    for it, lb in pts: b = max(b, lb); out[it] = b
    return out
def run(best, H, rho):
    worst = 0
    for t in sorted(best):
        if best[t] > cut: return ("fathom", t, worst)
        if t < H or t - 50 not in best: continue
        ratio = (cut - best[t]) / (cut - best[t - 50])
        worst = max(worst, ratio)
        if ratio > rho: return ("give up", t, ratio)
    return ("never", None, worst)
fath = [22, 24, 27, 29]; miss = [1, 5, 6, 9]
for k in fath + miss:
    b = load(f"cloud/trace/node{k}.log")
    r = [round((cut - b[t]) / (cut - b[t - 50]), 2) for t in range(150, 451, 50) if t in b and t - 50 in b and b[t] <= cut]
    print(f"node {k:2d} gap ratio per 50 iterations from 150: {r}")
for H in (150, 200, 250):
    for rho in (0.6, 0.65, 0.7):
        F = [run(load(f"cloud/trace/node{k}.log"), H, rho) for k in fath]
        M = [run(load(f"cloud/trace/node{k}.log"), H, rho) for k in miss]
        print(f"H {H} rho {rho}: fathoms {[f[0] for f in F]} (max ratio {max(f[2] for f in F):.2f}); misses give up at {[m[1] for m in M]}")
