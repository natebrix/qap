"""Replay a zone Z and a phase iteration cap T on the phases of a traced cut arm (default z010) [e].
  python3 replay.py [ARM]
Model, per phase of ARM on a node the base arm also visited (unpooled phases only, so the node is base's node):
  not run (shortfall > Z): base's cost below the node;
  run and fathomed within T: its iterations x W;
  otherwise: T x W + base's cost below (the pool is assumed to help the children not at all: pessimistic).
W = 1.07, the measured cost of a phase iteration in first-solve iterations. Output: iterations saved against base (110,926).
"""
import sys, collections
from summarize import load
arm = sys.argv[1] if len(sys.argv) > 1 else "z010"
_, b = load("base"); _, c = load(arm)
key = lambda e: frozenset(e["fixes"].split())
bset = collections.defaultdict(set)
for e in b: bset[e["depth"]].add(key(e))
below = collections.Counter()
for e in b:
    k = key(e)
    for d in range(1, e["depth"]):
        for a in bset[d]:
            if a < k: below[a] += e["iters"]
ph = [(e["cutoff"] - e["lb"], e["cit"] if e["cstatus"] == "fathomed" else 10**9, below[key(e)]) for e in c if e["cut"] and not e["kept"] and key(e) in bset[e["depth"]]]
W = 1.07; tot = sum(e["iters"] for e in b)
Ts = [100, 200, 300, 400, 600, 1000, 3400]
print("unpooled matched phases:", len(ph), " base iterations:", tot)
print("saved iterations (percent of base) | T = " + "  ".join("%6d" % t for t in Ts))
for Z in [5, 8, 10, 12, 15, 18.4, 25, 30, 40, 61.3]:
    row = []
    for T in Ts:
        s = 0
        for short, it, bl in ph:
            if short > Z: continue
            s += bl - (it * W if it <= T else T * W + bl)
        row.append("%5.1f%%" % (100.0 * s / tot))
    print("Z = %5.1f (cut.zone %.4f)          | " % (Z, Z / 6124) + "  ".join(row))
