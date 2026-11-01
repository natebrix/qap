#!/usr/bin/env python3
"""Depth-4 sample under the clear misses of all32.txt: for each parent (a
depth-3 node the DNN bound gives up on, 158-203 short at convergence), one
random free facility and four random free locations, as in the depth-4
study of ../NOTES.md section 3 (seed there 20261030).

  python3 make_sample.py > nodes.txt      # lines: <id> <parent id> <spec>
"""
import random, sys

N = 50
PARENTS = {10: "39:27,30:48,25:49", 13: "21:28,8:41,28:19", 2: "18:43,11:29,46:44", 7: "12:7,33:6,32:11"}
rng = random.Random(20261101)
k = 1
for pid, spec in PARENTS.items():
    fixes = [tuple(map(int, s.split(":"))) for s in spec.split(",")]
    locs = set(range(N)) - {l for l, _ in fixes}
    facs = set(range(N)) - {f for _, f in fixes}
    f = rng.choice(sorted(facs))
    for l in rng.sample(sorted(locs), 4):
        print("%d %d %s,%d:%d" % (k, pid, spec, l, f))
        k += 1
