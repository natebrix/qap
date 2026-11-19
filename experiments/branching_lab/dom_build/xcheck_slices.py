#!/usr/bin/env python3
"""Every X node of a key-on slice trace, recomputed with dominance/dom.py (original data): best pair min Delta >= 1.
  python3 xcheck_slices.py INSTANCE ON.trace.gz"""
import gzip, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dominance"))
from qapdata import read_dat
from dom import Inst
name, tr = sys.argv[1:3]
inst = Inst(*read_dat(name))
vals = []
for l in gzip.open(tr, "rt"):
    t = l.split()
    if t[7] != "X": continue
    pm, pr = inst.pair_test(dict(tuple(map(int, s.split(":"))) for s in t[2].split(",")))
    vals.append(pm)
print(f"{name} {os.path.basename(tr)}: {len(vals)} X nodes, dom.py best pair min Delta >= 1 on {sum(v >= 1 for v in vals)}; min {min(vals) if vals else None}")
