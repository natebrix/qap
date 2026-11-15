#!/usr/bin/env python3
"""Part 1 job list: the CP cut phase with a cutoff above the best known value on nodes of the
path to the best known solution (scripts/path.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import path as P
n, A, B = P.read_dat(os.path.join(P.ROOT, "qaplib/tho40t.dat"))
q0 = P.solution()
q = [x for x in (P.image(q0, fr, fc) for fr in (0, 1) for fc in (0, 1)) if x[26] == 2][0]
def job(name, order, d, iters):
    s = P.spec(q, order, d)
    print("p1_%s_d%02d\texperiments/dnn_cuts/dnncuts -x %s -c 241500 -grid 5:8 -side 1 -I 300 -k 1000 "
          "-cpsep 100:1000 -cpiters %d -T 1 -v qaplib/tho40t.dat" % (name, d, s, iters))
# the sharp ones first: nodes whose plain bound is within 1,000 of the optimum
for d in (10, 8, 9, 11, 12): job("warm", P.WARM, d, 2000)
for d in (10, 8): job("cold", P.COLD, d, 2000)
for d in (10, 8): job("index", P.INDEX, d, 2000)
for d in (6, 5, 4): job("warm", P.WARM, d, 1000)
for d in (6, 5, 4): job("index", P.INDEX, d, 1000)
for d in (3, 2, 1): job("warm", P.WARM, d, 1000)
for d in (3, 2, 1): job("index", P.INDEX, d, 600)
