#!/usr/bin/env python3
"""Are the best known solutions exchange-optimal (they must be, if they are optimal)?
  python3 check_bks.py > logs/check_bks.out"""
import os, sys, itertools
import numpy as np
from swaplib import read_dat, cost, swap_delta
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
for inst in ("tho40", "wil50", "sko42", "sko49", "tai50b", "nug30"):
    n, A, B = read_dat(os.path.join(ROOT, "qaplib", inst + ".dat"))
    t = open(os.path.join(ROOT, "QAPSOLN", inst + ".sln")).read().replace(",", " ").split()
    val = float(t[1]); p = [int(x) - 1 for x in t[2:2 + n]]
    q = [0] * n
    for i, x in enumerate(p): q[x] = i
    cands = [(name, s) for name, s in (("as given", p), ("inverted", q)) if abs(cost(A, B, s) - val) < 1e-6]
    if not cands:
        print("%-7s best known %.0f: neither the .sln permutation nor its inverse evaluates to it (%.0f, %.0f)" % (inst, val, cost(A, B, p), cost(A, B, q)))
        continue
    name, s = cands[0]
    deltas = [swap_delta(A, B, s, i, j) for i, j in itertools.combinations(range(n), 2)]
    # the solver symmetrizes and (for the t files) swaps A and B: the cost and hence every exchange gain is the same
    print("%-7s best known %.0f (%s): %d exchanges, smallest gain %+.0f, %d with gain 0, %d improving  -> %s" %
          (inst, val, name, len(deltas), min(deltas), sum(d == 0 for d in deltas), sum(d < 0 for d in deltas),
           "exchange-optimal" if min(deltas) >= 0 else "NOT exchange-optimal: not optimal"))
    if inst in ("tho40", "wil50"):
        print("        permutation (0-based, %s): %s" % (name, " ".join(map(str, s))))
