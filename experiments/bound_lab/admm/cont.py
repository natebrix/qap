#!/usr/bin/env python3
"""Same-node continuation: run the unscaled PRSM for N0 iterations, then continue from that state
(Y, Z in original units) for N1 more, unscaled (control) and with the congruence scaling p.
usage: cont.py inst.dat fixes|- N0 N1 F p[,p2...]"""
import sys, numpy as np, proto
A, B = proto.read_qap(sys.argv[1]); B = (B + B.T) / 2
fx = [] if sys.argv[2] == '-' else [tuple(int(v) for v in s.split(':')) for s in sys.argv[2].split(',')]
N0, N1, F = int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
A2, B2, C2, const = proto.reduce_node(A, B, fx)
r0 = proto.prsm(A2, B2, C2, const, F=F, p=0, iters=N0, every=50, tag='phase0')
for p in [float(x) for x in sys.argv[6].split(',')]:
    proto.prsm(A2, B2, C2, const, F=F, p=p, iters=N1, every=50, start=(r0['Y'], r0['Z']), tag='cont_p%g' % p)
