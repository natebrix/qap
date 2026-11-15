#!/usr/bin/env python3
"""Long-baseline dual extrapolation, certificate side only: Z(t) = Z_k2 + t (Z_k2 - Z_k1) from two stored
certificate multipliers of one PRSM run, certified like any Z (positive part of V'ZV projected out, then
Gilmore-Lawler). Is the certified value at the best t ahead of the trajectory's own value at k2?
usage: extrap.py inst.dat fixes|- p iters k1,k2,k3..."""
import sys, numpy as np, proto
A, B = proto.read_qap(sys.argv[1]); B = (B + B.T) / 2
fx = [] if sys.argv[2] == '-' else [tuple(int(v) for v in s.split(':')) for s in sys.argv[2].split(',')]
p = float(sys.argv[3]); iters = int(sys.argv[4]); ks = [int(x) for x in sys.argv[5].split(',')]
A2, B2, C2, const = proto.reduce_node(A, B, fx)
m = A2.shape[0]
L, free, V = proto.build(A2, B2, C2); scale = max(1.0, np.abs(L).max()); Ls = L / scale
snap = {k: None for k in ks}
r = proto.prsm(A2, B2, C2, const, F=0.21, p=p, iters=iters, every=50, snap=snap, tag='run')
hist = {h[0]: h[2] for h in r['hist']}
for a in range(len(ks)):
    for b in range(a + 1, len(ks)):
        k1, k2 = ks[a], ks[b]
        Z1, Z2 = snap[k1], snap[k2]
        row = []
        for t in [0, 0.25, 0.5, 1, 2, 4]:
            lb = proto.certify(Ls, V, Z2 + t * (Z2 - Z1), m) * scale + const
            row.append('t=%g: %.0f' % (t, lb))
        print('EXTRAP k1 %d k2 %d | own cert at k2 %.0f |' % (k1, k2, hist.get(k2, float('nan'))), ' | '.join(row), flush=True)
