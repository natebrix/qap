#!/usr/bin/env python3
"""Lagrangian terms for one side of a count split (disjunction c <= beta OR
c >= beta + 1) for fdrv -cont -extra FILE:CONST. c = sum_{j in Q+} z_j +
sum_{j in Q-} (1 - z_j), z_j = sum_{a in S} x_aj (S a location region, Q
signed facilities, local indices). On the side 'le' every feasible Y has
c - beta <= 0, so min <L + mu (c - beta), Y> is a lower bound on that side's
value for any mu >= 0 (likewise 'ge' with beta + 1 - c). The term is linear:
entries (0, (aj)) and ((aj), 0) get half the coefficient each.
  python3 count_term.py SPEC --region col<5 --plus 1,2 --minus 7 --beta 2 --side le --mu 20 --out F"""
import argparse
import numpy as np
from ynode import Node
ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('--region'); ap.add_argument('--pick', help='pick_split.py log: take region, Q+, Q-, beta from its ARGS line')
ap.add_argument('--plus', default=''); ap.add_argument('--minus', default='')
ap.add_argument('--beta', type=int); ap.add_argument('--side', choices=['le', 'ge'], required=True)
ap.add_argument('--mu', type=float, required=True); ap.add_argument('--out', required=True)
a = ap.parse_args()
if a.pick:
    toks = [l for l in open(a.pick) if l.startswith('ARGS')][0].split()[1:]
    kv = dict(zip(toks[0::2], toks[1::2]))
    a.region = kv['--region']; a.plus = kv.get('--plus', ''); a.minus = kv.get('--minus', ''); a.beta = int(kv['--beta'])
nd = Node(a.spec); m = nd.m; n1 = m * m + 1
kind, t = a.region.split('<'); t = int(t)
inS = (nd.rc[:, 0] < t) if kind == 'row' else (nd.rc[:, 1] < t)
plus = [int(v) for v in a.plus.split(',') if v]; minus = [int(v) for v in a.minus.split(',') if v]
s = np.zeros(m); s[plus] = 1; s[minus] = -1
sign = 1.0 if a.side == 'le' else -1.0
L = np.zeros((n1, n1))
for aa in np.where(inS)[0]:
    for j in np.where(s != 0)[0]:
        r = 1 + aa * m + j
        L[0, r] += sign * a.mu * s[j] / 2; L[r, 0] += sign * a.mu * s[j] / 2
nm = len(minus)
# objective + term = objective + sign*mu*(s.z) + mu*(nm - beta) [le]  or  + mu*(beta + 1 - nm) [ge]
add = a.mu * (nm - a.beta) if a.side == 'le' else a.mu * (a.beta + 1 - nm)
const = -add     # fdrv: shift - const
L.tofile(a.out)
print('%s:%.12g' % (a.out, const))
