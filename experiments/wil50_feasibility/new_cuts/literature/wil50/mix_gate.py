#!/usr/bin/env python3
"""Mixed 3-local polytope gate on a dumped wil50 iterate: the coupled
projection of RLT2 (both kinds of RLT2 equations on the same cubic variables).

Type F (two facilities f1, f2, one location s): configurations
(s1, s2, f3) = locations of f1 and f2 and the facility at s (f3 = f1 if
s1 = s, f2 if s2 = s, else any other facility). The three pair marginals are
DNN entries: Y_(s1 f1),(s2 f2); Y_(s1 f1),(s f3) (x_{s f1} when s1 = s);
Y_(s2 f2),(s f3). Type L swaps the roles (two locations, one facility).
Separation LP as in lmp_loc_gate.py (potentials in [-1, 1]).

  python3 -I mix_gate.py PREFIX --fixes SPEC [--nF 10] [--nL 6] [--nrand 4]
"""
import argparse, itertools, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

ap = argparse.ArgumentParser()
ap.add_argument('prefix'); ap.add_argument('--fixes', default='-'); ap.add_argument('--n', type=int, default=50)
ap.add_argument('--grid', default='5:10'); ap.add_argument('--nF', type=int, default=10)
ap.add_argument('--nL', type=int, default=6); ap.add_argument('--nrand', type=int, default=4)
ap.add_argument('--seed', type=int, default=1); ap.add_argument('--out', default=None)
ap.add_argument('--dat', default='/home/user/qap/qaplib/wil50.dat')
a = ap.parse_args()
n = a.n
fixes = [] if a.fixes == '-' else [tuple(map(int, t.split(':'))) for t in a.fixes.split(',')]
free_loc = [i for i in range(n) if i not in {x for x, _ in fixes}]
free_fac = [j for j in range(n) if j not in {y for _, y in fixes}]
m = len(free_loc); N = m * m + 1
Y = np.fromfile(a.prefix + '.Y.bin').reshape(N, N)
Y = (Y + Y.T) / 2
rows, cols = map(int, a.grid.split(':'))
tok = open(a.dat).read().split(); nn = int(tok[0]); vals = np.array(tok[1:], dtype=float)
Bflow = vals[nn * nn:2 * nn * nn].reshape(nn, nn)[np.ix_(free_fac, free_fac)]
lift = lambda loc, fac: 1 + loc * m + fac


def separate(typ, e1, e2, s):
    """typ F: e1, e2 facilities, s a location; typ L: e1, e2 locations, s a facility"""
    L = (lambda i, j: lift(i, j)) if typ == 'F' else (lambda i, j: lift(j, i))
    var = {}; ent = {}; rows_, cols_ = [], []; c = 0
    def v(k):
        if k not in var: var[k] = len(var)
        return var[k]
    for s1 in range(m):
        for s2 in range(m):
            if s1 == s2: continue
            if s1 == s: f3s = [e1]
            elif s2 == s: f3s = [e2]
            else: f3s = [f for f in range(m) if f != e1 and f != e2]
            for f3 in f3s:
                k12 = ('12', s1, s2); ent[k12] = (L(s1, e1), L(s2, e2))
                k13 = ('13', s1, f3); ent[k13] = (L(s1, e1), L(s, f3))
                k23 = ('23', s2, f3); ent[k23] = (L(s2, e2), L(s, f3))
                for kk in (k12, k13, k23):
                    rows_.append(c); cols_.append(v(kk))
                c += 1
    At = sp.csr_matrix((np.ones(len(rows_)), (rows_, cols_)), shape=(c, len(var)))
    keys = sorted(var, key=var.get)
    bvec = np.array([Y[ent[k]] for k in keys])
    res = linprog(-bvec, A_ub=At, b_ub=np.zeros(c), bounds=[(-1, 1)] * len(var), method='highs')
    return -res.fun

rng = np.random.default_rng(a.seed)
central = [k for k, loc in enumerate(free_loc) if 1 <= loc // cols <= 3 and 3 <= loc % cols <= 6]
iu = np.triu_indices(m, 1)
heavy = np.argsort(-Bflow[iu])[:a.nF]
tests = []
for h in heavy:
    tests.append(('F', int(iu[0][h]), int(iu[1][h]), int(rng.choice(central)), 'heavy'))
for _ in range(a.nrand):
    f1, f2 = rng.choice(m, 2, replace=False)
    tests.append(('F', int(f1), int(f2), int(rng.integers(m)), 'rand'))
for _ in range(a.nL):
    s1, s2 = rng.choice(central, 2, replace=False)
    tests.append(('L', int(s1), int(s2), int(rng.integers(m)), 'central'))
out = []
t0 = time.time()
for (typ, e1, e2, s, kind) in tests:
    val = separate(typ, e1, e2, s)
    out.append((typ, e1, e2, s, kind, float(val)))
    print(f"{typ} {kind:7s} ({e1},{e2};{s})  violation {val:.4f}  ({time.time() - t0:.0f}s)", flush=True)
for typ in ('F', 'L'):
    for kind in ('heavy', 'rand', 'central'):
        vs = [o[5] for o in out if o[0] == typ and o[4] == kind]
        if vs:
            print(f"{typ} {kind}: max {max(vs):.4f} median {float(np.median(vs)):.4f}")
if a.out:
    json.dump(out, open(a.out, 'w'), indent=1)
print("RESULT max %.4f" % max(o[5] for o in out))
