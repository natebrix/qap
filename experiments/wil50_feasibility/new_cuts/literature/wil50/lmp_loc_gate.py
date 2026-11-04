#!/usr/bin/env python3
"""Location-side 3-local marginal polytope (the half of projected RLT2 that
the step-8 triple trial did not test) on a dumped wil50 iterate.

For three free locations a, b, c the DNN pair blocks Y_(a.),(b.), Y_(a.),(c.),
Y_(b.),(c.) (m x m each, facility x facility) must be the pair marginals of one
distribution over injective facility triples (i, j, k) placed at (a, b, c).
Separation: max sum_pairs <U, block> s.t. U_ab[i,j] + U_ac[i,k] + U_bc[j,k] <= 0
for every injective (i, j, k), -1 <= U <= 1; a positive optimum is a valid cut
<G, Y> <= 0 (every permutation gives one vertex). The LP has 3 m(m-1)
variables and m(m-1)(m-2) rows (97,290 at m = 47); scipy HiGHS.

Triples: the location triples that carry the most objective (sum of
distance-weighted DNN mass on their pairs is not used; we take geometric
families): L-shapes and straight segments of adjacent cells in the grid, plus
random triples, reported separately.

  python3 -I lmp_loc_gate.py PREFIX --fixes SPEC [--ntri 12] [--nrand 6]
"""
import argparse, itertools, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

ap = argparse.ArgumentParser()
ap.add_argument('prefix'); ap.add_argument('--fixes', default='-'); ap.add_argument('--n', type=int, default=50)
ap.add_argument('--grid', default='5:10'); ap.add_argument('--ntri', type=int, default=12)
ap.add_argument('--nrand', type=int, default=6); ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--out', default=None)
a = ap.parse_args()
n = a.n
fixes = [] if a.fixes == '-' else [tuple(map(int, t.split(':'))) for t in a.fixes.split(',')]
free_loc = [i for i in range(n) if i not in {x for x, _ in fixes}]
m = len(free_loc); N = m * m + 1
Y = np.fromfile(a.prefix + '.Y.bin').reshape(N, N)
Y = (Y + Y.T) / 2
rows, cols = map(int, a.grid.split(':'))
loc_of = {loc: k for k, loc in enumerate(free_loc)}
blk = lambda la, lb: Y[1 + la * m:1 + (la + 1) * m, 1 + lb * m:1 + (lb + 1) * m]

# constraint matrix over injective (i, j, k): rows; variables U_ab (m*m), U_ac, U_bc
cfg = np.array([(i, j, k) for i in range(m) for j in range(m) for k in range(m) if i != j and i != k and j != k])
nc = len(cfg)
r = np.repeat(np.arange(nc), 3)
c = np.stack([cfg[:, 0] * m + cfg[:, 1], m * m + cfg[:, 0] * m + cfg[:, 2], 2 * m * m + cfg[:, 1] * m + cfg[:, 2]], 1).reshape(-1)
At = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(nc, 3 * m * m))
# variables with i == j etc. never appear (gangster): fix them to 0 via bounds
diagmask = np.zeros(3 * m * m, dtype=bool)
for t in range(3):
    diagmask[t * m * m + np.arange(m) * m + np.arange(m)] = True
bounds = [(0, 0) if diagmask[v] else (-1, 1) for v in range(3 * m * m)]


def separate(la, lb, lc):
    bvec = np.concatenate([blk(la, lb).reshape(-1), blk(la, lc).reshape(-1), blk(lb, lc).reshape(-1)])
    res = linprog(-bvec, A_ub=At, b_ub=np.zeros(nc), bounds=bounds, method='highs')
    return -res.fun


def cells(trip):
    return [loc_of[x] for x in trip if x in loc_of]

cand = []
for r0 in range(rows):
    for c0 in range(cols):
        x = r0 * cols + c0
        for (d1, d2) in (((0, 1), (0, 2)), ((1, 0), (2, 0)), ((0, 1), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1)), ((1, 0), (1, -1))):
            p1 = (r0 + d1[0], c0 + d1[1]); p2 = (r0 + d2[0], c0 + d2[1])
            if all(0 <= p[0] < rows and 0 <= p[1] < cols for p in (p1, p2)):
                t = (x, p1[0] * cols + p1[1], p2[0] * cols + p2[1])
                if len(cells(t)) == 3:
                    cand.append(t)
rng = np.random.default_rng(a.seed)
# geometric: the central part of the grid (where the DNN's mid regions are)
rng.shuffle(cand)
geo = sorted(cand, key=lambda t: sum(abs(x % cols - 4.5) + abs(x // cols - 2) for x in t))[:a.ntri]
rnd = []
while len(rnd) < a.nrand:
    t = tuple(sorted(rng.choice(free_loc, 3, replace=False).tolist()))
    rnd.append(t)
out = {'geo': [], 'rand': []}
t0 = time.time()
for kind, lst in (('geo', geo), ('rand', rnd)):
    for t in lst:
        la, lb, lc = cells(t)
        v = separate(la, lb, lc)
        out[kind].append((list(map(int, t)), float(v)))
        print(f"{kind:4s} cells {t}  violation {v:.4f}  ({time.time() - t0:.0f}s)", flush=True)
for kind in ('geo', 'rand'):
    vs = [v for _, v in out[kind]]
    print(f"{kind}: max {max(vs):.4f} median {float(np.median(vs)):.4f}")
if a.out:
    json.dump(out, open(a.out, 'w'), indent=1)
print("RESULT", json.dumps({k: max(v for _, v in out[k]) for k in out}))
