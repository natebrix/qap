#!/usr/bin/env python3
"""Local marginal polytope of seven regional variables: the memo's cheap gate for
the regional moment block (candidate 6) and Astra's moment separation (10b).

For a region S of free locations and seven free facilities T, the regional
moments p_i = sum_{a in S} Y_(a,i),(a,i) and W_ij = sum_{a,b in S} Y_(a,i),(b,j)
are the first and second moments of the 0/1 vector z_T = (facility i placed in
S), so (p, W) must lie in the convex hull of the admissible configurations
c in {0,1}^7 (those with sum c <= |S| and 7 - sum c <= m - |S|). The most
violated separating inequality  a.p + <B, W> <= beta  with |a_i|, |B_ij| <= 1
(so violations are on the triangle inequalities' scale) is an LP with 29
variables and at most 129 constraints; the triangles are the case |T| = 3.
Every inequality found is valid for every permutation, so a cut written out
(--dump-top) can go through dnntest/dnncuts -Ladd at a fixed multiplier via
make_ladd7.py exactly as the triple cuts did.

  python3 local7_separate.py PREFIX --fixes SPEC [--n 50] [--grid 5:10]
      [--subsets 300] [--seed 1] [--dump-top 0] [--out PREFIX.local7.json]
Reads PREFIX.Y.bin (raw doubles, (m^2+1)^2 row-major, dnntest/dnncuts -dumpY).
"""
import argparse, itertools, json, random, time
import numpy as np
from scipy.optimize import linprog

ap = argparse.ArgumentParser()
ap.add_argument('prefix'); ap.add_argument('--fixes', default='-'); ap.add_argument('--n', type=int, default=50)
ap.add_argument('--grid', default='5:10', help='rows:cols of the location grid (first matrix)')
ap.add_argument('--subsets', type=int, default=300, help='random 7-subsets per region')
ap.add_argument('--seed', type=int, default=1); ap.add_argument('--dump-top', type=int, default=0)
ap.add_argument('--out', default=None)
a = ap.parse_args()
n = a.n
fixes = [] if a.fixes == '-' else [tuple(map(int, t.split(':'))) for t in a.fixes.split(',')]
free_loc = [i for i in range(n) if i not in {x for x, _ in fixes}]      # first matrix (A): locations
free_fac = [j for j in range(n) if j not in {y for _, y in fixes}]      # second matrix (B): facilities
m = len(free_loc); N = m * m + 1
Y = np.fromfile(a.prefix + '.Y.bin')
assert Y.size == N * N, (Y.size, N)
Y = Y.reshape(N, N)
rows, cols = map(int, a.grid.split(':'))
# regions: row and column prefixes of the grid over the original location indices, in local indices
regions = []
for t in range(1, rows):
    regions.append((f'row<{t}', [k for k, loc in enumerate(free_loc) if loc // cols < t]))
for t in range(1, cols):
    regions.append((f'col<{t}', [k for k, loc in enumerate(free_loc) if loc % cols < t]))
regions = [(nm, S) for nm, S in regions if 0 < len(S) < m]
idx = lambda loc, fac: 1 + loc * m + fac

def moments(S):
    """p (m), W (m x m) of region S from Y"""
    I = np.array([idx(loc, f) for loc in S for f in range(m)]).reshape(len(S), m)
    p = np.zeros(m); W = np.zeros((m, m))
    for r in range(len(S)):
        p += Y[I[r], I[r]]
        for s in range(len(S)):
            W += Y[np.ix_(I[r], I[s])]
    return p, W

K = 7
configs = np.array(list(itertools.product([0, 1], repeat=K)))            # 128 x 7
pairs = [(i, j) for i in range(K) for j in range(i + 1, K)]            # 21
# LP: variables x = (a_1..a_7, B_12..B_67, beta); maximize a.p + B.W - beta
# s.t. a.c + B.(c c^T) - beta <= 0 for admissible c; |a|,|B| <= 1; beta free
def separate(p7, W7, ksz):
    adm = configs[(configs.sum(1) <= ksz) & (K - configs.sum(1) <= m - ksz)]
    A = np.hstack([adm, np.array([[c[i] * c[j] for i, j in pairs] for c in adm]), -np.ones((len(adm), 1))])
    obj = -np.concatenate([p7, [W7[i, j] for i, j in pairs], [-1.0]])
    bounds = [(-1, 1)] * (K + len(pairs)) + [(None, None)]
    lp = linprog(obj, A_ub=A, b_ub=np.zeros(len(adm)), bounds=bounds, method='highs')
    if lp.x is None: return None
    return -lp.fun, lp.x

def triangle_max(p, W, T):
    v = 0.0
    for i, j, k in itertools.combinations(T, 3):
        v = max(v, W[i, j] + W[i, k] - W[j, k] - p[i], W[i, j] + W[j, k] - W[i, k] - p[j],
                W[i, k] + W[j, k] - W[i, j] - p[k], p[i] + p[j] + p[k] - W[i, j] - W[i, k] - W[j, k] - 1)
    return v

rng = random.Random(a.seed)
summary = {'prefix': a.prefix, 'fixes': a.fixes, 'm': m, 'regions': {}, 'subsets_per_region': a.subsets}
top = []
t0 = time.time()
for name, S in regions:
    p, W = moments(S)
    frac = sorted(range(m), key=lambda i: -p[i] * (1 - p[i]))
    subsets = [tuple(sorted(frac[:K]))] + [tuple(sorted(rng.sample(range(m), K))) for _ in range(a.subsets)]
    viol = []; tri = []
    for T in subsets:
        p7 = p[list(T)]; W7 = W[np.ix_(T, T)]
        r = separate(p7, W7, len(S))
        if r is None: continue
        v, x = r
        viol.append(v); tri.append(triangle_max(p, W, T))
        if a.dump_top > 0:
            top.append((v, name, S, T, x))
    viol = np.array(viol); tri = np.array(tri)
    summary['regions'][name] = dict(size=len(S), subsets=len(viol), max=float(viol.max()), mean=float(viol.mean()),
        over_0_05=int((viol > 0.05).sum()), over_0_1=int((viol > 0.1).sum()), over_0_2=int((viol > 0.2).sum()),
        triangle_max=float(tri.max()), most_fractional_subset=dict(T=list(subsets[0]), violation=float(viol[0]), triangle=float(tri[0])),
        excess_over_triangles_max=float((viol - tri).max()))
    print(f"{name:7s} |S|={len(S):2d} subsets {len(viol):4d}  local7 max {viol.max():.4f} mean {viol.mean():.4f}  "
          f">0.05: {(viol>0.05).sum():4d}  >0.1: {(viol>0.1).sum():4d}  triangle max {tri.max():.4f}  "
          f"excess max {(viol-tri).max():.4f}  mostfrac {viol[0]:.4f}", flush=True)
allv = [r['max'] for r in summary['regions'].values()]
summary['overall'] = dict(max=max(allv), seconds=time.time() - t0)
print(f"overall max {max(allv):.4f} in {time.time()-t0:.0f} s")
# dump the top cuts as dense G (sum of <G,Y> <= beta) for make_ladd7.py
if a.dump_top > 0:
    top.sort(key=lambda t: -t[0])
    dumped = []
    for rank, (v, name, S, T, x) in enumerate(top[:a.dump_top]):
        acoef = x[:K]; bcoef = x[K:K + len(pairs)]; beta = float(x[-1])
        G = np.zeros((N, N))
        for i in range(K):
            for loc in S: G[idx(loc, T[i]), idx(loc, T[i])] += acoef[i]
        for q, (i, j) in enumerate(pairs):
            for loc in S:
                for loc2 in S:
                    r, c = idx(loc, T[i]), idx(loc2, T[j])
                    G[r, c] += bcoef[q] / 2; G[c, r] += bcoef[q] / 2
        chk = float((G * Y).sum() - beta)
        path = f"{a.prefix}.local7_{rank}"
        G.tofile(path + '.G.bin')
        dumped.append(dict(rank=rank, region=name, T=list(T), violation=float(v), check=chk, beta=beta, file=path + '.G.bin'))
    summary['dumped'] = dumped
    print(f"dumped {len(dumped)} cuts, violations {dumped[-1]['violation']:.4f}..{dumped[0]['violation']:.4f}")
out = a.out or a.prefix + '.local7.json'
json.dump(summary, open(out, 'w'), indent=1)
print("summary ->", out)
