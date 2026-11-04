#!/usr/bin/env python3
"""LOCATION-triple version of lmp3_sep.py (LMP3_loc): compact location triples
(perimeter 4, L-shapes) instead of heavy facility triples; same LP on the
transposed blocks (pairs of facilities at the three locations), same exact
right-hand side.  Original docstring:
Separate facility-triple consistency cuts (LMP3_fac) for the heaviest-flow
triples of a node: the dual of lmp3.py's min-L1 LP gives potentials pi on the
three pair blocks with sum over the pair entries of any injective placement
<= 0; pi is rounded to multiples of 1/1024 and the right-hand side recomputed
exactly as the max over all placements, so each cut <pi, blocks(Y)> <= h is
valid for every permutation independently of the LP. Writes a .cb file.
  python3 lmp3_sep.py SPEC Y.bin --k 40 --out F.cb"""
import argparse, itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from ynode import Node
import cutio

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y'); ap.add_argument('--k', type=int, default=40)
    ap.add_argument('--out', required=True); ap.add_argument('--seconds', type=float, default=1500)
    a = ap.parse_args(); t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m; B = nd.A.max() - nd.A; T = nd.T.transpose(1, 0, 3, 2)
    trip = np.array(list(itertools.combinations(range(m), 3)))
    w = B[trip[:, 0], trip[:, 1]] + B[trip[:, 0], trip[:, 2]] + B[trip[:, 1], trip[:, 2]]
    order = np.random.default_rng(0).permutation(np.where(w == w.max())[0])[:a.k]
    idx = np.array([(x, y, z) for x in range(m) for y in range(m) for z in range(m) if x != y and x != z and y != z])
    pairs = [(x, y) for x in range(m) for y in range(m) if x != y]; pid = {q: t for t, q in enumerate(pairs)}; npr = len(pairs)
    P1 = np.array([pid[(x, y)] for x, y, z in idx]); P2 = np.array([pid[(x, z)] for x, y, z in idx]); P3 = np.array([pid[(y, z)] for x, y, z in idx])
    nz = len(idx)
    rows = np.concatenate([P1, npr + P2, 2 * npr + P3]); cols = np.tile(np.arange(nz), 3)
    M = sp.csr_matrix((np.ones(3 * nz), (rows, cols)), shape=(3 * npr, nz))
    I = sp.identity(3 * npr, format='csr'); Aeq = sp.hstack([M, I, -I], format='csr')
    c = np.concatenate([np.zeros(nz), np.ones(6 * npr)])
    PA = np.array(pairs)
    cuts = []; viols = []
    for t in order:
        if time.time() - t0 > a.seconds: break
        i, k, p = trip[t]
        rhs = np.concatenate([T[PA[:, 0], i, PA[:, 1], k], T[PA[:, 0], i, PA[:, 1], p], T[PA[:, 0], k, PA[:, 1], p]])
        res = linprog(c, A_eq=Aeq, b_eq=rhs, bounds=(0, None), method='highs')
        if res.status != 0 or res.fun < 1e-4: continue
        pi = np.round(np.asarray(res.eqlin.marginals) * 1024) / 1024     # dual potentials, |pi| <= 1
        # exact right-hand side: max over placements of pi1(x,y) + pi2(x,z) + pi3(y,z)
        h = float((pi[P1] + pi[npr + P2] + pi[2 * npr + P3]).max())
        v = float(pi @ rhs - h)
        if v < 1e-4: continue
        r = np.concatenate([1 + i * m + PA[:, 0], 1 + i * m + PA[:, 0], 1 + k * m + PA[:, 0]])
        cc = np.concatenate([1 + k * m + PA[:, 1], 1 + p * m + PA[:, 1], 1 + p * m + PA[:, 1]])
        g = pi.copy(); keep = g != 0
        cuts.append((int(i) * 10000 + int(k) * 100 + int(p), h, r[keep], cc[keep], g[keep]))
        viols.append((v, nd.loc[i], nd.loc[k], nd.loc[p], int(w[t])))
        print('location triple (%2d,%2d,%2d) score %d: LP %.4f  rounded cut violation %.4f  h %.4f  nnz %d  (%.0f s)' % (
            nd.loc[i], nd.loc[k], nd.loc[p], w[t], res.fun, v, h, keep.sum(), time.time() - t0), flush=True)
    cutio.write_cb(a.out, cuts)
    vs = np.array([v[0] for v in viols]) if viols else np.zeros(1)
    chk = cutio.perm_check(cuts[:10], m, 300) if cuts else 0
    print('wrote %d cuts to %s: violation sum %.4f max %.4f; values on Y reproduce: %s; perm check (10 cuts x 300 perms) max %.3g (<= 0)' % (
        len(cuts), a.out, vs.sum(), vs.max(), np.allclose([cutio.value(cu, nd.Y) for cu in cuts], [v[0] for v in viols], atol=1e-6) if cuts else True, chk))

if __name__ == '__main__':
    main()
