#!/usr/bin/env python3
"""Facility-side diagnostics of a dumped DNN iterate of a wil50 node.
  python3 ydiag.py SPEC Y.bin [--primal P]
Everything is computed from Y (ynode.py conventions). Violations are of valid
inequalities (each proved in ../NOTES.md); positive = violated."""
import argparse, itertools, sys
import numpy as np
from ynode import Node

def tri_viol(p, W):
    """max violation of the 4 Boolean triangle types over all triples of a 0/1
    vector's moments (p, W), and the count above 1e-3"""
    m = len(p); best = -1e9; cnt = 0
    Wd = W.copy(); np.fill_diagonal(Wd, p)
    for i in range(m):
        # apex i: W_ij + W_ik - W_jk - p_i
        v = W[i][:, None] + W[i][None, :] - W - p[i]
        v[i, :] = -9; v[:, i] = -9; np.fill_diagonal(v, -9)
        best = max(best, v.max()); cnt += int((v > 1e-3).sum() // 2)
    iu = list(itertools.combinations(range(m), 3)) if m <= 60 else []
    a = np.array(iu)
    s = p[a[:, 0]] + p[a[:, 1]] + p[a[:, 2]] - W[a[:, 0], a[:, 1]] - W[a[:, 0], a[:, 2]] - W[a[:, 1], a[:, 2]] - 1
    best = max(best, s.max()); cnt += int((s > 1e-3).sum())
    return best, cnt

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y'); ap.add_argument('--primal', type=float)
    a = ap.parse_args()
    nd = Node(a.spec).load(a.Y); m = nd.m; T = nd.T; x = nd.x
    D = nd.delta()
    print('node %s m=%d  objective from Y %.4f%s' % (a.spec, m, nd.objective(), '' if a.primal is None else ' (log primal %.4f)' % a.primal))
    print('x: max %.3f, entries > 0.5: %d, row max median %.3f' % (x.max(), (x > 0.5).sum(), np.median(x.max(0))))
    iu = np.triu_indices(m, 1)
    print('Delta: min %.4f mean %.3f (random placement mean %.3f)' % (D[iu].min(), D[iu].mean(), nd.A[np.triu_indices(m, 1)].mean()))
    # 1. parity/dispersion on triples: Delta_jk + Delta_jl + Delta_kl >= 4
    trip = np.array(list(itertools.combinations(range(m), 3)))
    s3 = D[trip[:, 0], trip[:, 1]] + D[trip[:, 0], trip[:, 2]] + D[trip[:, 1], trip[:, 2]]
    wB = nd.B[trip[:, 0], trip[:, 1]] + nd.B[trip[:, 0], trip[:, 2]] + nd.B[trip[:, 1], trip[:, 2]]
    print('triples: min sum of expected distances %.4f (valid >= 4); %d of %d below 4; below 4.05: %d; heavy (flow sum >= 24) min %.3f' % (
        s3.min(), (s3 < 4 - 1e-6).sum(), len(s3), (s3 < 4.05).sum(), s3[wB >= 24].min()))
    # P(distance = 1) for pairs, and P(odd)
    adj = (nd.A == 1).astype(float); odd = (nd.A % 2 == 1).astype(float)
    P1 = np.einsum('ab,ajbl->jl', adj, T, optimize=True); Po = np.einsum('ab,ajbl->jl', odd, T, optimize=True)
    print('pairs: max P(adjacent) %.3f; heavy pairs (B=9) mean P(adj) %.3f, mean Delta %.3f; B=0 pairs mean Delta %.3f' % (
        P1[iu].max(), P1[iu][nd.B[iu] == 9].mean(), D[iu][nd.B[iu] == 9].mean(), D[iu][nd.B[iu] == 0].mean()))
    # 2. checkerboard region (black = (r + c) even): regional triangles
    col = (nd.rc.sum(1) % 2 == 0)
    def region_moments(mask):
        idx = np.where(mask)[0]
        p = x[idx].sum(0)
        W = T[np.ix_(idx, range(m), idx, range(m))].sum(axis=(0, 2))
        return p, W
    pc, Wc = region_moments(col)
    v, c = tri_viol(pc, Wc)
    so = pc[:, None] + pc[None, :] - 2 * Wc   # P(different colours) = P(odd distance)
    print('checkerboard region (%d free cells): max triangle violation %.4f, %d violated > 1e-3 (P_odd check: %.2e)' % (
        col.sum(), v, c, np.abs(so - Po)[iu].max()))
    # 3. the 13 threshold regions, for scale
    out = []
    for t in range(1, 5):
        out.append(('row<%d' % t, nd.rc[:, 0] < t))
    for t in range(1, 10):
        out.append(('col<%d' % t, nd.rc[:, 1] < t))
    tv = []
    for name, mask in out:
        p, W = region_moments(mask); v, c = tri_viol(p, W); tv.append((name, int(mask.sum()), v, c))
    print('threshold regions (name, capacity, max violation, count):', ' '.join('%s:%d:%.4f:%d' % t for t in tv))
    # 4. facility-side regions: location triples, facility sets
    rng = np.random.default_rng(0)
    def fac_region(Tset):
        Tset = np.array(sorted(Tset))
        q = x[:, Tset].sum(1)
        U = T[np.ix_(range(m), Tset, range(m), Tset)].sum(axis=(1, 3))
        return q, U
    res = []
    for k in [k for k in (10, 20, 23) if k < m]:
        for rep in range(3):
            S = rng.choice(m, k, replace=False)
            q, U = fac_region(S); v, c = tri_viol(q, U); res.append((k, v, c))
    print('facility regions (random sets; size, max violation, count):', ' '.join('%d:%.4f:%d' % r for r in res))
    np.save(a.Y.replace('.Y.bin', '') + '.Delta.npy', D)

if __name__ == '__main__':
    main()
