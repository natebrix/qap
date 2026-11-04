#!/usr/bin/env python3
"""Facility-side distance cuts on a dumped iterate. The objective sees Y only
through the expected distances Delta_jl = sum_ab A_ab Y_(aj),(bl) between free
facilities j, l (and the linear x terms), so valid inequalities on Delta are
the facility-side inequalities that matter for it. For a set T of facilities
and weights H on its pairs, sum_{j<l in T} H_jl d(pi_j, pi_l) >= h(H) :=
min over injective placements of T into the node's free locations; in Y,
  - sum_{j<l} H_jl sum_{a != b} A_ab Y_(aj),(bl) <= -h(H).
Families:
  parity  |T| = 3, H = 1, h = min sum of three pairwise distances (4 on the
          grid: it is bipartite, so the three distances have even sum, each >= 1)
  hull3   |T| = 3, H from an LP that separates (Delta_12, Delta_13, Delta_23)
          from the convex hull of the achievable distance triples
  sub4    |T| = 4, H = B_TT (the flows), h by enumeration of all placements
  python3 sep_dist.py SPEC Y.bin --fam parity|hull3|sub4 [--top K] [--out F.cb]"""
import argparse, itertools, sys, time
import numpy as np
from ynode import Node
import cutio

def lifted_pairs(nd, pairs_w):
    """entries of - sum_(j,l) w_jl Delta_jl as (r, c, g)"""
    m = nd.m; A = nd.A
    a, b = np.nonzero(A)
    R = []; C = []; G = []
    for (j, l), w in pairs_w:
        R.append(1 + a * m + j); C.append(1 + b * m + l); G.append(-w * A[a, b])
    return np.concatenate(R), np.concatenate(C), np.concatenate(G)

def achievable_triples(nd):
    A = nd.A.astype(int); m = nd.m
    D = set()
    for a in range(m):
        for b in range(m):
            if b == a: continue
            dab = A[a, b]
            for c in range(m):
                if c == a or c == b: continue
                D.add((dab, A[a, c], A[b, c]))
    return np.array(sorted(D), dtype=float)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--fam', default='parity'); ap.add_argument('--top', type=int, default=500)
    ap.add_argument('--out'); ap.add_argument('--min-viol', type=float, default=1e-3)
    ap.add_argument('--check', type=int, default=200, help='random permutations for the validity check')
    a = ap.parse_args()
    t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m
    Dl = nd.delta(); B = nd.B
    cuts = []; found = []
    if a.fam in ('parity', 'hull3'):
        Dset = achievable_triples(nd)
        h3 = Dset.sum(1).min()
        trip = np.array(list(itertools.combinations(range(m), 3)))
        v12 = Dl[trip[:, 0], trip[:, 1]]; v13 = Dl[trip[:, 0], trip[:, 2]]; v23 = Dl[trip[:, 1], trip[:, 2]]
        s3 = v12 + v13 + v23
        print('achievable distance triples: %d, min sum %g; DNN: min sum %.4f, %d triples below %g' % (
            len(Dset), h3, s3.min(), (s3 < h3 - a.min_viol).sum(), h3))
        if a.fam == 'parity':
            order = np.argsort(s3)
            for t in order[:a.top]:
                v = h3 - s3[t]
                if v < a.min_viol: break
                j, k, l = trip[t]
                found.append((v, (int(j), int(k), int(l)), [((j, k), 1.0), ((j, l), 1.0), ((k, l), 1.0)], h3))
        else:
            from scipy.optimize import linprog
            # LP: max t  s.t.  t <= H.(d - Delta) for all d in D,  -1 <= H <= 1
            order = np.argsort(s3)[:max(a.top * 4, 2000)]
            for t in order:
                dv = np.array([v12[t], v13[t], v23[t]])
                M = Dset - dv
                # variables H1 H2 H3 tt; minimize -tt; constraints tt - M H <= 0
                Aub = np.hstack([-M, np.ones((len(M), 1))])
                res = linprog([0, 0, 0, -1], A_ub=Aub, b_ub=np.zeros(len(M)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs')
                if res.status != 0: continue
                H = res.x[:3]; viol = -res.fun
                if viol < a.min_viol: continue
                H = np.round(H * 64) / 64          # rational coefficients, exact right-hand side below
                hH = (Dset @ H).min(); viol = hH - H @ dv
                if viol < a.min_viol: continue
                j, k, l = trip[t]
                found.append((viol, (int(j), int(k), int(l)), [((j, k), H[0]), ((j, l), H[1]), ((k, l), H[2])], hH))
    elif a.fam == 'sub4':
        locs = range(m); A = nd.A
        # all ordered 4-placements: distance 6-vectors (d01 d02 d03 d12 d13 d23)
        P = np.array([p for p in itertools.permutations(range(m), 4)], dtype=np.int16) if m <= 30 else None
        if P is None:
            # m large: sample the distinct distance vectors by placements within a radius
            rng = np.random.default_rng(0)
            Pl = []
            for a0 in range(m):
                near = [b for b in range(m) if b != a0 and A[a0, b] <= 4]
                for rest in itertools.permutations(near, 3):
                    Pl.append((a0,) + rest)
            P = np.array(Pl, dtype=np.int32)
        Dv = np.stack([A[P[:, 0], P[:, 1]], A[P[:, 0], P[:, 2]], A[P[:, 0], P[:, 3]], A[P[:, 1], P[:, 2]], A[P[:, 1], P[:, 3]], A[P[:, 2], P[:, 3]]], 1)
        Dv = np.unique(Dv, axis=0).astype(float)
        print('distinct 6-distance vectors (placements within radius 4 of the first): %d' % len(Dv))
        # candidate 4-subsets: heaviest by flow weight sum
        quads = []
        for T in itertools.combinations(range(m), 4):
            pass
        rng = np.random.default_rng(1)
        cand = set()
        # greedy: each heavy triangle extended by the best fourth facility
        trip = np.array(list(itertools.combinations(range(m), 3)))
        wt = B[trip[:, 0], trip[:, 1]] + B[trip[:, 0], trip[:, 2]] + B[trip[:, 1], trip[:, 2]]
        for t in np.argsort(-wt)[:3000]:
            j, k, l = trip[t]
            for q in range(m):
                if q in (j, k, l): continue
                cand.add(tuple(sorted((j, k, l, q))))
        cand = sorted(cand)
        print('candidate 4-subsets: %d' % len(cand))
        for T in cand:
            pairs = list(itertools.combinations(T, 2))
            H = np.array([B[p] for p in pairs]); dv = np.array([Dl[p] for p in pairs])
            if H.sum() == 0: continue
            hH = (Dv @ H).min()       # over the sampled placements: an upper estimate of h, used only to screen
            viol = hH - H @ dv
            if viol > a.min_viol:
                found.append((viol, T, list(zip(pairs, H)), hH))
        found.sort(key=lambda z: -z[0]); found = found[:a.top]
        # exact right-hand sides for the kept ones: full enumeration over all placements
        exact = []
        for viol, T, pw, hH in found:
            H = np.array([w for p, w in pw])
            best = np.inf
            for a0 in range(m):
                for a1 in range(m):
                    if a1 == a0: continue
                    for a2 in range(m):
                        if a2 in (a0, a1): continue
                        a3 = np.array([b for b in range(m) if b not in (a0, a1, a2)])
                        v = H[0] * A[a0, a1] + H[1] * A[a0, a2] + H[3] * A[a1, a2] + H[2] * A[a0, a3] + H[4] * A[a1, a3] + H[5] * A[a2, a3]
                        best = min(best, v.min())
            dv = np.array([Dl[p] for p, w in pw])
            exact.append((best - H @ dv, T, pw, best))
        found = [e for e in exact if e[0] > a.min_viol]
    found.sort(key=lambda z: -z[0])
    found = found[:a.top]
    if found:
        vs = np.array([f[0] for f in found])
        print('%s: %d cuts kept, violation %.4f..%.4f (sum %.3f), weighted by sum |H|: max %.4f' % (
            a.fam, len(found), vs.min(), vs.max(), vs.sum(), max(f[0] / sum(abs(w) for p, w in f[2]) for f in found)))
        for f in found[:8]:
            print('   viol %.4f  T %s  H %s  h %.4f' % (f[0], f[1], [round(float(w), 3) for p, w in f[2]], f[3]))
    else:
        print('%s: nothing violated by more than %g' % (a.fam, a.min_viol))
    if a.out and found:
        for k, (v, T, pw, h) in enumerate(found):
            r, c, g = lifted_pairs(nd, pw)
            key = hash((a.fam,) + tuple(T) + tuple(round(float(w), 6) for p, w in pw)) & ((1 << 62) - 1)
            cuts.append((key, -h, r, c, g))
        cutio.write_cb(a.out, cuts)
        # check: the separated violation reproduces on Y, and no permutation violates
        rv = [cutio.value(cu, nd.Y) for cu in cuts[:20]]
        print('wrote %d cuts to %s; values on Y (first 20) match: %s; perm check over %d permutations: max %.3g (must be <= 0)' % (
            len(cuts), a.out, np.allclose(rv, [f[0] for f in found[:20]], atol=1e-6), a.check, cutio.perm_check(cuts[:50], m, a.check)))
    print('time %.1f s' % (time.time() - t0))

if __name__ == '__main__':
    main()
