#!/usr/bin/env python3
"""Facility-triple consistency (the 3-local marginal polytope on the facility
side, LMP3_fac; RLT2's projection on a triple): is there a distribution z over
injective placements (a, b, c) of facilities (i, k, p) whose pair marginals are
the DNN's blocks Y_(ai),(bk), Y_(ai),(cp), Y_(bk),(cp)? LP: minimise the L1
deviation  sum |marginal(z) - Y-block|  (0 = consistent). Compares the heaviest
flow triples with random ones.  python3 lmp3.py SPEC Y.bin [--k 5]"""
import argparse, itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from ynode import Node

def inconsistency(T, i, k, p, m):
    idx = [(a, b, c) for a in range(m) for b in range(m) for c in range(m) if a != b and a != c and b != c]
    nz = len(idx); A_ = np.array(idx)
    pairs = [(a, b) for a in range(m) for b in range(m) if a != b]; pid = {q: t for t, q in enumerate(pairs)}; npr = len(pairs)
    rows = np.concatenate([[pid[(a, b)] for a, b, c in idx], [npr + pid[(a, c)] for a, b, c in idx], [2 * npr + pid[(b, c)] for a, b, c in idx]])
    cols = np.tile(np.arange(nz), 3)
    M = sp.csr_matrix((np.ones(3 * nz), (rows, cols)), shape=(3 * npr, nz))
    rhs = np.concatenate([[T[a, i, b, k] for a, b in pairs], [T[a, i, c, p] for a, c in pairs], [T[b, k, c, p] for b, c in pairs]])
    # M z + s+ - s- = rhs, minimise sum s+ + s-
    I = sp.identity(3 * npr, format='csr')
    Aeq = sp.hstack([M, I, -I], format='csr')
    c = np.concatenate([np.zeros(nz), np.ones(6 * npr)])
    res = linprog(c, A_eq=Aeq, b_eq=rhs, bounds=(0, None), method='highs')
    return res.fun if res.status == 0 else np.nan

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y'); ap.add_argument('--k', type=int, default=4); ap.add_argument('--loc', action='store_true', help='location triples (LMP3_loc): roles swapped, compact vs random triples')
    a = ap.parse_args(); t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m; B = nd.B
    Tt = nd.T.transpose(1, 0, 3, 2) if a.loc else nd.T     # [j, a, l, b]: location triples as 'facilities'
    W = (nd.A.max() - nd.A) if a.loc else B                 # compact location triples score high
    trip = np.array(list(itertools.combinations(range(m), 3)))
    w = W[trip[:, 0], trip[:, 1]] + W[trip[:, 0], trip[:, 2]] + W[trip[:, 1], trip[:, 2]]
    rng = np.random.default_rng(0)
    heavy = trip[rng.permutation(np.where(w == w.max())[0])[:a.k]] if a.loc else trip[np.argsort(-w)[:a.k]]
    rnd = trip[rng.choice(len(trip), a.k, replace=False)]
    for name, S in ((('compact-locations' if a.loc else 'heaviest-flow'), heavy), ('random', rnd)):
        for (i, k, p) in S:
            v = inconsistency(Tt, i, k, p, m)
            lab = (nd.loc if a.loc else nd.fac)
            print('%-17s triple (%2d,%2d,%2d) %s %2d: L1 inconsistency of the three pair blocks %.4f (%.0f s)' % (
                name, lab[i], lab[k], lab[p], 'perimeter' if a.loc else 'flow sum', (nd.A if a.loc else B)[i, k] + (nd.A if a.loc else B)[i, p] + (nd.A if a.loc else B)[k, p], v, time.time() - t0), flush=True)

if __name__ == '__main__':
    main()
