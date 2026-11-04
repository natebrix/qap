#!/usr/bin/env python3
"""Cutting-plane estimates of DNN + family on a small instance (nug12 root):
solve the DNN (SCS) with the current cut pool, separate the K most violated
cuts of the family at the primal Y, repeat. Each value is the optimum of the
DNN plus finitely many valid cuts, so it is a lower estimate of DNN + family
and an upper bound... no: it is a relaxation value, valid as a lower bound on
the QAP, and at most the value of DNN + the whole family.

Families:
  box1   Jünger-Kaibel 1-box inequalities: for P (side-0 set) and Q (side-1
         set), q = #{j in Q : pi^-1(j) in P} is an integer, so
         E[(q - b)(q - b - 1)] >= 0, i.e.
         2b sum_{a in P, j in Q} Y_(aj),(aj) - sum_{a != c in P, j != k in Q} Y_(aj),(ck) <= b(b+1).
         Exhaustive over P in the candidate list and all Q with |Q| >= 3.
  lmp1   3-local marginal polytope on side-1 triples (Astra's triples): LP
         separation per triple, cut sum u . (pair blocks) <= 0.
  lmp0   the same on side-0 triples.
  tri    regional triangles (threshold regions of the grid), added in full.
Usage: cp_small.py inst.dat rows cols family [rounds] [K] [eps]
"""
import sys, json, itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from relax_small import *


def Ymat(yv, N):
    return yv.reshape(N, N)   # symmetric


def box_separate(Y, n, Plist, K, minQ=3):
    N = n * n + 1
    X = Y[0, 1:].reshape(n, n)          # X[a, j]
    out = []
    Qmask = np.array([[(s >> j) & 1 for j in range(n)] for s in range(1 << n)], dtype=float)
    Qsize = Qmask.sum(1)
    keepQ = Qsize >= minQ
    Qm = Qmask[keepQ]
    Qids = np.where(keepQ)[0]
    Yb = Y[1:, 1:].reshape(n, n, n, n)  # Yb[a, j, c, k] = Y_(aj),(ck)
    for P in Plist:
        P = list(P)
        p = X[P, :].sum(0)                          # p_j = sum_{a in P} x_aj
        Wfull = Yb[np.ix_(P, range(n), P, range(n))].sum(axis=(0, 2))   # sum_{a,c in P} Y_(aj),(ck)
        # off-diagonal: a != c and j != k; the a == c, j != k terms are gangster 0,
        # the a != c, j == k terms are gangster 0, a == c, j == k is p_j
        W = Wfull - np.diag(np.diag(Wfull))
        mu = Qm @ p
        Eq2 = mu + np.einsum('qj,jk,qk->q', Qm, W, Qm)
        b = np.floor(mu)
        viol = -(Eq2 - (2 * b + 1) * mu + b * (b + 1))
        idx = np.argsort(-viol)[:K]
        for t in idx:
            if viol[t] <= 1e-6:
                break
            out.append((viol[t], ('box', tuple(P), int(Qids[t]), int(b[t]))))
    out.sort(key=lambda r: -r[0])
    return out[:K]


def box_row(n, P, qid, b):
    N = n * n + 1
    Q = [j for j in range(n) if (qid >> j) & 1]
    coefs = {}
    for a in P:
        for j in Q:
            i1 = 1 + a * n + j
            coefs[i1 * N + i1] = coefs.get(i1 * N + i1, 0) + 2 * b
            for c in P:
                if c == a: continue
                for k in Q:
                    if k == j: continue
                    i2 = 1 + c * n + k
                    coefs[i1 * N + i2] = coefs.get(i1 * N + i2, 0) - 1.0
    return coefs, float(b * (b + 1))


def lmp_separate(Y, n, side, K, triples_list=None):
    """3-local marginal polytope on triples of side `side` indices. For side 1:
    triple (j1,j2,j3) of side-1 indices; configurations = injective (i1,i2,i3)."""
    N = n * n + 1
    cfg = list(itertools.permutations(range(n), 3))
    nc = len(cfg)
    pairs_idx = {}
    def var(pr, u, v):
        key = (pr, u, v)
        if key not in pairs_idx:
            pairs_idx[key] = len(pairs_idx)
        return pairs_idx[key]
    # constraint matrix A' u <= 0 : rows = configurations
    rows, cols = [], []
    for c, (u1, u2, u3) in enumerate(cfg):
        for pr, (s, t) in enumerate(((u1, u2), (u1, u3), (u2, u3))):
            rows.append(c); cols.append(var(pr, s, t))
    nu = len(pairs_idx)
    At = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nc, nu))
    keys = sorted(pairs_idx, key=pairs_idx.get)
    out = []
    T = triples_list or list(itertools.combinations(range(n), 3))
    for tr in T:
        e = tr
        def ent(pr, s, t):
            ea, eb = [(e[0], e[1]), (e[0], e[2]), (e[1], e[2])][pr]
            if side == 1:   # e's are side-1 indices; s, t side-0 images
                return 1 + s * n + ea, 1 + t * n + eb
            else:           # e's are side-0 indices; s, t side-1 images
                return 1 + ea * n + s, 1 + eb * n + t
        bvec = np.array([Y[ent(*k)] for k in keys])
        res = linprog(-bvec, A_ub=At, b_ub=np.zeros(nc), bounds=[(-1, 1)] * nu, method='highs')
        v = -res.fun
        if v > 1e-6:
            coefs = {}
            for kk, uu in zip(keys, res.x):
                if abs(uu) > 1e-9:
                    a, bb = ent(*kk)
                    coefs[a * N + bb] = coefs.get(a * N + bb, 0) + uu
            out.append((v, coefs))
    out.sort(key=lambda r: -r[0])
    return out[:K]


def mix_separate(Y, n, side, K):
    """Mixed 3-local polytope: two indices f1 < f2 of side `side` and one
    index s of the other side; configurations (s1, s2, f3) = images of f1, f2
    and the preimage-side occupant of s, consistent with a permutation. Its
    three marginal families are RLT2 equations of BOTH kinds on the same
    cubic variables (see NOTES): a coupled projection of RLT2."""
    N = n * n + 1
    if side == 1:
        lift = lambda s0, s1: 1 + s0 * n + s1      # (side-0 index, side-1 index)
    else:
        lift = lambda s0, s1: 1 + s1 * n + s0      # roles swapped: f's are side-0 indices
    out = []
    # configuration list is the same for every (f1, f2, s) up to relabelling: build per triple
    for f1, f2 in itertools.combinations(range(n), 2):
        for s in range(n):
            var = {}
            def v(key):
                if key not in var: var[key] = len(var)
                return var[key]
            rows, cols = [], []
            ent = {}
            c = 0
            for s1 in range(n):
                for s2 in range(n):
                    if s1 == s2: continue
                    if s1 == s: f3s = [f1]
                    elif s2 == s: f3s = [f2]
                    else: f3s = [f for f in range(n) if f not in (f1, f2)]
                    for f3 in f3s:
                        k12 = ('12', s1, s2); ent[k12] = (lift(s1, f1), lift(s2, f2))
                        k13 = ('13', s1, f3); ent[k13] = (lift(s1, f1), lift(s, f3))
                        k23 = ('23', s2, f3); ent[k23] = (lift(s2, f2), lift(s, f3))
                        for kk in (k12, k13, k23):
                            rows.append(c); cols.append(v(kk))
                        c += 1
            nv = len(var)
            At = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(c, nv))
            keys = sorted(var, key=var.get)
            bvec = np.array([Y[ent[k]] for k in keys])
            res = linprog(-bvec, A_ub=At, b_ub=np.zeros(c), bounds=[(-1, 1)] * nv, method='highs')
            val = -res.fun
            if val > 1e-6:
                coefs = {}
                for kk, uu in zip(keys, res.x):
                    if abs(uu) > 1e-9:
                        a_, b_ = ent[kk]
                        coefs[a_ * N + b_] = coefs.get(a_ * N + b_, 0) + uu
                out.append((val, coefs))
    out.sort(key=lambda r: -r[0])
    return out[:K]


def rows_to_sparse(rowlist, N):
    data, ri, ci, h = [], [], [], []
    for r, (coefs, rhs) in enumerate(rowlist):
        for k, v in coefs.items():
            ri.append(r); ci.append(k); data.append(v)
        h.append(rhs)
    return sp.csr_matrix((data, (ri, ci)), shape=(len(rowlist), N * N)), np.array(h)


def main():
    dat, rows_, cols_, fam = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    rounds = int(sys.argv[5]) if len(sys.argv) > 5 else 6
    K = int(sys.argv[6]) if len(sys.argv) > 6 else 200
    eps = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-6
    n, A, B, C = read_qaplib(dat)
    N = n * n + 1
    L = lifted_cost(n, A, B, C)
    regions = grid_regions(n, rows_, cols_)
    pool = []
    tags = ["dnn"]
    if "tri" in fam.split("+"):
        tags = ["tri"]
    Plist = [S for _, S in regions]
    if "boxall" in fam:
        Plist = [[a for a in range(n) if (s >> a) & 1] for s in range(1, (1 << n) - 1)
                 if 2 <= bin(s).count("1") <= n - 2]
    hist = []
    for r in range(rounds + 1):
        extra = [rows_to_sparse(pool, N)] if pool else None
        t0 = time.time()
        val, st, t, yv = solve(n, L, tags, regions, solver="SCS", extra=extra, eps=eps)
        Y = Ymat(yv, N); Y = (Y + Y.T) / 2
        new = []
        ts = time.time()
        fams = fam.split("+")
        if "box1" in fams or "boxall" in fams:
            for v, (kind, P, qid, b) in box_separate(Y, n, Plist, K):
                new.append((v, box_row(n, list(P), qid, b)))
        if "lmp1" in fams:
            for v, coefs in lmp_separate(Y, n, 1, K):
                new.append((v, (coefs, 0.0)))
        if "lmp0" in fams:
            for v, coefs in lmp_separate(Y, n, 0, K):
                new.append((v, (coefs, 0.0)))
        if "mix1" in fams:
            for v, coefs in mix_separate(Y, n, 1, K):
                new.append((v, (coefs, 0.0)))
        if "mix0" in fams:
            for v, coefs in mix_separate(Y, n, 0, K):
                new.append((v, (coefs, 0.0)))
        new.sort(key=lambda x: -x[0])
        new = new[:K]
        maxv = new[0][0] if new else 0.0
        rec = dict(round=r, value=val, status=st, solve_s=round(t, 1), sep_s=round(time.time() - ts, 1),
                   pool=len(pool), new=len(new), maxviol=maxv)
        print(json.dumps(rec), flush=True)
        hist.append(rec)
        pool += [c for _, c in new]
        if not new:
            break
    print("RESULT", json.dumps(dict(dat=dat, fam=fam, hist=hist)))


if __name__ == "__main__":
    main()
