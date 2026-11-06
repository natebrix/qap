#!/usr/bin/env python3
"""Nonlinear-constraint candidates at a dumped wil50 DNN iterate (README.md).

  python3 nlc_analyze.py SPEC Y.bin [--label L] [--kmax 5]

SPEC is the node ("cell:facility,..." on qaplib/wil50.dat, first index = grid
cell, second = facility), Y.bin a -dumpY file (raw doubles, (m^2+1)^2 row-major,
lifted index of x_(a,i) = 1 + a*m + i with a the local free cell and i the
local free facility, both in original order). Prints, for each candidate family,
its violation at the iterate (positive = violated):

  1. feasibility residuals of the iterate (so that small violations can be
     read against them)
  2. spectrum of Y (how far from rank one)
  3. majorization of E[X B X'] by lambda(B) and of E[X' A X] by lambda(A)
     (implied by the DNN: a check of the proof)
  4. the partial transpose Y_xx^{T_B}: -I <= Y^{T_B} <= I for every permutation
     lift (an orthogonal involution); the operator norm of the realignment
     E[X (x) X] (<= 1, implied by double stochasticity)
  5. k-partition clique (pigeonhole) inequalities on the rows, the columns and
     the mid5 super-cells: any t facilities with t > #parts have at least f(t)
     same-part pairs (f by DP over the parts' capacities)
  6. spread inequalities: sum over pairs of a k-set of E d_ij <= Dmax(k), the
     largest total pairwise distance of k free cells (exact by enumeration)
  7. regional triangle max violation per threshold (identifies the iterate)
"""
import argparse
import itertools
import sys
import time

import numpy as np

ROWS, COLS = 5, 10
WIL = '/home/user/qap/qaplib/wil50.dat'


def read_qaplib(path):
    tok = open(path).read().split()
    n = int(tok[0]); v = list(map(float, tok[1:1 + 2 * n * n]))
    A = np.array(v[:n * n]).reshape(n, n)
    B = np.array(v[n * n:2 * n * n]).reshape(n, n)
    return n, A, B


def parse_spec(spec):
    if spec in ('-', ''):
        return []
    return [tuple(map(int, t.split(':'))) for t in spec.split(',')]


class Node:
    def __init__(self, spec, path=WIL):
        n, A, B = read_qaplib(path)
        A = (A + A.T) / 2; B = (B + B.T) / 2
        fx = parse_spec(spec)
        self.n, self.fixes = n, fx
        fl = {a for a, b in fx}; ff = {b for a, b in fx}
        self.loc = [a for a in range(n) if a not in fl]
        self.fac = [b for b in range(n) if b not in ff]
        self.m = m = len(self.loc)
        self.A = A[np.ix_(self.loc, self.loc)]
        self.B = B[np.ix_(self.fac, self.fac)]
        self.r = np.array([a // COLS for a in self.loc])
        self.c = np.array([a % COLS for a in self.loc])
        self.Dr = np.abs(self.r[:, None] - self.r[None, :]).astype(float)
        self.Dc = np.abs(self.c[:, None] - self.c[None, :]).astype(float)
        C = np.zeros((m, m))
        for k, pk in fx:
            C += np.outer(A[self.loc, k], B[self.fac, pk]) + np.outer(A[k, self.loc], B[pk, self.fac])
        self.C = C
        self.shift = sum(A[i, k] * B[pi, pk] for i, pi in fx for k, pk in fx)

    def load(self, path):
        m = self.m; n1 = m * m + 1
        Y = np.fromfile(path, dtype=np.float64)
        assert Y.size == n1 * n1, (Y.size, n1)
        self.Y = Y = Y.reshape(n1, n1)
        self.x = Y[0, 1:].reshape(m, m)                  # [a, i]
        self.T = Y[1:, 1:].reshape(m, m, m, m)           # [a, i, b, j] = Y_(ai),(bj)
        return self


def residuals(nd):
    m, Y, T, x = nd.m, nd.Y, nd.T, nd.x
    d = np.einsum('aiai->ai', T)
    out = {}
    out['min entry'] = float(Y.min())
    out['max |diag - row0|'] = float(np.abs(d - x).max())
    out['max |row sums of x - 1|'] = float(max(np.abs(x.sum(0) - 1).max(), np.abs(x.sum(1) - 1).max()))
    # lifted row sums: sum_b Y_(ai),(bj) = x_ai for every j (b over all cells, gangster included)
    rs = T.sum(axis=2)                                   # [a, i, j]
    out['max |lifted row sums - x|'] = float(np.abs(rs - x[:, :, None]).max())
    cs = T.sum(axis=3)                                   # [a, i, b]
    out['max |lifted col sums - x|'] = float(np.abs(cs - x[:, :, None]).max())
    g1 = np.abs(T[np.arange(m), :, np.arange(m), :])     # [a, i, j]: same cell, i != j
    g1 = g1 - np.einsum('aij->aij', np.einsum('aii->ai', g1)[:, :, None] * np.eye(m)[None, :, :])
    out['max gangster (same cell)'] = float(np.abs(g1).max())
    g2 = T[:, np.arange(m), :, np.arange(m)]              # [i, a, b]: same facility, a != b
    g2 = g2 - np.einsum('iab->iab', np.einsum('iaa->ia', g2)[:, :, None] * np.eye(m)[None, :, :])
    out['max gangster (same facility)'] = float(np.abs(g2).max())
    return out


def majorization(nd):
    """lambda(E[X B X']) majorized by lambda(B) and lambda(E[X' A X]) by lambda(A):
    returns the largest excess of a top-k partial sum (<= 0 means satisfied)"""
    T = nd.T
    M = np.einsum('ij,aibj->ab', nd.B, T, optimize=True)   # E[X B X'] on the cells
    N = np.einsum('ab,aibj->ij', nd.A, T, optimize=True)   # E[X' A X] on the facilities
    res = {}
    for name, Mat, Ref in (('E[XBX] vs B', M, nd.B), ('E[XAX] vs A', N, nd.A)):
        Mat = (Mat + Mat.T) / 2
        lm = np.sort(np.linalg.eigvalsh(Mat))[::-1]
        lr = np.sort(np.linalg.eigvalsh(Ref))[::-1]
        top = np.cumsum(lm) - np.cumsum(lr)
        bot = np.cumsum(lr[::-1]) - np.cumsum(lm[::-1])
        res[name] = dict(max_top_excess=float(top.max()), max_bottom_excess=float(bot.max()),
                         trace_diff=float(lm.sum() - lr.sum()),
                         frob_M=float(np.linalg.norm(Mat)), frob_ref=float(np.linalg.norm(Ref)))
    res['objective <L,Y>+shift'] = float((nd.B * N).sum() + (nd.C * nd.x).sum() + nd.shift)
    return res, N


def partial_transpose(nd):
    m, T = nd.m, nd.T
    PT = np.ascontiguousarray(T.transpose(0, 3, 2, 1)).reshape(m * m, m * m)   # [(a,i),(b,j)] = Y_(aj),(bi)
    asym = float(np.abs(PT - PT.T).max())
    PT = (PT + PT.T) / 2
    t0 = time.time()
    ev = np.linalg.eigvalsh(PT)
    out = dict(lam_max=float(ev[-1]), lam_min=float(ev[0]), trace=float(np.trace(PT)),
               asym=asym, n_above_1=int((ev > 1 + 1e-6).sum()), n_below_m1=int((ev < -1 - 1e-6).sum()),
               top5=[float(v) for v in ev[-5:][::-1]], bottom5=[float(v) for v in ev[:5]],
               eig_seconds=time.time() - t0)
    # the realignment E[X (x) X]: rows (a,b), columns (i,j)
    R = np.ascontiguousarray(T.transpose(0, 2, 1, 3)).reshape(m * m, m * m)
    sv = np.linalg.svd(R, compute_uv=False)
    out['realign_sigma_max'] = float(sv[0]); out['realign_sigma_2'] = float(sv[1])
    return out, ev


def spectrum_Y(nd):
    ev = np.sort(np.linalg.eigvalsh((nd.Y + nd.Y.T) / 2))[::-1]
    return dict(top8=[float(v) for v in ev[:8]], n_above_1e_2=int((ev > 1e-2).sum()),
                n_above_1e_3=int((ev > 1e-3).sum()), min=float(ev[-1]), trace=float(ev.sum()))


def same_part_matrix(nd, parts):
    """P_ij = probability that free facilities i, j sit in the same part (zero diagonal)"""
    m, T = nd.m, nd.T
    P = np.zeros((m, m))
    for cells in parts:
        if not cells:
            continue
        idx = np.array(cells)
        P += T[np.ix_(idx, np.arange(m), idx, np.arange(m))].sum(axis=(0, 2))
    np.fill_diagonal(P, 0.0)
    return P


def f_pairs(caps, t):
    """min over integer n_p in [0, cap_p], sum n_p = t of sum_p C(n_p, 2) (DP)"""
    INF = 10 ** 9
    best = {0: 0}
    for cap in caps:
        nb = {}
        for s, v in best.items():
            for k in range(0, cap + 1):
                if s + k > t:
                    break
                w = v + k * (k - 1) // 2
                if w < nb.get(s + k, INF):
                    nb[s + k] = w
        best = nb
    return best.get(t, None)


def min_subset_sum(P, t, starts=40, seed=1):
    """heuristic min over t-sets of sum_{i<j in set} P_ij: greedy from each of the
    `starts` smallest pairs plus random, then 1-swap local search"""
    m = P.shape[0]
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(m, 1)
    order = np.argsort(P[iu])
    best = (np.inf, None)
    seeds = [(iu[0][k], iu[1][k]) for k in order[:starts]]
    for s in range(starts):
        if s < len(seeds):
            S = list(seeds[s])
        else:
            S = list(rng.choice(m, 2, replace=False))
        inS = np.zeros(m, bool); inS[S] = True
        while len(S) < t:
            gain = P[:, S].sum(1); gain[inS] = np.inf
            k = int(np.argmin(gain)); S.append(k); inS[k] = True
        val = P[np.ix_(S, S)].sum() / 2
        improved = True
        while improved:
            improved = False
            for pos in range(len(S)):
                i = S[pos]
                rest = [j for j in S if j != i]
                base = P[np.ix_(rest, rest)].sum() / 2
                add = P[:, rest].sum(1); add[inS] = np.inf; add[i] = P[i, rest].sum()
                k = int(np.argmin(add))
                if base + add[k] < val - 1e-12:
                    inS[i] = False; inS[k] = True; S[pos] = k; val = base + add[k]; improved = True
        if val < best[0]:
            best = (val, sorted(S))
    return best


def max_subset_sum(D, k, starts=40, seed=2):
    return tuple(-v if isinstance(v, float) else v for v in min_subset_sum(-D, k, starts, seed))


def dmax(Dcells, k, limit=12_000_000):
    """exact max over k-sets of free cells of the total pairwise distance"""
    m = Dcells.shape[0]
    from math import comb
    if comb(m, k) > limit:
        return None
    iu = np.triu_indices(m, 1)
    best = -np.inf
    chunk = []
    def flush(chunk):
        nonlocal best
        C = np.array(chunk)
        tot = np.zeros(len(C))
        for p in range(k):
            for q in range(p + 1, k):
                tot += Dcells[C[:, p], C[:, q]]
        best = max(best, float(tot.max()))
    for comb_ in itertools.combinations(range(m), k):
        chunk.append(comb_)
        if len(chunk) >= 500000:
            flush(chunk); chunk = []
    if chunk:
        flush(chunk)
    return best


def thresholds(nd):
    regs = {}
    for t in range(1, ROWS):
        regs[f'row<{t}'] = [a for a in range(nd.m) if nd.r[a] < t]
    for t in range(1, COLS):
        regs[f'col<{t}'] = [a for a in range(nd.m) if nd.c[a] < t]
    return regs


def triangle_max(nd, cells):
    m, T = nd.m, nd.T
    idx = np.array(cells)
    W = T[np.ix_(idx, np.arange(m), idx, np.arange(m))].sum(axis=(0, 2))
    p = np.diag(W).copy()
    best = -np.inf
    for i in range(m):
        for j in range(i + 1, m):
            wij = W[i, j]
            wik = W[i, j + 1:]; wjk = W[j, j + 1:]; pk = p[j + 1:]
            if len(pk) == 0:
                continue
            v = np.max([
                (wij + wik - wjk - p[i]).max(), (wij + wjk - wik - p[j]).max(),
                (wik + wjk - wij - pk).max(), (p[i] + p[j] + pk - wij - wik - wjk - 1).max()])
            best = max(best, float(v))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('ybin')
    ap.add_argument('--label', default='')
    ap.add_argument('--kmax', type=int, default=5)
    ap.add_argument('--skip-eig', action='store_true')
    a = ap.parse_args()
    nd = Node(a.spec).load(a.ybin)
    m = nd.m
    print(f'# {a.label} spec {a.spec} m {m} file {a.ybin}')
    t0 = time.time()
    print('## 1. residuals')
    for k, v in residuals(nd).items():
        print(f'  {k:32s} {v: .3e}')
    print('## 2. spectrum of Y')
    for k, v in spectrum_Y(nd).items():
        print(f'  {k:16s} {v}')
    print('## 3. majorization (implied; excess <= 0 up to residuals)')
    maj, Delta = majorization(nd)
    for k, v in maj.items():
        print(f'  {k:16s} {v}')
    if not a.skip_eig:
        print('## 4. partial transpose Y_xx^{T_B} (valid: -1 <= lambda <= 1) and realignment')
        pt, ev = partial_transpose(nd)
        for k, v in pt.items():
            print(f'  {k:18s} {v}')
        np.save(a.ybin.rsplit('/', 1)[-1].replace('.Y.bin', '') + '.ptspec.npy', ev) if False else None
    print('## 5. k-partition clique (pigeonhole) inequalities: slack = min same-part pairs - f(t)')
    partitions = {
        'rows': [[a_ for a_ in range(m) if nd.r[a_] == r] for r in range(ROWS)],
        'cols': [[a_ for a_ in range(m) if nd.c[a_] == c] for c in range(COLS)],
        'mid5 super-cells': [[a_ for a_ in range(m) if rg(nd.r[a_]) == R and cg(nd.c[a_]) == Cc]
                             for R in range(3) for Cc in range(4)],
    }
    for name, parts in partitions.items():
        parts = [p for p in parts if p]
        caps = [len(p) for p in parts]
        P = same_part_matrix(nd, parts)
        np_ = len(parts)
        line = []
        for t in range(np_ + 1, min(np_ + 6, m) + 1):
            ft = f_pairs(caps, t)
            if ft is None:
                continue
            val, S = min_subset_sum(P, t)
            line.append(f't={t}: f={ft} min={val:.3f} slack={val - ft:+.3f}')
        print(f'  {name} ({np_} parts, caps {caps}):')
        for l in line:
            print('    ' + l)
    print('## 6. spread inequalities: slack = Dmax(k) - max over k-sets of sum E d_ij')
    for name, Dc_, Dm in (('Manhattan', nd.A, Delta),
                          ('rows only', nd.Dr, np.einsum('ab,aibj->ij', nd.Dr, nd.T, optimize=True)),
                          ('cols only', nd.Dc, np.einsum('ab,aibj->ij', nd.Dc, nd.T, optimize=True))):
        Dm = (Dm + Dm.T) / 2; np.fill_diagonal(Dm, 0.0)
        for k in range(3, a.kmax + 1):
            dm = dmax(Dc_, k)
            if dm is None:
                continue
            val, S = max_subset_sum(Dm, k)
            print(f'  {name:10s} k={k}: Dmax={dm:.0f} max E sum={val:.3f} slack={dm - val:+.3f} set={S}')
    print('## 7. regional triangle max violation per threshold')
    for name, cells in thresholds(nd).items():
        if 0 < len(cells) < m:
            print(f'  {name:7s} cap {len(cells):2d} max {triangle_max(nd, cells):+.4f}')
    print(f'# done in {time.time() - t0:.1f} s')


def rg(r):
    return 0 if r < 2 else (1 if r < 3 else 2)


def cg(c):
    return 0 if c < 4 else (1 if c < 5 else (2 if c < 6 else 3))


if __name__ == '__main__':
    main()
