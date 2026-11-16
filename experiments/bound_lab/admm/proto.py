#!/usr/bin/env python3
"""Prototype of the DNN PRSM (as latest/dnn.C: gamma1 = gamma = 1, beta = F m/3, objective scaled to
max|L| = 1, exact eigendecomposition) with a diagonal congruence scaling of the lifted variable:

    Yt = D Y D,  D = diag(d0, d1 (x) d2)   =>   Yt = Vt Rt Vt',  Vt = orth(D V),  Rt psd,
    min <D^-1 L D^-1, Yt>  s.t.  0 <= Yt_ab <= d_a d_b, Yt_00 = d0^2, gangster zeros.

The R-step stays a plain PSD projection (in the basis Vt); in the original variables the penalty is
entrywise beta d_a^2 d_b^2. The certificate is taken in the ORIGINAL variables from Z = D Zt D
(V'ZV = T'(Vt' Zt Vt)T <= 0 with D V = Vt T), by the Gilmore-Lawler polyhedral part of dnn.C, after
projecting out any positive part of V'ZV (so every value printed is a valid bound, up to float64
rounding: this is a prototype, not Rump-certified).

usage: proto.py inst.dat [-x i:j,...] [-B F] [-p P] [-d0 X] [-i N] [-k N] [-mode row|max]
"""
import sys, time, argparse
import numpy as np
from scipy.optimize import linear_sum_assignment


def read_qap(path):
    t = open(path).read().split()
    n = int(t[0]); v = np.array(t[1:1 + 2 * n * n], float)
    return v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)


def reduce_node(A, B, fixes):
    """Fix facility i at location j (the solver's convention: first matrix distances between
    locations? we only need consistency): objective sum_ik sum_jl A[i,k] B[p(i),p(k)]. Returns
    (A', B', C', const) of the subproblem over the free indices."""
    n = A.shape[0]
    fi = [f[0] for f in fixes]; fj = [f[1] for f in fixes]
    ri = [i for i in range(n) if i not in fi]; rj = [j for j in range(n) if j not in fj]
    const = 0.0
    for a, (i, j) in enumerate(fixes):
        for b, (k, l) in enumerate(fixes):
            const += A[i, k] * B[j, l]
    C = np.zeros((len(ri), len(rj)))
    for (k, l) in fixes:
        C += np.outer(A[ri, k], B[rj, l]) + np.outer(A[k, ri], B[l, rj])
    return A[np.ix_(ri, ri)], B[np.ix_(rj, rj)], C, const


def build(A, B, C):
    m = A.shape[0]; N = m * m
    Q = np.kron(A, B); Q = (Q + Q.T) / 2
    L = np.zeros((N + 1, N + 1)); L[1:, 1:] = Q
    L[0, 1:] = C.reshape(-1) / 2; L[1:, 0] = C.reshape(-1) / 2
    I = np.repeat(np.arange(m), m); J = np.tile(np.arange(m), m)
    gang = ((I[:, None] == I[None, :]) | (J[:, None] == J[None, :])) & ~np.eye(N, dtype=bool)
    free = np.ones((N + 1, N + 1), bool); free[1:, 1:] = ~gang; free[0, 0] = False
    e = np.ones(m) / np.sqrt(m)
    Vb = np.linalg.qr(np.column_stack([e, np.eye(m)[:, :m - 1]]))[0][:, 1:]
    K = (m - 1) ** 2 + 1
    V = np.zeros((N + 1, K)); V[0, 0] = 1 / np.sqrt(2); V[1:, 0] = 1 / (m * np.sqrt(2)); V[1:, 1:] = np.kron(Vb, Vb)
    return L, free, V


def gl_bound(M, m):
    """min <M, Y> over the RLT-1 column-block polytope (dnn.C header, item 2): M00 + LAP(c)."""
    N = m * m
    c = np.zeros((m, m))
    Mx = M[1:, 1:].reshape(m, m, m, m)   # [k, l, i, j]
    for i in range(m):
        ki = [k for k in range(m) if k != i]
        for j in range(m):
            lj = [l for l in range(m) if l != j]
            sub = Mx[np.ix_(ki, lj, [i], [j])][:, :, 0, 0]
            r, cc = linear_sum_assignment(sub)
            a = 1 + i * m + j
            c[i, j] = 2 * M[0, a] + M[a, a] + sub[r, cc].sum()
    r, cc = linear_sum_assignment(c)
    return M[0, 0] + c[r, cc].sum()


def certify(L, V, Z, m, box_free=None):
    """Valid bound from a symmetric Z (original, unscaled-objective units): project V'ZV to nsd."""
    W = V.T @ Z @ V; W = (W + W.T) / 2
    w, U = np.linalg.eigh(W)
    pos = w > 0
    if pos.any():
        Up = U[:, pos]
        Z = Z - V @ ((Up * w[pos]) @ Up.T) @ V.T
    M = L + Z; M = (M + M.T) / 2
    if box_free is not None:
        return M[0, 0] + np.minimum(M[box_free], 0).sum()
    return gl_bound(M, m)


def scaling(A, B, C, p, d0, mode='row'):
    m = A.shape[0]
    if p == 0:
        d1 = np.ones(m); d2 = np.ones(m)
    else:
        As = (np.abs(A) + np.abs(A.T)) / 2; Bs = (np.abs(B) + np.abs(B.T)) / 2
        if mode == 'row':
            r1 = np.sqrt((As ** 2).sum(1)); r2 = np.sqrt((Bs ** 2).sum(1))
        elif mode == 'sum':
            r1 = As.sum(1); r2 = Bs.sum(1)
        else:
            r1 = As.max(1); r2 = Bs.max(1)
        r1 = np.maximum(r1, 1e-3 * r1.max()); r2 = np.maximum(r2, 1e-3 * r2.max())
        d1 = (r1 / np.exp(np.log(r1).mean())) ** p; d2 = (r2 / np.exp(np.log(r2).mean())) ** p
    return np.concatenate([[d0], np.kron(d1, d2)])


def prsm(A, B, C, const, F=0.21, p=0.0, d0=1.0, iters=1500, every=10, mode='row', verbose=True, cutoff=None,
         start=None, d=None, box=False, tag='', snap=None):
    m = A.shape[0]
    L, free, V = build(A, B, C)
    scale = max(1.0, np.abs(L).max()); Ls = L / scale
    if d is None:
        d = scaling(A, B, C, p, d0, mode)
    Vt = np.linalg.qr(d[:, None] * V)[0]
    Lt = Ls / np.outer(d, d)
    Ub = np.outer(d, d)
    Lt_max = np.abs(Lt).max()
    beta = F * m / 3.0
    n1 = L.shape[0]
    Y = np.zeros((n1, n1)); Y[0, 0] = d[0] ** 2
    Z = np.zeros((n1, n1))
    if start is not None:
        Y = start[0] * Ub; Z = start[1] / Ub
    best = -np.inf; hist = []
    t0 = time.process_time()
    for it in range(1, iters + 1):
        W = Vt.T @ (Y + Z / beta) @ Vt; W = (W + W.T) / 2
        w, U = np.linalg.eigh(W)
        pos = w > 0
        Up = U[:, pos]
        R = (Up * w[pos]) @ Up.T
        VRV = Vt @ R @ Vt.T
        check = it % every == 0 or it == iters
        if check or (snap is not None and it in snap):
            Zc = Z + beta * (Y - VRV)
            if snap is not None and it in snap: snap[it] = Zc * Ub
        Z = Z + beta * (Y - VRV)                 # gamma1 = 1
        Yn = np.clip(VRV - (Lt + Z) / beta, 0.0, Ub)
        Yn[~free] = 0.0; Yn[0, 0] = d[0] ** 2
        res = Yn - VRV
        Z = Z + beta * res                        # gamma = 1
        dY = np.linalg.norm((Yn - Y) / Ub); Y = Yn
        if check:
            lb = certify(Ls, V, Zc * Ub, m, free if box else None) * scale + const
            best = max(best, lb)
            pres = np.linalg.norm(res / Ub)
            obj = (Lt * Y).sum() * scale + const
            hist.append((it, best, lb, obj, pres, dY, int(pos.sum())))
            if verbose:
                print('it %5d lb %.4f (cert %.4f) primal %.4f pres %.2e dY %.2e npos %d cpu %.1f %s' %
                      (it, best, lb, obj, pres, dY, pos.sum(), time.process_time() - t0, tag), flush=True)
            if cutoff is not None and best > cutoff:
                break
    return dict(hist=hist, Y=Y / Ub, Z=Z * Ub, best=best)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('file'); ap.add_argument('-x', default='')
    ap.add_argument('-B', type=float, default=0.21); ap.add_argument('-p', type=float, default=0.0)
    ap.add_argument('-d0', type=float, default=1.0); ap.add_argument('-i', type=int, default=1500)
    ap.add_argument('-k', type=int, default=10); ap.add_argument('-mode', default='row')
    ap.add_argument('-box', action='store_true')
    a = ap.parse_args()
    A, B = read_qap(a.file)
    B = (B + B.T) / 2
    fixes = [tuple(int(v) for v in s.split(':')) for s in a.x.split(',') if s]
    A2, B2, C2, const = reduce_node(A, B, fixes)
    print('instance %s m=%d const=%.1f F=%g p=%g d0=%g mode=%s' % (a.file, A2.shape[0], const, a.B, a.p, a.d0, a.mode))
    prsm(A2, B2, C2, const, F=a.B, p=a.p, d0=a.d0, iters=a.i, every=a.k, mode=a.mode, box=a.box)
