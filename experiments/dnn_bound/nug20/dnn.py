"""DNN (doubly nonnegative) SDP relaxation of a QAP subproblem

    min  tr(A Y B Y') + <C, Y> + const   over permutation matrices Y (m x m)

via the facially-reduced formulation of Zhao-Karisch-Rendl-Wolkowicz / Oliveira-Wolkowicz-Xu:

    min <L, Yh>  s.t.  Yh = Vh R Vh',  R psd,  Yh_00 = 1,  Yh_J = 0 (gangster),  0 <= Yh <= 1

with Vh orthonormal ((m^2+1) x ((m-1)^2+1)).  Solved by the OWX ADMM.  A VALID lower bound is
computed from the dual multiplier Z via the Lagrangian dual (projecting Vh' Z Vh onto the
negative semidefinite cone), which never needs the ADMM to have converged.
"""
import numpy as np
import time


def build_lifted(A, B, C):
    m = A.shape[0]
    N = m * m
    # index (i,j) -> i*m + j ; tr(A Y B Y') = sum_{ijkl} A_ik B_jl y_ij y_kl = x' (A kron B) x
    Q = np.kron(A, B)
    L = np.zeros((N + 1, N + 1))
    L[1:, 1:] = Q
    L[0, 1:] = C.reshape(-1) / 2.0
    L[1:, 0] = C.reshape(-1) / 2.0
    # gangster mask (entries forced to zero) within the x-block
    I = np.repeat(np.arange(m), m)   # facility of index
    J = np.tile(np.arange(m), m)     # location of index
    same_fac = I[:, None] == I[None, :]
    same_loc = J[:, None] == J[None, :]
    diag = np.eye(N, dtype=bool)
    gang = (same_fac | same_loc) & ~diag
    free = np.ones((N + 1, N + 1), dtype=bool)
    free[1:, 1:] = ~gang
    free[0, 0] = False      # fixed to 1
    # Vh
    e = np.ones(m) / np.sqrt(m)
    # orthonormal basis of e-perp
    Vb = np.linalg.qr(np.column_stack([e, np.eye(m)[:, :m - 1]]))[0][:, 1:]   # m x (m-1)
    assert np.allclose(Vb.T @ np.ones(m), 0) and np.allclose(Vb.T @ Vb, np.eye(m - 1))
    K = (m - 1) ** 2 + 1
    Vh = np.zeros((N + 1, K))
    Vh[0, 0] = 1 / np.sqrt(2)
    Vh[1:, 0] = 1 / (m * np.sqrt(2))
    Vh[1:, 1:] = np.kron(Vb, Vb)
    assert np.allclose(Vh.T @ Vh, np.eye(K), atol=1e-10)
    return L, free, Vh


def check_lifting(A, B, C, L, Vh, rng, trials=3):
    """Permutations are feasible for the relaxation and the objective matches."""
    m = A.shape[0]
    for _ in range(trials):
        p = rng.permutation(m)
        X = np.zeros((m, m)); X[np.arange(m), p] = 1
        xh = np.concatenate([[1.0], X.reshape(-1)])
        Y = np.outer(xh, xh)
        val = xh @ L @ xh
        true = np.sum(A * B[np.ix_(p, p)]) + C[np.arange(m), p].sum()
        assert abs(val - true) < 1e-8, (val, true)
        # in range of Vh
        assert np.allclose(Vh @ (Vh.T @ xh), xh, atol=1e-9)


def proj_psd(M, sign=+1.0):
    w, U = np.linalg.eigh((M + M.T) / 2)
    if sign > 0:
        w = np.maximum(w, 0.0)
    else:
        w = np.minimum(w, 0.0)
    return (U * w) @ U.T


def dual_bound(L, free, Vh, Z):
    """Valid lower bound on the DNN value from any symmetric Z (see module docstring).
    Returns (lb, Ztilde)."""
    W = Vh.T @ Z @ Vh
    Zt = Z - Vh @ proj_psd(W, +1.0) @ Vh.T      # Vh' Zt Vh = negative part  (psd-negative)
    M = L + Zt
    M = (M + M.T) / 2
    lb = M[0, 0] + np.minimum(M[free], 0.0).sum()
    return lb, Zt


def dnn_admm(A, B, C, const=0.0, beta=None, gamma=1.618, maxit=20000, check_every=25,
             fathom_at=None, tol=1e-5, time_limit=None, verbose=False, stall_window=None):
    """Returns dict with certified lower bound 'lb' (includes const), primal 'obj', iterations, time."""
    m = A.shape[0]
    L, free, Vh = build_lifted(A, B, C)
    N1 = L.shape[0]
    scale = max(1.0, np.abs(L).max())
    Ls = L / scale
    if beta is None:
        beta = m / 3.0  # OWX-style default; objective is normalised by max|L|
    t0 = time.process_time()
    Y = np.zeros((N1, N1)); Y[0, 0] = 1.0
    Z = np.zeros((N1, N1))
    best_lb = -np.inf
    hist = []
    it = 0
    stall = 0
    last_lbs = []
    status = 'maxit'
    while it < maxit:
        # R-step
        W = Vh.T @ (Y + Z / beta) @ Vh
        R = proj_psd(W, +1.0)
        VRV = Vh @ R @ Vh.T
        # Y-step
        Ynew = VRV - (Ls + Z) / beta
        Ynew = np.clip(Ynew, 0.0, 1.0)
        Ynew[~free] = 0.0
        Ynew[0, 0] = 1.0
        Ynew = (Ynew + Ynew.T) / 2
        # Z-step
        resid = Ynew - VRV
        Z = Z + gamma * beta * resid
        dY = np.linalg.norm(Ynew - Y)
        Y = Ynew
        it += 1
        if it % check_every == 0 or it == maxit:
            lb, _ = dual_bound(Ls, free, Vh, Z)
            lb = lb * scale + const
            obj = np.sum(Ls * Y) * scale + const
            pres = np.linalg.norm(resid)
            best_lb = max(best_lb, lb)
            hist.append((it, lb, obj, pres, dY))
            if verbose:
                print('it %5d lb %.3f obj %.3f pres %.2e dY %.2e  t=%.1f' % (it, lb, obj, pres, dY, time.process_time() - t0), flush=True)
            if fathom_at is not None and best_lb > fathom_at:
                status = 'fathomed'
                break
            if pres < tol and dY < tol:
                status = 'converged'
                break
            if time_limit is not None and time.process_time() - t0 > time_limit:
                status = 'timelimit'
                break
            if stall_window is not None:
                last_lbs.append(lb)
                if len(last_lbs) > stall_window:
                    last_lbs.pop(0)
                    if max(last_lbs) - min(last_lbs) < 0.05 and abs(obj - lb) < 0.5:
                        status = 'stalled'
                        break
    # final bound from the last Z (rigorous), plus the best seen
    lb_final, _ = dual_bound(Ls, free, Vh, Z)
    lb_final = lb_final * scale + const
    best_lb = max(best_lb, lb_final)
    return dict(lb=best_lb, lb_last=lb_final, obj=np.sum(Ls * Y) * scale + const, it=it,
                cpu=time.process_time() - t0, status=status, hist=hist, m=m)
