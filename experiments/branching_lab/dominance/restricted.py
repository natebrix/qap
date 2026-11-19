"""Step 4: the pair test restricted to completions that could beat the cutoff.
For a node N, a pair (i, k) of its fixed rows and mu >= 0:
   m(mu) = min over completions p of N of  Delta(p) + mu * (G(p) - UB)
where G(p) <= cost(p) is the Gilmore-Lawler linearisation of cost on N (fixed-fixed constant, fixed-free and
diagonal terms exact, free-free terms bounded below row by row by a minimal scalar product). If m(mu) >= 1 for
some mu, every completion with cost(p) < UB has Delta(p) >= 1 (Delta >= m - mu (G - UB) >= m - mu (cost - UB) > m),
so N holds no solution below UB that is not strictly dominated: dropping N is as valid as fathoming it.
mu -> infinity is GLB fathoming (G >= UB), so the test is only worth anything where GLB alone does not fathom.
Requires A symmetric (tai-b, tho); B is symmetrised, which leaves cost unchanged.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

def glb_parts(A, B, F, UB):
    n = len(A); Bs = (B + B.T) / 2.0; Af = A.astype(float)
    fr = np.array(list(F.keys())); fc = np.array(list(F.values()))
    mask = np.ones(n, bool); mask[fr] = False; R = np.nonzero(mask)[0]
    mask = np.ones(n, bool); mask[fc] = False; C = np.nonzero(mask)[0]
    const = float((Af[np.ix_(fr, fr)] * Bs[np.ix_(fc, fc)]).sum())
    lin = 2 * Af[np.ix_(R, fr)] @ Bs[np.ix_(fc, C)] + np.outer(np.diag(Af)[R], np.diag(Bs)[C])
    m = len(R)
    if m > 1:
        a = Af[np.ix_(R, R)].copy(); a[np.eye(m, dtype=bool)] = np.inf; a = np.sort(a, 1)[:, :m - 1]
        b = Bs[np.ix_(C, C)].copy(); b[np.eye(m, dtype=bool)] = -np.inf; b = -np.sort(-b, 1)[:, :m - 1]
        lin = lin + a @ b.T
    return const, lin, R, C

def pair_lin(A, B, F, i, k, R, C):
    a, b = F[i], F[k]
    D0 = (A[i, i] - A[k, k]) * (B[a, a] - B[b, b]) + (A[i, k] - A[k, i]) * (B[a, b] - B[b, a])
    const = float(D0)
    for l, c in F.items():
        if l in (i, k): continue
        const += float((A[i, l] - A[k, l]) * (B[a, c] - B[b, c]) + (A[l, i] - A[l, k]) * (B[c, a] - B[c, b]))
    M = np.outer(A[i, R] - A[k, R], B[a, C] - B[b, C]) + np.outer(A[R, i] - A[R, k], B[C, a] - B[C, b])
    return const, M.astype(float)

def restricted_test(A, B, F, UB, pairs, mus=(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)):
    """returns (best m over pairs and mu, its pair and mu, GLB value of N)"""
    gc, G, R, C = glb_parts(A, B, F, UB)
    r, c = linear_sum_assignment(G); glb = gc + G[r, c].sum()
    best = (-np.inf, None, None)
    for (i, k) in pairs:
        dc, D = pair_lin(A, B, F, i, k, R, C)
        for mu in mus:
            M = D + mu * G
            r, c = linear_sum_assignment(M)
            v = dc + M[r, c].sum() + mu * (gc - UB)
            if v > best[0]: best = (v, (i, k), mu)
    return best[0], best[1], best[2], glb
