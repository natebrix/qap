"""Exchange dominance tests for QAP nodes.

cost(p) = sum_{u,v} A[u,v] B[p[u],p[v]].  A node is a partial assignment F = {row: col}.
Pair exchange of fixed rows i, k (p[i]=a, p[k]=b): p' = p with p'[i]=b, p'[k]=a.
Delta(p) = cost(p) - cost(p') = D0 + sum_{l != i,k} D_l(p[l]), with
  D_l(c) = (A[i,l]-A[k,l]) (B[a,c]-B[b,c]) + (A[l,i]-A[l,k]) (B[c,a]-B[c,b])
  D0     = (A[i,i]-A[k,k]) (B[a,a]-B[b,b]) + (A[i,k]-A[k,i]) (B[a,b]-B[b,a]).
Min over completions = D0 + sum_{fixed l} D_l(F[l]) + LAP over free rows x free cols of D_l(c).
If min >= 1 (integer data) every completion is strictly worse than its exchange image: no optimum in F.
When A or B is symmetric, the LAP matrix is rank one, outer(x, y), and the LAP is a sort
(rearrangement: x ascending against y descending).
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

def pair_min_general(A, B, F, i, k):
    a, b = F[i], F[k]
    n = len(A)
    rows = [l for l in range(n) if l not in F]
    used = set(F.values()); cols = [c for c in range(n) if c not in used]
    D0 = (A[i, i] - A[k, k]) * (B[a, a] - B[b, b]) + (A[i, k] - A[k, i]) * (B[a, b] - B[b, a])
    const = int(D0)
    for l, c in F.items():
        if l in (i, k): continue
        const += int((A[i, l] - A[k, l]) * (B[a, c] - B[b, c]) + (A[l, i] - A[l, k]) * (B[c, a] - B[c, b]))
    if not rows:
        return const
    R = np.array(rows); C = np.array(cols)
    M = np.outer(A[i, R] - A[k, R], B[a, C] - B[b, C]) + np.outer(A[R, i] - A[R, k], B[C, a] - B[C, b])
    r, c = linear_sum_assignment(M)
    return const + int(M[r, c].sum())

class Inst:
    def __init__(self, A, B):
        self.A = np.asarray(A, dtype=np.int64); self.B = np.asarray(B, dtype=np.int64)
        self.n = len(A)
        self.Asym = bool((self.A == self.A.T).all()); self.Bsym = bool((self.B == self.B.T).all())
        # rank one form: if A symmetric, M = outer(alpha, beta+beta'); if B symmetric, outer(alpha+alpha', beta)
        self.rank1 = self.Asym or self.Bsym

    def pair_mins(self, F):
        """min Delta over completions for every pair of fixed rows (i<k by order in F). Returns
        (array of mins, list of pairs)."""
        A, B, n = self.A, self.B, self.n
        fr = np.array(list(F.keys()), dtype=np.int64); fc = np.array(list(F.values()), dtype=np.int64)
        d = len(fr)
        if d < 2: return np.zeros(0, dtype=np.int64), []
        mask = np.ones(n, bool); mask[fr] = False; R = np.nonzero(mask)[0]
        mask = np.ones(n, bool); mask[fc] = False; C = np.nonzero(mask)[0]
        I, K = np.triu_indices(d, 1)
        i, k, a, b = fr[I], fr[K], fc[I], fc[K]
        D0 = (A[i, i] - A[k, k]) * (B[a, a] - B[b, b]) + (A[i, k] - A[k, i]) * (B[a, b] - B[b, a])
        # fixed l contributions (including l = i, k gives zero? check: l=i: (A[i,i]-A[k,i])(B[a,a]-B[b,a]) not zero) -> mask them
        Ai = A[i][:, fr] - A[k][:, fr]           # P x d : A[i,l]-A[k,l]
        Ai2 = A[fr][:, i].T - A[fr][:, k].T      # P x d : A[l,i]-A[l,k]
        Bc = B[a][:, fc] - B[b][:, fc]           # B[a,c]-B[b,c]
        Bc2 = B[fc][:, a].T - B[fc][:, b].T      # B[c,a]-B[c,b]
        T = Ai * Bc + Ai2 * Bc2
        P = len(I); T[np.arange(P), I] = 0; T[np.arange(P), K] = 0
        const = D0 + T.sum(1)
        if len(R) == 0:
            return const, list(zip(i, k))
        if self.Asym:
            x = A[i][:, R] - A[k][:, R]
            y = (B[a][:, C] - B[b][:, C]) + (B[C][:, a].T - B[C][:, b].T)
        elif self.Bsym:
            x = (A[i][:, R] - A[k][:, R]) + (A[R][:, i].T - A[R][:, k].T)
            y = B[a][:, C] - B[b][:, C]
        else:
            out = np.empty(P, dtype=np.int64)
            for p in range(P):
                M = np.outer(A[i[p], R] - A[k[p], R], B[a[p], C] - B[b[p], C]) + \
                    np.outer(A[R, i[p]] - A[R, k[p]], B[C, a[p]] - B[C, b[p]])
                r, c = linear_sum_assignment(M); out[p] = M[r, c].sum()
            return const + out, list(zip(i, k))
        lap = (np.sort(x, 1) * np.sort(y, 1)[:, ::-1]).sum(1)
        return const + lap, list(zip(i, k))

    def pair_test(self, F):
        m, pairs = self.pair_mins(F)
        if len(m) == 0: return None, None
        j = int(np.argmax(m)); return int(m[j]), pairs[j]

    def onefixed(self, F, ks=None):
        """One-fixed test: k fixed at a, i free. For each free col c, min of Delta_(i,k) over completions of
        F + {i: c}. N dominated by (k, i) if that min >= 1 for every c. Returns best (min over c) over
        (k,i), and the pair. Only rank-one instances (sort form)."""
        A, B, n = self.A, self.B, self.n
        if not self.rank1:
            return self.onefixed_general(F)
        fr = np.array(list(F.keys()), dtype=np.int64); fc = np.array(list(F.values()), dtype=np.int64)
        mask = np.ones(n, bool); mask[fr] = False; R = np.nonzero(mask)[0]
        mask = np.ones(n, bool); mask[fc] = False; C = np.nonzero(mask)[0]
        m = len(R)
        if m < 1 or len(fr) < 1: return None, None
        best, arg = None, None
        for kk in (range(len(fr)) if ks is None else ks):
            k, a = fr[kk], fc[kk]
            for ii in range(m):
                i = R[ii]
                # subnode: i -> c for each c in C. pair (i at c, k at a). Exchange gives i->a, k->c.
                # D0(c) = (A[i,i]-A[k,k])(B[c,c]-B[a,a]) + (A[i,k]-A[k,i])(B[c,a]-B[a,c])
                D0 = (A[i, i] - A[k, k]) * (B[C, C] - B[a, a]) + (A[i, k] - A[k, i]) * (B[C, a] - B[a, C])
                # fixed l != k: (A[i,l]-A[k,l])(B[c,F[l]]-B[a,F[l]]) + (A[l,i]-A[l,k])(B[F[l],c]-B[F[l],a])
                o = np.arange(len(fr)) != kk; fl, flc = fr[o], fc[o]
                cf = ((B[C][:, flc] - B[a, flc]) * (A[i, fl] - A[k, fl])).sum(1) + \
                     ((B[flc][:, C].T - B[flc, a]) * (A[fl, i] - A[fl, k])).sum(1)
                Rr = np.delete(R, ii)
                # remaining free rows Rr, cols C minus c
                if self.Asym:
                    x = A[i, Rr] - A[k, Rr]
                    Y = (B[C][:, C] - B[a, C][None, :]) + (B[C][:, C].T - B[C, a][None, :])  # row c: y over C
                else:
                    x = (A[i, Rr] - A[k, Rr]) + (A[Rr, i] - A[Rr, k])
                    Y = B[C][:, C] - B[a, C][None, :]
                xs = np.sort(x)
                # remove column c from each row of Y: set to +inf-ish then sort descending, drop the first
                Yd = Y.astype(float); Yd[np.arange(m), np.arange(m)] = -np.inf
                Ys = -np.sort(-Yd, 1)[:, :m - 1]
                lap = (Ys * xs[None, :]).sum(1) if m > 1 else np.zeros(m)
                tot = D0 + cf + lap
                v = tot.min()
                if best is None or v > best: best, arg = v, (int(k), int(i))
        return best, arg

    def onefixed_general(self, F):
        """same as onefixed, by one LAP per (k, i, c) (for instances with neither matrix symmetric)"""
        n = self.n; best, arg = None, None
        cols = [c for c in range(n) if c not in F.values()]
        for k in F:
            for i in range(n):
                if i in F: continue
                v = min(pair_min_general(self.A, self.B, {**F, i: c}, i, k) for c in cols)
                if best is None or v > best: best, arg = v, (k, i)
        return best, arg

def onefixed_fast(inst, F):
    """Vectorised one-fixed test for rank-one instances (A or B symmetric). Same value as Inst.onefixed:
    max over (k fixed, i free) of min over free cols c of [min Delta_(i,k) over completions of F + {i: c}].
    Returns (value, (k, i)), value >= 1 meaning dominated."""
    A, B, n = inst.A, inst.B, inst.n
    fr = np.array(list(F.keys()), dtype=np.int64); fc = np.array(list(F.values()), dtype=np.int64)
    mask = np.ones(n, bool); mask[fr] = False; R = np.nonzero(mask)[0]
    mask = np.ones(n, bool); mask[fc] = False; C = np.nonzero(mask)[0]
    m = len(R)
    if m < 1 or len(fr) < 1: return None, None
    Af = A.astype(np.float64); Bf = B.astype(np.float64)
    ARR = Af[np.ix_(R, R)]; BCC = Bf[np.ix_(C, C)]
    eye = np.eye(m, dtype=bool)
    best, arg = -np.inf, None
    for kk in range(len(fr)):
        k, a = fr[kk], fc[kk]
        o = np.arange(len(fr)) != kk; fl, flc = fr[o], fc[o]
        # D0[c, i]
        D0 = np.outer(Bf[C, C] - Bf[a, a], Af[R, R] - Af[k, k]) + np.outer(Bf[C, a] - Bf[a, C], Af[R, k] - Af[k, R])
        cf = (Bf[np.ix_(C, flc)] - Bf[a, flc]) @ (Af[np.ix_(R, fl)] - Af[k, fl]).T + \
             (Bf[np.ix_(flc, C)].T - Bf[flc, a]) @ (Af[np.ix_(fl, R)].T - Af[fl, k]).T
        if m == 1:
            tot = D0 + cf
        else:
            if inst.Asym:
                X = ARR - Af[k, R][None, :]
                Y = (BCC - Bf[a, C][None, :]) + (BCC.T - Bf[C, a][None, :])
            else:
                X = (ARR - Af[k, R][None, :]) + (ARR.T - Af[R, k][None, :])
                Y = BCC - Bf[a, C][None, :]
            X = X.copy(); X[eye] = np.inf; Xs = np.sort(X, 1)[:, :m - 1]
            Y = Y.copy(); Y[eye] = -np.inf; Ys = -np.sort(-Y, 1)[:, :m - 1]
            tot = D0 + cf + Ys @ Xs.T
        v = tot.min(0)            # over c, per i
        j = int(np.argmax(v))
        if v[j] > best: best, arg = float(v[j]), (int(k), int(R[j]))
    return best, arg
