#!/usr/bin/env python3
"""wil50's flow matrix B (second matrix): distribution, spectrum, dense
subsets, heavy triangles; compared with a random symmetric 0-9 matrix."""
import sys, itertools
import numpy as np
from qapio import read_qaplib
n, A, B = read_qaplib('/home/user/qap/qaplib/wil50.dat')
iu = np.triu_indices(n, 1); w = B[iu]
print('flows: mean %.3f sd %.3f, zeros %d of %d pairs; row sums %d..%d (mean %.1f, sd %.1f)' % (
    w.mean(), w.std(), (w == 0).sum(), len(w), B.sum(1).min(), B.sum(1).max(), B.sum(1).mean(), B.sum(1).std()))
ev = np.linalg.eigvalsh(B)
rng = np.random.default_rng(1)
R = np.triu(rng.integers(0, 10, (n, n)), 1).astype(float); R = R + R.T
evr = np.linalg.eigvalsh(R)
print('spectrum B: top 6', np.round(ev[::-1][:6], 1), 'bottom 3', np.round(ev[:3], 1))
print('spectrum of a random 0-9 matrix: top 6', np.round(evr[::-1][:6], 1), 'bottom 3', np.round(evr[:3], 1))
# centered: B - mean
Bc = B - w.mean() * (1 - np.eye(n))
evc = np.linalg.eigvalsh(Bc); Rc = R - R[iu].mean() * (1 - np.eye(n)); evrc = np.linalg.eigvalsh(Rc)
print('centered B: top 5', np.round(evc[::-1][:5], 1), ' random centered top 5', np.round(evrc[::-1][:5], 1),
      ' semicircle edge 2*sd*sqrt(n) = %.1f' % (2 * w.std() * np.sqrt(n)))
# correlation of rows: is there block structure?
C = np.corrcoef(B + np.diag(B.sum(1) / (n - 1)))
print('row-correlation off-diagonal: mean %.3f, max %.3f, min %.3f' % (C[iu].mean(), C[iu].max(), C[iu].min()))
# heavy triangles
for th in (9, 8, 7):
    H = B >= th
    tri = int(np.trace(np.linalg.matrix_power(H.astype(int), 3)) // 6)
    pr = (w >= th).mean()
    print('triangles with all flows >= %d: %d (random expectation %.0f)' % (th, tri, pr ** 3 * n * (n - 1) * (n - 2) / 6))
# densest k-subsets by greedy peeling + local search, vs random
def dens(S): S = list(S); return B[np.ix_(S, S)].sum() / 2
def best_k(k, M, tries=200, seed=0):
    r = np.random.default_rng(seed); best = None
    for t in range(tries):
        S = set(r.choice(n, k, replace=False).tolist())
        improved = True
        while improved:
            improved = False
            for i in list(S):
                for j in range(n):
                    if j in S: continue
                    d = M[j, list(S - {i})].sum() - M[i, list(S - {i})].sum()
                    if d > 1e-9: S.remove(i); S.add(j); improved = True; break
                if improved: break
        v = M[np.ix_(list(S), list(S))].sum() / 2
        if best is None or v > best[0]: best = (v, sorted(S))
    return best
for k in (4, 5, 6, 8, 10, 15):
    vb, Sb = best_k(k, B, tries=60)
    vr, _ = best_k(k, R, tries=60)
    print('densest %2d-subset: B %.0f (mean flow %.2f) %s | random matrix %.0f (%.2f)' % (
        k, vb, vb / (k * (k - 1) / 2), Sb, vr, vr / (k * (k - 1) / 2)))
