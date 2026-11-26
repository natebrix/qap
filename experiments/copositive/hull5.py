#!/usr/bin/env python3
"""hull5.py Y.bin [TOP]: how far is a dumped iterate from the exact hull on small sets of assignments?

Every copositive cut on a principal submatrix Y_S of 0/1 variables is valid for conv{yy'} over the feasible 0/1
vectors y on S, so the distance of Y_S from that hull bounds what any such cut (Horn, Burer-Dong, Hildebrand's
T-matrices, ...) can be violated by. Here:
 1. atomic triangles, all triples (a, b, c) of assignments: Y_ab + Y_ac - Y_bc - x_a (and the fourth type
    x_a + x_b + x_c - Y_ab - Y_ac - Y_bc - 1), the largest violations;
 2. 5-sets grown from the most violated triangles (two more indices, each the one that raises the hull distance
    most among the 40 most correlated): the hull distance t = min over convex combinations of feasible yy' of the
    largest entry difference (an LP), against the triangle violation it started from."""
import sys, itertools, numpy as np
from scipy.optimize import linprog

path = sys.argv[1]; TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 30
Y = np.fromfile(path); N = int(round(len(Y) ** .5)); Y = Y.reshape(N, N); Y = np.maximum((Y + Y.T) / 2, 0); m = int(round((N - 1) ** .5))
A = Y[1:, 1:].astype(np.float32); x = np.diag(A).copy(); n = len(x)
fi, lo = np.divmod(np.arange(n), m)
conf = ((fi[:, None] == fi[None, :]) | (lo[:, None] == lo[None, :])) & ~np.eye(n, dtype=bool)
def lab(S): return ' '.join(f"{t//m}:{t%m}" for t in S)

best = []                                              # (violation, a, b, c)
A0 = A.copy(); np.fill_diagonal(A0, 0)
for a in range(n):
    if x[a] < 1e-6: continue
    r = A0[a]; T = r[:, None] + r[None, :] - A0 - x[a]   # b, c; needs b != c != a
    T[a, :] = -9; T[:, a] = -9; np.fill_diagonal(T, -9)
    k = int(T.argmax()); best.append((float(T.flat[k]), a, k // n, k % n))
best.sort(reverse=True)
print(f"{path}: m={m}; largest atomic triangle violations Y_ab + Y_ac - Y_bc - x_a (one per a):")
for v, a, b, c in best[:5]: print(f"  {v:+.4f}  a={lab([a])} b={lab([b])} c={lab([c])}  x_a={x[a]:.3f} Y_ab={A[a,b]:.3f} Y_ac={A[a,c]:.3f} Y_bc={A[b,c]:.3f}")

def hull_dist(S):
    S = list(S); k = len(S); V = [y for y in itertools.product((0, 1), repeat=k)
                                 if not any(y[i] and y[j] and conf[S[i], S[j]] for i in range(k) for j in range(i))]
    iu = [(i, j) for i in range(k) for j in range(i, k)]
    M = np.array([[y[i] * y[j] for (i, j) in iu] for y in V], float).T      # entries x vertices
    tgt = np.array([A[S[i], S[j]] for (i, j) in iu], float)
    nv = len(V); c = np.zeros(nv + 1); c[-1] = 1
    Aub = np.block([[M, -np.ones((len(iu), 1))], [-M, -np.ones((len(iu), 1))]]); bub = np.concatenate([tgt, -tgt])
    Aeq = np.concatenate([np.ones(nv), [0]])[None]; r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1], bounds=[(0, None)] * (nv + 1), method='highs')
    return r.fun

Ah = A0 / np.sqrt(np.outer(x, x) + 1e-30)
print(f"5-sets grown from the {TOP} most violated triangles: triangle violation -> hull distance of the triple -> of the 5-set")
out = []
for v, a, b, c in best[:TOP]:
    S = [a, b, c]; t3 = hull_dist(S)
    for _ in range(2):
        cand = np.argsort(-(Ah[S].sum(0)))[:60]; cand = [int(u) for u in cand if u not in S][:40]
        vals = [hull_dist(S + [u]) for u in cand]; S.append(cand[int(np.argmax(vals))])
    t5 = hull_dist(S); out.append((t5, t3, v, S))
out.sort(reverse=True)
for t5, t3, v, S in out[:8]: print(f"  triangle {v:+.4f}  triple {t3:.4f}  5-set {t5:.4f}  {lab(S)}")
print(f"  largest 5-set hull distance {out[0][0]:.4f}; median {np.median([o[0] for o in out]):.4f}; largest triangle violation {best[0][0]:+.4f}")
