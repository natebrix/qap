#!/usr/bin/env python3
"""Exact check of the K_{3,3} example (Proposition "not implied") in dnn_plus.tex.

n = 6, A = B = adjacency of K_{3,3}. Builds the lifted matrix Y of the proof in rationals and
checks: feasibility for the DNN relaxation (Y_00 = 1, diagonal = row 0, gangster zeros,
0 <= Y <= 1, row/column sum equations, PSD), its objective (9), the optimum by enumeration (10),
the triangle inequalities on Y's own entries, and the violated regional inequality (1/8).
PSD is exact: the characteristic polynomial is factored over the rationals (sympy).
"""
from fractions import Fraction as F
from itertools import permutations, combinations
import sympy as sp

n = 6
grp = lambda t: t // 3                      # groups {0,1,2}, {3,4,5} for facilities and locations
A = [[1 if grp(i) != grp(j) else 0 for j in range(n)] for i in range(n)]
idx = lambda i, a: 1 + i * n + a            # facility i at location a
N = n * n + 1
W = lambda i, j: F(1, 8) if grp(i) == grp(j) else F(1, 4)   # (1 + Q_ij) / 4, i != j

Y = [[F(0)] * N for _ in range(N)]
Y[0][0] = F(1)
for i in range(n):
    for a in range(n):
        u = idx(i, a)
        Y[0][u] = Y[u][0] = Y[u][u] = F(1, 6)
        for j in range(n):
            for b in range(n):
                if i != j and a != b:
                    Y[u][idx(j, b)] = W(i, j) / 6 if grp(a) == grp(b) else (F(1, 2) - W(i, j)) / 9

# feasibility
assert all(Y[u][v] == Y[v][u] and 0 <= Y[u][v] <= 1 for u in range(N) for v in range(N))
for i in range(n):
    for a in range(n):
        for b in range(n):
            if a != b:
                assert Y[idx(i, a)][idx(i, b)] == 0 and Y[idx(a, i)][idx(b, i)] == 0   # gangster
for v in range(N):                           # sum_a Y_(ia),v = Y_0,v and sum_i Y_(ia),v = Y_0,v
    for t in range(n):
        assert sum(Y[idx(t, a)][v] for a in range(n)) == Y[0][v]
        assert sum(Y[idx(i, t)][v] for i in range(n)) == Y[0][v]
M = sp.Matrix(N, N, lambda u, v: sp.Rational(Y[u][v].numerator, Y[u][v].denominator))
x = sp.symbols('x')
roots = sp.roots(sp.Poly(M.charpoly(x).as_expr(), x))
assert sum(roots.values()) == N and all(r >= 0 for r in roots), roots
print("eigenvalues (value: multiplicity):", dict(sorted(roots.items(), key=lambda kv: -kv[0])))

# objective and optimum
obj = sum(A[i][j] * A[a][b] * Y[idx(i, a)][idx(j, b)] for i in range(n) for j in range(n)
          for a in range(n) for b in range(n))
opt = min(sum(A[i][j] * A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in permutations(range(n)))
print("objective of Y:", obj, " optimum:", opt)
assert obj == 9 and opt == 10

# triangle inequalities on Y's own entries (lifted variables with distinct facilities and locations)
worst = None
lifted = [(i, a) for i in range(n) for a in range(n)]
for (u, v, w) in combinations(lifted, 3):
    if len({u[0], v[0], w[0]}) < 3 or len({u[1], v[1], w[1]}) < 3:
        continue
    a_, b_, c_ = idx(*u), idx(*v), idx(*w)
    ab, ac, bc = Y[a_][b_], Y[a_][c_], Y[b_][c_]
    for apex, (p, q, r) in ((a_, (ab, ac, bc)), (b_, (ab, bc, ac)), (c_, (ac, bc, ab))):
        viol = p + q - r - Y[apex][apex]
        worst = viol if worst is None else max(worst, viol)
    viol = Y[a_][a_] + Y[b_][b_] + Y[c_][c_] - ab - ac - bc - 1
    worst = max(worst, viol)
print("largest violation of a triangle inequality on Y's entries:", worst)
assert worst < 0

# the regional inequality: S = first group of locations, facilities 0, 1, 2
S = [0, 1, 2]
p = [sum(Y[idx(i, a)][idx(i, a)] for a in S) for i in range(3)]
Wr = lambda i, j: sum(Y[idx(i, a)][idx(j, b)] for a in S for b in S)
viol = sum(p) - Wr(0, 1) - Wr(0, 2) - Wr(1, 2) - 1
print("violation of the regional inequality:", viol)
assert viol == F(1, 8)

# objective = 10 + 8 g on permutations, g the inequality's slack
for q in permutations(range(n)):
    z = [1 if q[i] in S else 0 for i in range(3)]
    g = 1 - sum(z) + z[0] * z[1] + z[0] * z[2] + z[1] * z[2]
    assert sum(A[i][j] * A[q[i]][q[j]] for i in range(n) for j in range(n)) == 10 + 8 * g
print("objective = 10 + 8 g on all 720 permutations: OK")
