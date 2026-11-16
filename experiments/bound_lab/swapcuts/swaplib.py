#!/usr/bin/env python3
"""Exchange-optimality cuts: the derivation, shared by the stage scripts.

Cost convention (the solver's, qap.C): C(p) = sum_{k,l} A[k,l] B[p(k),p(l)],
p maps an index of the first matrix (a "facility") to an index of the second
(a "location"); `-x i:a` fixes p(i) = a.

Exchange of i (at a) and j (at b): p'(i) = b, p'(j) = a, the rest unchanged.
  delta = C(p') - C(p)
        = sum_{k not in {i,j}} g_k(p(k))
          + (A[i,j] - A[j,i]) (B[b,a] - B[a,b])
          + (A[i,i] - A[j,j]) (B[b,b] - B[a,a]),
  g_k(c) = (A[i,k] - A[j,k]) (B[b,c] - B[a,c]) + (A[k,i] - A[k,j]) (B[c,b] - B[c,a]).
With i fixed at a by the node and x_jb = 1 this is linear in the other free
assignments: delta = c0 + sum_{free k != j, free c != b} g_k(c) x_kc, c0 = the two
last terms plus sum over the other fixed facilities k of g_k(p(k)).
"""
import numpy as np


def read_dat(path):
    t = open(path).read().split()
    n = int(t[0])
    v = np.array(t[1:1 + 2 * n * n], float)
    return n, v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)


def cost(A, B, p):
    p = np.asarray(p)
    return float((A * B[np.ix_(p, p)]).sum())


def swap_terms(A, B, i, a, j, b):
    """c_pair (the constant of the pair itself) and G[k, c] = g_k(c) for every k, c
    (rows i and j are meaningless and set to 0)."""
    n = A.shape[0]
    dr = B[b, :] - B[a, :]            # B[b,c] - B[a,c]
    dc = B[:, b] - B[:, a]            # B[c,b] - B[c,a]
    G = np.outer(A[i, :] - A[j, :], dr) + np.outer(A[:, i] - A[:, j], dc)
    G[i, :] = 0.0
    G[j, :] = 0.0
    cpair = (A[i, j] - A[j, i]) * (B[b, a] - B[a, b]) + (A[i, i] - A[j, j]) * (B[b, b] - B[a, a])
    return cpair, G


def swap_delta(A, B, p, i, j):
    """C(p with i and j exchanged) - C(p) by the formula."""
    a, b = p[i], p[j]
    cpair, G = swap_terms(A, B, i, a, j, b)
    n = len(p)
    return cpair + sum(G[k, p[k]] for k in range(n) if k != i and k != j)


def node_cut(A, B, fixed, i, j, b):
    """The cut for fixed facility i (at fixed[i]), free j, free location b:
    returns (c0, g) with g a dict-free dense array over ALL (k, c) that is zero
    outside free k != j, free c != b.  Valid: x_jb * (c0 + <g, x>) >= 0."""
    n = A.shape[0]
    a = fixed[i]
    cpair, G = swap_terms(A, B, i, a, j, b)
    c0 = cpair + sum(G[k, c] for k, c in fixed.items() if k != i)
    g = G.copy()
    used = set(fixed.values())
    for k in fixed:
        g[k, :] = 0.0
    for c in used:
        g[:, c] = 0.0
    g[j, :] = 0.0
    g[:, b] = 0.0
    return c0, g


def fixed_pair_cut(A, B, fixed, i, i2):
    """Two fixed facilities exchanged: l(x) = c0 + <g, x> >= 0, linear in x."""
    a, b = fixed[i], fixed[i2]
    cpair, G = swap_terms(A, B, i, a, i2, b)
    c0 = cpair + sum(G[k, c] for k, c in fixed.items() if k != i and k != i2)
    g = G.copy()
    for k in fixed:
        g[k, :] = 0.0
    for c in set(fixed.values()):
        g[:, c] = 0.0
    return c0, g
