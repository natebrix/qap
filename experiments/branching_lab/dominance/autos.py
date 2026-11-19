"""Automorphisms of a symmetric matrix M (permutations g with M[g u, g v] = M[u, v]) by backtracking."""
import numpy as np
def automorphisms(M, limit=1000):
    n = len(M); out = []
    rowsig = [tuple(sorted(M[u])) for u in range(n)]
    def bt(g, used):
        u = len(g)
        if u == n: out.append(list(g)); return len(out) >= limit
        for v in range(n):
            if v in used or rowsig[v] != rowsig[u] or M[v, v] != M[u, u]: continue
            if all(M[g[w], v] == M[w, u] for w in range(u)):
                g.append(v); used.add(v)
                if bt(g, used): return True
                g.pop(); used.discard(v)
        return False
    bt([], set()); return out
def optimal_orbit(A, B, p):
    """known optimum p and its images under Aut(A) (p o g) and Aut(B) (h o p)"""
    GA = automorphisms(A) if (A == A.T).all() else [list(range(len(A)))]
    GB = automorphisms(B) if (B == B.T).all() else [list(range(len(A)))]
    S = set()
    for g in GA:
        for h in GB:
            S.add(tuple(np.array(h)[p[np.array(g)]]))
    return [np.array(s) for s in S], len(GA), len(GB)
