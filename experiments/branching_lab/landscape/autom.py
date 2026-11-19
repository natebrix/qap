"""autom.py: automorphisms of a symmetric matrix (permutations s with M[s][:, s] == M), by backtracking with
row-multiset refinement; used for the starting points of the grid instances (every image of the best known
solution under the symmetries of either matrix is also optimal)."""
import numpy as np

def automorphisms(M, limit=64):
    n = len(M)
    key = [tuple(sorted(M[i])) for i in range(n)]
    cand = [[j for j in range(n) if key[j] == key[i]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i]))
    out = []; s = [-1] * n; used = [False] * n

    def rec(k):
        if len(out) >= limit: return
        if k == n:
            out.append(np.array(s)); return
        i = order[k]
        for j in cand[i]:
            if used[j]: continue
            ok = True
            for kk in range(k):
                ii = order[kk]
                if M[i, ii] != M[j, s[ii]] or M[ii, i] != M[s[ii], j]: ok = False; break
            if ok:
                s[i] = j; used[j] = True; rec(k + 1); used[j] = False; s[i] = -1
    rec(0)
    return out
