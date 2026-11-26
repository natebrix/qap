#!/usr/bin/env python3
"""cop5.py Y.bin [--cap N]: is a DNN iterate (dnntest/dnncuts -dumpY: the (1+m^2)^2 doubles of Y, index 0 the
constant, 1 + i*m + j the assignment i -> j) locally completely positive?

A copositive matrix K supported on an index set S gives the cut <K, Y_S> >= 0, valid because the exact QAP is a
linear program over the completely positive cone (Burer; Povh-Rendl). Such a cut can be violated only if the
principal submatrix Y_S is doubly nonnegative but not completely positive, which needs |S| >= 5.

1. Scale: Yh = D^-1/2 Y D^-1/2, D = diag(Y) (cones are invariant; Yh has unit diagonal, Yh_0a = sqrt(x_a)).
2. Diagonal dominance: a nonnegative symmetric matrix whose comparison matrix is PSD is completely positive
   (Drew, Johnson, Loewy 1994), so every Y_S with Yh_S diagonally dominant admits no copositive cut at all. Reported:
   the largest k such that every principal submatrix of order k is diagonally dominant, without and with index 0.
3. Horn cuts: for 5-cycles (a1..a5) of the graph {Yh_ab >= tau}, the least value of d'(H o Yh_S)d over the simplex
   (H the Horn matrix in that cyclic order, o the entrywise product). Negative = a violated cut D H D.
"""
import sys, itertools, numpy as np

def cop_min(M):
    """M: (B,5,5) symmetric. min over the simplex of d'Md, by the stationary points of every support."""
    B = M.shape[0]; best = M[:, np.arange(5), np.arange(5)].min(1).copy()
    for r in range(2, 6):
        for I in itertools.combinations(range(5), r):
            sub = M[:, I][:, :, I]
            ok = np.abs(np.linalg.det(sub)) > 1e-12
            u = np.zeros((B, r)); u[ok] = np.linalg.solve(sub[ok], np.ones((ok.sum(), r, 1)))[..., 0]
            s = u.sum(1); pos = ok & ((u > 0).all(1) | (u < 0).all(1)) & (np.abs(s) > 1e-12)
            val = np.where(pos, 1.0 / np.where(pos, s, 1.0), np.inf)
            best = np.minimum(best, val)
    return best

HORN = np.array([[1, -1, 1, 1, -1], [-1, 1, -1, 1, 1], [1, -1, 1, -1, 1], [1, 1, -1, 1, -1], [-1, 1, 1, -1, 1]], float)

def cycles5(adj, cap):
    """5-cycles of a graph (list of neighbour sets), each once (least vertex first, second < last), at most cap"""
    out = []
    for a in range(len(adj)):
        na = [v for v in adj[a] if v > a]
        for b in na:
            for c in adj[b]:
                if c <= a or c == b: continue
                for d in adj[c]:
                    if d <= a or d == b or d == c: continue
                    for e in adj[d]:
                        if e <= b or e == c or e == d or a not in adj[e] or e == a: continue
                        out.append((a, b, c, d, e))
                        if len(out) >= cap: return out
    return out

def main():
    path = sys.argv[1]; cap = int(sys.argv[sys.argv.index('--cap') + 1]) if '--cap' in sys.argv else 300000
    Y = np.fromfile(path); N = int(round(len(Y) ** .5)); Y = Y.reshape(N, N); Y = (Y + Y.T) / 2
    m = int(round((N - 1) ** .5)); d = np.diag(Y).copy()
    print(f"{path}: N={N} m={m} Y00={Y[0,0]:.6f} min entry {Y.min():.2e} least eigenvalue {np.linalg.eigvalsh(Y)[0]:.2e}")
    x = Y[0, 1:]; print(f"x: max {x.max():.4f} median {np.median(x):.4f} entries above 0.5: {(x > .5).sum()}, above 0.2: {(x > .2).sum()}")
    live = np.where(d > 1e-7)[0]; Yh = np.maximum(Y[np.ix_(live, live)], 0); s = np.sqrt(d[live]); Yh = Yh / np.outer(s, s)
    n = len(live); off = Yh.copy(); np.fill_diagonal(off, 0)
    a0 = off[1:, 1:]; pos = a0[a0 > 0]
    print(f"scaled off-diagonal entries between assignments: nonzero {len(pos)//2}, median {np.median(pos):.4f}, 99% {np.quantile(pos,.99):.4f}, 99.99% {np.quantile(pos,.9999):.4f}, max {pos.max():.4f}")
    print(f"scaled entries to the constant (sqrt x): median {np.median(off[0,1:]):.4f}, max {off[0,1:].max():.4f}")
    for name, M in (("without the constant", off[1:, 1:]), ("with the constant", off)):
        srt = -np.sort(-M, axis=1); cs = np.cumsum(srt, axis=1)       # cs[:, k-2] = sum of the k-1 largest in the row
        k = 1 + int((cs <= 1).all(0).sum())
        print(f"diagonal dominance ({name}): every principal submatrix of order <= {k} is diagonally dominant, hence completely positive")
    for tau in (0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1):
        adj = [set(np.nonzero(off[a] >= tau)[0].tolist()) for a in range(n)]
        ne = sum(len(v) for v in adj) // 2
        if ne > 60000: print(f"tau {tau}: {ne} edges, skipped"); break
        cyc = cycles5(adj, cap)
        if not cyc: print(f"tau {tau}: {ne} edges, no 5-cycle"); continue
        C = np.array(cyc); best = np.inf; arg = None
        for lo in range(0, len(C), 20000):
            S = C[lo:lo + 20000]; M = Yh[S[:, :, None], S[:, None, :]] * HORN
            v = cop_min(M); i = v.argmin()
            if v[i] < best: best, arg = v[i], S[i]
        lab = ['const' if live[t] == 0 else f"{(live[t]-1)//m}:{(live[t]-1)%m}" for t in arg]
        print(f"tau {tau}: {ne} edges, {len(C)} 5-cycles{' (cap)' if len(C) >= cap else ''}, least Horn value {best:+.4f} at {' '.join(lab)}", flush=True)

main()
