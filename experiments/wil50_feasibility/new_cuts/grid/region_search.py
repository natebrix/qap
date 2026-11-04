#!/usr/bin/env python3
"""Free-region separation of Boolean triangles: is any cell region much more
violated than the 13 prefix regions (and the parity regions)?

For a triple (i, j, k) and a triangle type, the violation is a quadratic
function of the region's 0/1 indicator s over the free cells:
   p_i(s) = sum_a s_a Y_(ai),(ai),  W_ij(s) = s' Y_(.i),(.j) s.
Alternating search: fix the region, take the most violated triples
(vectorized over all triples); fix a triple, improve the region by 1-flip
local search; restart from every named region. Any region is valid (z_i^S is
0/1 for every S), so whatever is found is a valid cut.

  python3 region_search.py SPEC YFILE [--top 30]
"""
import argparse, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from analyze_Y import tri_viol, region_pW, regions_of


def viol_of_region(Y4, s, i, j, k, typ):
    p = np.array([np.dot(s, Y4[:, t, :, t].diagonal()) for t in (i, j, k)])
    Wij = s @ Y4[:, i, :, j] @ s; Wik = s @ Y4[:, i, :, k] @ s; Wjk = s @ Y4[:, j, :, k] @ s
    if typ == 3:
        return p.sum() - Wij - Wik - Wjk - 1
    if typ == 0:   # apex i: W_ij + W_ik - W_jk - p_i  (delta_jk - delta_ij - delta_ik)/2
        return Wij + Wik - Wjk - p[0]
    if typ == 1:   # apex j
        return Wij + Wjk - Wik - p[1]
    return Wik + Wjk - Wij - p[2]


def improve_region(Y4, s, i, j, k, typ, max_pass=20):
    s = s.copy(); best = viol_of_region(Y4, s, i, j, k, typ)
    for _ in range(max_pass):
        improved = False
        for a in np.random.permutation(len(s)):
            s[a] = 1 - s[a]
            if 0 < s.sum() < len(s):
                v = viol_of_region(Y4, s, i, j, k, typ)
                if v > best + 1e-9:
                    best = v; improved = True; continue
            s[a] = 1 - s[a]
        if not improved:
            break
    return best, s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y'); ap.add_argument('--top', type=int, default=20)
    a = ap.parse_args()
    np.random.seed(1)
    nd = gl.Node(a.spec); Y = nd.load_Y(a.Y); Y4, _ = nd.blocks(Y)
    regs = regions_of(nd)
    t0 = time.time()
    print('region  start-max  -> best after free-region search (top triples)')
    allbest = []
    for name in ['row<2', 'row<3', 'col<4', 'col<5', 'col<6', 'checker', 'rowpar', 'colpar']:
        if name not in regs:
            continue
        S = regs[name]; p, W = region_pW(nd, Y4, S)
        V, (I, J, K) = tri_viol(p, W)
        flat = np.argsort(V.ravel())[::-1][:a.top]
        s0 = np.zeros(nd.m); s0[S] = 1
        gains = []
        for f in flat:
            t, typ = divmod(f, 4)
            v0 = V.ravel()[f]
            v1, s1 = improve_region(Y4, s0, I[t], J[t], K[t], typ)
            gains.append((v0, v1, int(s1.sum())))
            allbest.append((v1, name, I[t], J[t], K[t], typ, s1))
        g = np.array([(x, y) for x, y, _ in gains])
        print(f'{name:8s} start max {g[:, 0].max():.4f} mean(top{a.top}) {g[:, 0].mean():.4f} -> searched max {g[:, 1].max():.4f} mean {g[:, 1].mean():.4f}  sizes {sorted(set(z for _, _, z in gains))[:6]}')
    allbest.sort(key=lambda z: -z[0])
    v, name, i, j, k, typ, s = allbest[0]
    cells = [nd.free_loc[c] for c in np.where(s > 0)[0]]
    grid = np.zeros((nd.rows, nd.cols), dtype=int)
    for c in cells:
        grid[c // nd.cols, c % nd.cols] = 1
    print(f'best found {v:.4f} from {name}, triple {i},{j},{k} type {typ}; region:')
    print(grid)
    print(f'{time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
