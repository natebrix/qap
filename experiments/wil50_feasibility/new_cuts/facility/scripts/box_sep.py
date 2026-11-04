#!/usr/bin/env python3
"""Box (count-integrality) inequalities of any support on a region's 0/1
vector z (Juenger-Kaibel boxes, Padberg clique inequalities with switchings):
for signs s in {-1, 0, 1}^m, c = sum_{s_i = 1} z_i + sum_{s_i = -1} (1 - z_i)
is an integer, so E[(c - b)(c - b - 1)] >= 0 for every integer b, i.e.
Var(c) >= f (1 - f), f = frac(E c), in the moments (p, W) of z:
  violation(s) = f (1 - f) - s' (W - p p') s     (W_ii = p_i).
|support| = 3 are the triangles, 5 the 5-point count inequalities. Heuristic
separation by local search over s. Applied to the 13 threshold regions
(location side, s over facilities) and to facility regions (s over locations).
  python3 box_sep.py SPEC Y.bin [--starts N]"""
import argparse, itertools, time
import numpy as np
from ynode import Node

def viol(s, p, Sig):
    mu = s @ p; f = mu - np.floor(mu)
    return f * (1 - f) - s @ Sig @ s

def search(p, W, rng, starts=200, maxsupp=None):
    m = len(p); Wd = W.copy(); np.fill_diagonal(Wd, p)
    Sig = Wd - np.outer(p, p)
    best = (-1, None)
    for st in range(starts):
        k = rng.integers(3, 12)
        s = np.zeros(m); idx = rng.choice(m, k, replace=False); s[idx] = rng.choice([-1, 1], k)
        cur = viol(s, p, Sig)
        improved = True
        while improved:
            improved = False
            for i in rng.permutation(m):
                for v in (-1.0, 0.0, 1.0):
                    if v == s[i]: continue
                    old = s[i]; s[i] = v
                    if maxsupp and np.count_nonzero(s) > maxsupp: s[i] = old; continue
                    nv = viol(s, p, Sig)
                    if nv > cur + 1e-12: cur = nv; improved = True
                    else: s[i] = old
        if cur > best[0]: best = (cur, s.copy())
    return best

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--starts', type=int, default=150); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args(); t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m; x = nd.x; T = nd.T
    rng = np.random.default_rng(a.seed)
    regs = [('row<%d' % t, nd.rc[:, 0] < t) for t in range(1, 5)] + [('col<%d' % t, nd.rc[:, 1] < t) for t in range(1, 10)]
    regs += [('checker', nd.rc.sum(1) % 2 == 0), ('rowpar', nd.rc[:, 0] % 2 == 0), ('colpar', nd.rc[:, 1] % 2 == 0)]
    print('location regions (s over facilities): region capacity | best box violation, support size; best with support <= 3 (triangle), <= 5')
    for name, mask in regs:
        idx = np.where(mask)[0]
        p = x[idx].sum(0); W = T[np.ix_(idx, range(m), idx, range(m))].sum(axis=(0, 2))
        v, s = search(p, W, rng, a.starts)
        v3, _ = search(p, W, rng, a.starts // 3, 3); v5, _ = search(p, W, rng, a.starts // 3, 5)
        print('  %-8s %2d | %.4f (support %2d)  tri %.4f  five %.4f' % (name, len(idx), v, np.count_nonzero(s), v3, v5))
    print('facility regions (random sets; s over locations):')
    for k in (5, 10, 20, 23):
        for rep in range(2):
            Tset = np.sort(rng.choice(m, k, replace=False))
            q = x[:, Tset].sum(1); U = T[np.ix_(range(m), Tset, range(m), Tset)].sum(axis=(1, 3))
            v, s = search(q, U, rng, a.starts // 2)
            print('  |T| = %2d  best %.4f (support %d)' % (k, v, np.count_nonzero(s)))
    print('time %.1f s' % (time.time() - t0))

if __name__ == '__main__':
    main()
