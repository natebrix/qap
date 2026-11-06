#!/usr/bin/env python3
"""Box (count) inequalities on the exact rows and columns of the grid ("line
regions"), the components of the k-partition pigeonhole inequalities
(README.md section 3): for a line S, a facility set Q and b >= 1,
  b sum_{i in Q} p_i - sum_{i<j in Q} W_ij <= b (b + 1) / 2,
with p, W the moments of the free facilities' membership in S. Separated by the
same add/drop local search as dnncut.C's type 8 (cp.box_starts starts per
region and b). Prints, for the lines and for the 13 thresholds, the K most
violated cuts' violation sum and maximum, next to the K most violated regional
triangles' (the pool's family at the same iterate), so that the new family's
"violation mass" can be read against a family that is known to pay (triangles)
and one that is known not to (threshold boxes, new_cuts/box_pool).

  python3 line_boxes.py SPEC Y.bin [--K 200] [--starts 40] [--qmax 24] [--bmax 7]
"""
import argparse
import itertools
import time

import numpy as np

from nlc_analyze import Node, ROWS, COLS, thresholds


def moments(nd, cells):
    m, T = nd.m, nd.T
    idx = np.array(cells)
    W = T[np.ix_(idx, np.arange(m), idx, np.arange(m))].sum(axis=(0, 2))
    p = np.diag(W).copy()
    W = W.copy(); np.fill_diagonal(W, 0.0)
    return p, W


def box_value(p, W, Q, b):
    Q = list(Q)
    return b * p[Q].sum() - W[np.ix_(Q, Q)].sum() / 2 - b * (b + 1) / 2


def box_local(p, W, m, b, qmax, inset):
    """best-improvement add/drop from the set `inset` (bool mask), sizes in [3, qmax]"""
    inS = inset.copy()
    while True:
        Q = np.flatnonzero(inS)
        gain_add = b * p - W[:, Q].sum(1); gain_add[inS] = -np.inf
        gain_drop = -(b * p - W[:, Q].sum(1)); gain_drop[~inS] = -np.inf
        ia = int(np.argmax(gain_add)); idp = int(np.argmax(gain_drop))
        ga = gain_add[ia] if len(Q) < qmax else -np.inf
        gd = gain_drop[idp] if len(Q) > 3 else -np.inf
        if max(ga, gd) <= 1e-12:
            return inS
        if ga >= gd:
            inS[ia] = True
        else:
            inS[idp] = False


def separate_boxes(p, W, m, starts, qmax, bmax, seed):
    rng = np.random.default_rng(seed)
    qmax = min(qmax, m - 1)
    order = np.argsort(-p)
    found = {}
    for b in range(1, bmax + 1):
        if b + 2 > qmax:
            break
        for st in range(starts):
            inS = np.zeros(m, bool)
            if st == 0:
                inS[order[:b + 2]] = True
            else:
                hi = min(qmax, 3 * b + 6); k0 = b + 1 + int(rng.random() * (hi - b))
                w = p * (1 - p) + 1e-3
                keys = np.log(np.maximum(rng.random(m), 1e-300)) / w
                inS[np.argsort(-keys)[:min(k0, m)]] = True
            inS = box_local(p, W, m, b, qmax, inS)
            Q = tuple(np.flatnonzero(inS))
            v = box_value(p, W, Q, b)
            if v > 1e-3:
                key = (b, Q)
                if key not in found or found[key] < v:
                    found[key] = v
    return found


def triangles(p, W, m):
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            wij = W[i, j]
            k = np.arange(j + 1, m)
            if len(k) == 0:
                continue
            wik = W[i, k]; wjk = W[j, k]; pk = p[k]
            vals = np.stack([wij + wik - wjk - p[i], wij + wjk - wik - p[j], wik + wjk - wij - pk,
                             p[i] + p[j] + pk - wij - wik - wjk - 1])
            out.extend(vals[vals > 1e-4].tolist())
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec'); ap.add_argument('ybin')
    ap.add_argument('--K', type=int, default=200)
    ap.add_argument('--starts', type=int, default=40)
    ap.add_argument('--qmax', type=int, default=24)
    ap.add_argument('--bmax', type=int, default=7)
    ap.add_argument('--label', default='')
    a = ap.parse_args()
    nd = Node(a.spec).load(a.ybin)
    m = nd.m
    print(f'# {a.label} spec {a.spec} m {m} file {a.ybin}')
    t0 = time.time()
    lines = {}
    for r in range(ROWS):
        lines[f'row={r}'] = [c for c in range(m) if nd.r[c] == r]
    for c_ in range(COLS):
        lines[f'col={c_}'] = [c for c in range(m) if nd.c[c] == c_]
    thr = thresholds(nd)
    fam = {'lines': lines, 'thresholds': thr}
    tri_all = []
    for name, regs in fam.items():
        allv = []
        print(f'## boxes on the {name}')
        for rn, cells in regs.items():
            if not 0 < len(cells) < m:
                continue
            p, W = moments(nd, cells)
            found = separate_boxes(p, W, m, a.starts, a.qmax, a.bmax, seed=hash(rn) & 0xffff)
            vals = sorted(found.values(), reverse=True)
            sizes = [len(Q) for (b, Q) in found]
            bs = [b for (b, Q) in found]
            print(f'  {rn:6s} cap {len(cells):2d}: violated {len(found):4d}  max {vals[0] if vals else 0:.4f}'
                  f'  sum(top 20) {sum(vals[:20]):.3f}  median |Q| {np.median(sizes) if sizes else 0:.0f}'
                  f'  b range {min(bs) if bs else 0}-{max(bs) if bs else 0}')
            allv.extend(vals)
            if name == 'thresholds':
                tri_all.extend(triangles(p, W, m).tolist())
        allv.sort(reverse=True)
        print(f'  => {name}: distinct violated boxes {len(allv)}, top {a.K}: sum {sum(allv[:a.K]):.3f}'
              f' max {allv[0] if allv else 0:.4f}, median of top {a.K} {np.median(allv[:a.K]) if allv else 0:.4f}')
    tri_all = np.array(sorted(tri_all, reverse=True))
    print(f'## threshold triangles at this iterate: violated {len(tri_all)}, top {a.K}: sum {tri_all[:a.K].sum():.3f}'
          f' max {tri_all[0] if len(tri_all) else 0:.4f}, top 1000 sum {tri_all[:1000].sum():.3f}')
    print(f'# done in {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
