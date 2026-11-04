#!/usr/bin/env python3
"""Exact minimum total pairwise Manhattan distance of k distinct cells of a
grid (optionally restricted to an allowed cell set), the right-hand side of
the clique dispersion inequality

    sum_{i<j in T} E_Y[d(loc_i, loc_j)] >= Dmin(k)     (|T| = k).

Method: for a cell set with row counts R_r and column counts C_c,
sum_{pairs} d = sum_{r<r'} R_r R_r' (r'-r) + sum_{c<c'} C_c C_c' (c'-c),
so Dmin(k) = min f(R) + g(C) over margins (R, C) realisable by a 0/1 matrix
supported on the allowed cells (checked by max-flow). Margins are enumerated
in increasing f and g; exact.

  python3 dispersion.py [--rows 5 --cols 10] [--kmax 20] [--forbid a,b,...]
"""
import argparse, heapq, itertools
import numpy as np


def pair_sum(counts):
    s = 0
    for a in range(len(counts)):
        for b in range(a + 1, len(counts)):
            s += counts[a] * counts[b] * (b - a)
    return s


def compositions(k, parts, cap):
    """all vectors of `parts` nonnegative ints <= cap[i] summing to k"""
    def rec(i, left):
        if i == parts - 1:
            if left <= cap[i]:
                yield (left,)
            return
        for v in range(min(left, cap[i]) + 1):
            for rest in rec(i + 1, left - v):
                yield (v,) + rest
    yield from rec(0, k)


def feasible(R, C, allowed):
    """0/1 matrix with row sums R, column sums C, support in allowed (rows x cols bool)"""
    rows, cols = allowed.shape
    # max flow: source -> row r (cap R_r), row r -> col c (cap 1 if allowed), col c -> sink (cap C_c)
    # simple augmenting paths (small graph)
    flow = np.zeros((rows, cols), dtype=int)
    rf = np.zeros(rows, dtype=int); cf = np.zeros(cols, dtype=int)
    total = sum(R)
    for _ in range(total):
        # BFS from rows with residual source capacity
        prev_r = {}; prev_c = {}
        frontier = [r for r in range(rows) if rf[r] < R[r]]
        for r in frontier:
            prev_r[r] = None
        found = None
        q = list(frontier)
        while q and found is None:
            nq = []
            for r in q:
                for c in range(cols):
                    if allowed[r, c] and flow[r, c] == 0 and c not in prev_c:
                        prev_c[c] = r
                        if cf[c] < C[c]:
                            found = c; break
                        for r2 in range(rows):
                            if flow[r2, c] == 1 and r2 not in prev_r:
                                prev_r[r2] = c; nq.append(r2)
                if found is not None:
                    break
            q = nq
        if found is None:
            return False
        c = found; cf[c] += 1
        while True:
            r = prev_c[c]; flow[r, c] = 1
            pc = prev_r[r]
            if pc is None:
                rf[r] += 1; break
            flow[r, pc] = 0; c = pc
    return True


def dmin(k, allowed):
    rows, cols = allowed.shape
    rcap = allowed.sum(axis=1); ccap = allowed.sum(axis=0)
    if k > allowed.sum():
        return None
    Rs = sorted((pair_sum(R), R) for R in compositions(k, rows, rcap))
    Cs = sorted((pair_sum(C), C) for C in compositions(k, cols, ccap))
    # best-first over (i, j) pairs
    heap = [(Rs[0][0] + Cs[0][0], 0, 0)]; seen = {(0, 0)}
    while heap:
        v, i, j = heapq.heappop(heap)
        if feasible(Rs[i][1], Cs[j][1], allowed):
            return v, Rs[i][1], Cs[j][1]
        for i2, j2 in ((i + 1, j), (i, j + 1)):
            if i2 < len(Rs) and j2 < len(Cs) and (i2, j2) not in seen:
                seen.add((i2, j2)); heapq.heappush(heap, (Rs[i2][0] + Cs[j2][0], i2, j2))
    return None


def brute(k, allowed):
    cells = [(r, c) for r in range(allowed.shape[0]) for c in range(allowed.shape[1]) if allowed[r, c]]
    best = None
    for S in itertools.combinations(cells, k):
        s = sum(abs(a[0] - b[0]) + abs(a[1] - b[1]) for a, b in itertools.combinations(S, 2))
        best = s if best is None or s < best else best
    return best


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--rows', type=int, default=5); ap.add_argument('--cols', type=int, default=10)
    ap.add_argument('--kmax', type=int, default=20); ap.add_argument('--forbid', default='')
    ap.add_argument('--check', action='store_true', help='compare with brute force on small cases')
    a = ap.parse_args()
    allowed = np.ones((a.rows, a.cols), dtype=bool)
    for t in filter(None, a.forbid.split(',')):
        allowed[int(t) // a.cols, int(t) % a.cols] = False
    if a.check:
        rng = np.random.default_rng(1)
        for trial in range(40):
            r, c = rng.integers(2, 5), rng.integers(2, 6)
            al = rng.random((r, c)) > 0.25
            for k in range(2, min(6, al.sum()) + 1):
                v = dmin(k, al)[0]; b = brute(k, al)
                assert v == b, (r, c, al, k, v, b)
        print('check: dmin == brute force on 40 random masked grids, k = 2..6')
    for k in range(2, a.kmax + 1):
        v, R, C = dmin(k, allowed)
        print(k, v, 'rows', R, 'cols', C)
