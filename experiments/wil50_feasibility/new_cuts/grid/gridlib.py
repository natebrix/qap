#!/usr/bin/env python3
"""Shared helpers for the grid/ experiments (wil50: first matrix = Manhattan
distance on a 5 x 10 grid, cells row-major, 10 per row).

Conventions (as in dnncuts/dnntest and stronger_bounds/local7):
  spec "i:j,..." on qaplib/wil50.dat: i = location (grid cell, A index),
  j = facility (flow, B index). Free locations / facilities keep their
  original order; the lifted index of (local location a, local facility i)
  is 1 + a*m + i.
"""
import itertools
import numpy as np

ROOT = '/home/user/qap'
ROWS, COLS = 5, 10


def read_qaplib(path):
    tok = open(path).read().split()
    n = int(tok[0]); v = list(map(int, tok[1:1 + 2 * n * n]))
    A = np.array(v[:n * n], dtype=float).reshape(n, n)
    B = np.array(v[n * n:2 * n * n], dtype=float).reshape(n, n)
    return n, A, B


def write_qaplib(path, A, B):
    n = len(A)
    with open(path, 'w') as f:
        f.write(f'{n}\n\n')
        for M in (A, B):
            for r in M:
                f.write(' '.join(str(int(round(x))) for x in r) + '\n')
            f.write('\n')


def cell_rc(a):
    return a // COLS, a % COLS


def grid_dist(rows=ROWS, cols=COLS):
    n = rows * cols
    r = np.arange(n) // cols; c = np.arange(n) % cols
    Dr = np.abs(r[:, None] - r[None, :]).astype(float)
    Dc = np.abs(c[:, None] - c[None, :]).astype(float)
    return Dr, Dc


def parse_spec(spec):
    if spec in ('-', ''):
        return []
    return [tuple(map(int, t.split(':'))) for t in spec.split(',')]


class Node:
    """A node of wil50 (or any grid instance): free locations, free facilities,
    and the lifted objective pieces split into row and column parts."""

    def __init__(self, spec, A=None, B=None, rows=ROWS, cols=COLS):
        if A is None:
            n, A, B = read_qaplib(ROOT + '/qaplib/wil50.dat')
        self.n = n = len(A); self.A = A; self.B = B
        self.rows, self.cols = rows, cols
        self.fixes = parse_spec(spec)
        fl = {x for x, _ in self.fixes}; ff = {y for _, y in self.fixes}
        self.free_loc = [a for a in range(n) if a not in fl]
        self.free_fac = [i for i in range(n) if i not in ff]
        self.m = m = len(self.free_loc)
        self.N = m * m + 1
        Dr, Dc = grid_dist(rows, cols)
        self.Dr_full, self.Dc_full = Dr, Dc
        L = np.array(self.free_loc); Fc = np.array(self.free_fac)
        self.Dr = Dr[np.ix_(L, L)]; self.Dc = Dc[np.ix_(L, L)]
        self.F = B[np.ix_(Fc, Fc)]
        # linear terms from fixed facilities: for free location a and free facility i,
        # 2 * sum_fixed (A[a, x] B[i, y])  (A, B symmetric)
        self.lin_r = np.zeros((m, m)); self.lin_c = np.zeros((m, m))
        for x, y in self.fixes:
            self.lin_r += 2 * np.outer(Dr[L, x], B[Fc, y])
            self.lin_c += 2 * np.outer(Dc[L, x], B[Fc, y])
        const_r = const_c = 0.0
        for (x1, y1), (x2, y2) in itertools.product(self.fixes, repeat=2):
            const_r += Dr[x1, x2] * B[y1, y2]; const_c += Dc[x1, x2] * B[y1, y2]
        self.const_r, self.const_c = const_r, const_c
        self.loc_r = np.array([a // cols for a in self.free_loc])
        self.loc_c = np.array([a % cols for a in self.free_loc])

    def idx(self, a, i):
        return 1 + a * self.m + i

    def perm_cost(self, perm, part='both'):
        """perm[i] = local location of local facility i"""
        perm = np.asarray(perm)
        tot = 0.0
        for D, lin, c0, nm in ((self.Dr, self.lin_r, self.const_r, 'row'), (self.Dc, self.lin_c, self.const_c, 'col')):
            if part not in ('both', nm):
                continue
            tot += (self.F * D[np.ix_(perm, perm)]).sum() + lin[perm, np.arange(self.m)].sum() + c0
        return tot

    def load_Y(self, path):
        Y = np.fromfile(path)
        assert Y.size == self.N * self.N, (Y.size, self.N)
        return Y.reshape(self.N, self.N)

    def blocks(self, Y):
        """Y4[a, i, b, j] = Y_(ai),(bj) and x[a, i] = Y_0,(ai)"""
        m = self.m
        Y4 = Y[1:, 1:].reshape(m, m, m, m)
        x = Y[0, 1:].reshape(m, m)
        return Y4, x

    def expected_dist(self, Y4):
        """E[d_row(i,j)], E[d_col(i,j)] over facility pairs (m x m)"""
        Er = np.einsum('ab,aibj->ij', self.Dr, Y4, optimize=True)
        Ec = np.einsum('ab,aibj->ij', self.Dc, Y4, optimize=True)
        return Er, Ec

    def objective_parts(self, Y):
        Y4, x = self.blocks(Y)
        Er, Ec = self.expected_dist(Y4)
        dg = np.einsum('aiai->ai', Y4)
        r = (self.F * Er).sum() + (self.lin_r * dg).sum() + self.const_r
        c = (self.F * Ec).sum() + (self.lin_c * dg).sum() + self.const_c
        return r, c


def region_moments(Y4, S):
    """p[i] = sum_{a in S} Y_(ai),(ai), W[i,j] = sum_{a,b in S} Y_(ai),(bj)"""
    S = np.asarray(S)
    sub = Y4[S][:, :, S, :]
    W = sub.sum(axis=(0, 2))
    p = np.einsum('aiai->i', Y4[np.ix_(S, np.arange(Y4.shape[1]), S, np.arange(Y4.shape[3]))]) if False else \
        np.array([sum(Y4[a, i, a, i] for a in S) for i in range(Y4.shape[1])])
    return p, W
