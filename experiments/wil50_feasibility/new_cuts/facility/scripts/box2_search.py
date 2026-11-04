#!/usr/bin/env python3
"""Jünger-Kaibel box inequalities with BOTH sides free: for a location set P
and a facility set Q the count c = x(P x Q) = |pi(Q) cap P| is an integer, so
Var(c) >= f (1 - f), f = frac(E c), with E c = s'Xt, E c^2 = s_a s_b t_j t_l
Y_(aj),(bl) (s, t indicators of P, Q). Alternating local search: for fixed P
the problem in Q is a QUBO over the regional moments of P, and vice versa.
Reports violation in triangle units (half the box measure).
  python3 box2_search.py SPEC Y.bin [--starts N]"""
import argparse, time
import numpy as np
from ynode import Node

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--starts', type=int, default=40); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args(); t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m; X = nd.x; T = nd.T
    rng = np.random.default_rng(a.seed)
    def val(s, t):
        mu = s @ X @ t
        e2 = np.einsum('a,ajbl,b->jl', s, T, s, optimize=True)
        e2 = t @ e2 @ t
        f = mu - np.floor(mu)
        return 0.5 * (f * (1 - f) - (e2 - mu * mu))
    def qubo_t(s, t):   # optimise t for fixed s by 1-flips
        W = np.einsum('a,ajbl,b->jl', s, T, s, optimize=True); p = s @ X
        def g(t):
            mu = p @ t; f = mu - np.floor(mu); return 0.5 * (f * (1 - f) - (t @ W @ t - mu * mu))
        cur = g(t); imp = True
        while imp:
            imp = False
            for j in rng.permutation(m):
                t[j] = 1 - t[j]; v = g(t)
                if v > cur + 1e-12: cur, imp = v, True
                else: t[j] = 1 - t[j]
        return t, cur
    def qubo_s(s, t):
        W = np.einsum('j,ajbl,l->ab', t, T, t, optimize=True); p = X @ t
        def g(s):
            mu = p @ s; f = mu - np.floor(mu); return 0.5 * (f * (1 - f) - (s @ W @ s - mu * mu))
        cur = g(s); imp = True
        while imp:
            imp = False
            for b in rng.permutation(m):
                s[b] = 1 - s[b]; v = g(s)
                if v > cur + 1e-12: cur, imp = v, True
                else: s[b] = 1 - s[b]
        return s, cur
    res = []
    thr = [nd.rc[:, 0] < r for r in range(1, 5)] + [nd.rc[:, 1] < c for c in range(1, 10)]
    for st in range(a.starts):
        if st < len(thr): s = thr[st].astype(float)        # start from a threshold region
        else: s = (rng.random(m) < rng.uniform(0.1, 0.6)).astype(float)
        t = (rng.random(m) < rng.uniform(0.05, 0.5)).astype(float)
        prev = -9
        for it in range(10):
            t, v = qubo_t(s, t)
            s, v = qubo_s(s, t)
            if v <= prev + 1e-9: break
            prev = v
        isthr = any(np.array_equal(s, z.astype(float)) for z in thr)
        res.append((v, int(s.sum()), int(t.sum()), st < len(thr), isthr))
    res.sort(key=lambda z: -z[0])
    print('box inequalities, P and Q free (alternating QUBO search, %d starts): best violation (triangle units) %.4f' % (a.starts, res[0][0]))
    for r in res[:10]:
        print('   viol %.4f  |P| %2d  |Q| %2d  started from a threshold %s  P still a threshold %s' % r)
    print('time %.0f s' % (time.time() - t0))

if __name__ == '__main__':
    main()
