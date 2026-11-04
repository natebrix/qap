#!/usr/bin/env python3
"""Iterated 2-exchange local search for a wil50 node (ynode conventions):
pi[a] = local facility on local location a. Prints the best value and writes
the permutation (original indices, full: location -> facility) to --out."""
import argparse, time
import numpy as np
from ynode import Node

def value(nd, pi):
    return float((nd.A * nd.B[np.ix_(pi, pi)]).sum() + nd.C[np.arange(nd.m), pi].sum() + nd.shift)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('--seconds', type=float, default=60)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--out')
    a = ap.parse_args()
    nd = Node(a.spec); m = nd.m; A, B, C = nd.A, nd.B, nd.C
    rng = np.random.default_rng(a.seed)
    def ls(pi):
        cur = value(nd, pi)
        improved = True
        while improved:
            improved = False
            for r in rng.permutation(m):
                pr = pi[r]; ps = pi                      # candidate partners s = 0..m-1
                Bpi = B[:, pi]                           # Bpi[f, k] = B[f, pi_k]
                full = ((A[r][None, :] - A) * (Bpi[ps] - Bpi[pr][None, :])).sum(1)
                d = 2 * (full + 2 * A[r] * B[pr, ps]) + C[r, ps] + C[np.arange(m), pr] - C[r, pr] - C[np.arange(m), ps]
                d[r] = 0
                t = int(np.argmin(d))
                if d[t] < -1e-9:
                    pi[r], pi[t] = pi[t], pi[r]; cur += d[t]; improved = True
        return pi, value(nd, pi)
    t0 = time.time(); best = None
    pi = rng.permutation(m)
    pi, cur = ls(pi)
    best = (cur, pi.copy())
    while time.time() - t0 < a.seconds:
        q = best[1].copy()
        for _ in range(3):
            i, j = rng.choice(m, 2, replace=False); q[i], q[j] = q[j], q[i]
        q, v = ls(q)
        if v < best[0]: best = (v, q.copy()); print('%.1f s: %.0f' % (time.time() - t0, v), flush=True)
    v, pi = best
    print('best %.0f' % v)
    if a.out:
        full = {}
        for la, lf in enumerate(pi): full[nd.loc[la]] = nd.fac[lf]
        for x, y in nd.fixes: full[x] = y
        open(a.out, 'w').write(' '.join(str(full[x]) for x in range(nd.n)) + '\n')

if __name__ == '__main__':
    main()
