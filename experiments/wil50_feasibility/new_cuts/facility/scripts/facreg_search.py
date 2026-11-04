#!/usr/bin/env python3
"""Most violated facility-region triangles: for a location triple (a, b, c)
and a facility set T (indicator t), with q_a = sum_{j in T} x_aj and
U_ab = sum_{j,l in T} Y_(aj),(bl) (the regional moments of w^T, the indicator
of 'location a holds a facility of T'), the four Boolean triangles are valid
for every permutation; for fixed (a, b, c) each is a quadratic in t, so the
best T is a 47-variable QUBO, searched here by multi-start 1-flip local search.
  python3 facreg_search.py SPEC Y.bin [--triples N]"""
import argparse, itertools, time
import numpy as np
from ynode import Node

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y')
    ap.add_argument('--triples', type=int, default=200); ap.add_argument('--starts', type=int, default=20)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args(); t0 = time.time()
    nd = Node(a.spec).load(a.Y); m = nd.m; x = nd.x; T = nd.T
    rng = np.random.default_rng(a.seed)
    A = nd.A
    trip = list(itertools.combinations(range(m), 3))
    # half the sample: compact triples (pairwise distance <= 2), half random
    comp = [t for t in trip if A[t[0], t[1]] <= 2 and A[t[0], t[2]] <= 2 and A[t[1], t[2]] <= 2]
    sel = [comp[i] for i in rng.choice(len(comp), min(len(comp), a.triples // 2), replace=False)]
    sel += [trip[i] for i in rng.choice(len(trip), a.triples - len(sel), replace=False)]
    best_all = []
    for (p, q, r) in sel:
        # Ypq[j, l] = Y_(pj),(ql)
        Ypq, Ypr, Yqr = T[p, :, q, :], T[p, :, r, :], T[q, :, r, :]
        # type 3: sum_j t_j (x_pj + x_qj + x_rj) - t'(sym(Ypq + Ypr + Yqr)) t - 1   (U_ab counts ordered pairs j,l in T)
        lin3 = x[p] + x[q] + x[r]; Q3 = Ypq + Ypr + Yqr
        # apex p: U_pq + U_pr - U_qr - q_p ; apex q: U_pq + U_qr - U_pr - q_q ; apex r: U_pr + U_qr - U_pq - q_r
        forms = [(lin3, Q3, -1.0), (-x[p], -(Ypq + Ypr - Yqr), 0.0), (-x[q], -(Ypq + Yqr - Ypr), 0.0), (-x[r], -(Ypr + Yqr - Ypq), 0.0)]
        for typ, (lin, Q, c0) in enumerate(forms):
            Qs = (Q + Q.T) / 2
            # violation(t) = lin.t - t'Qs t + c0   (U_ab with j, l over T x T, both orders: t' Q t)
            def f(t): return lin @ t - t @ Q @ t + c0
            best = -9
            for s in range(a.starts):
                t = (rng.random(m) < rng.uniform(0.1, 0.6)).astype(float)
                cur = f(t); improved = True
                while improved:
                    improved = False
                    # gain of flipping j: d = lin_j (1-2t_j) - [(1-2t_j)(2 (Qs t)_j) + Q_jj ... ] computed exactly
                    for j in rng.permutation(m):
                        t[j] = 1 - t[j]; v = f(t)
                        if v > cur + 1e-12: cur = v; improved = True
                        else: t[j] = 1 - t[j]
                best = max(best, cur)
            best_all.append((best, typ, (p, q, r), A[p, q] + A[p, r] + A[q, r]))
    best_all.sort(key=lambda z: -z[0])
    v = np.array([b[0] for b in best_all])
    print('facility-region triangles, %d location triples x 4 types (QUBO over T): max %.4f; > 0.01: %d; > 0.05: %d; time %.0f s' % (
        len(sel), v.max(), (v > 0.01).sum(), (v > 0.05).sum(), time.time() - t0))
    for b in best_all[:8]:
        print('   viol %.4f type %d locations %s (perimeter %d)' % b)

if __name__ == '__main__':
    main()
