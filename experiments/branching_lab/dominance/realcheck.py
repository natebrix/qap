"""Checks on tai60b's own data (n = 60, B asymmetric) for flagged nodes of the replay:
(1) the sort (rank-one) value equals the general LAP value (scipy, both rank-one terms kept separately);
(2) the LAP's minimising completion p, evaluated directly: cost(p) - cost(swap(p)) equals the value;
(3) 200 random completions and 20 local-search (2-exchange descent inside the node) completions all have
    Delta >= the value (so >= 1).  usage: realcheck.py <replay.tsv.gz> <count>"""
import sys, gzip, numpy as np
from scipy.optimize import linear_sum_assignment
from qapdata import read_dat, cost
from dom import Inst, pair_min_general
A, B = read_dat('tai60b'); inst = Inst(A, B); n = 60
rng = np.random.default_rng(11)
fl = []
with gzip.open(sys.argv[1], 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        r = dict(zip(hdr, ln.rstrip('\n').split('\t')))
        if r['pairmin'] != '' and int(r['pairmin']) >= 1: fl.append(r)
B_ = [r for r in fl if r['outcome'] == 'B']; others = [r for r in fl if r['outcome'] != 'B']
pick = B_ + [others[j] for j in rng.choice(len(others), min(len(others), int(sys.argv[2])), replace=False)]
def delta(p, i, k):
    q = p.copy(); q[[i, k]] = q[[k, i]]; return cost(A, B, p) - cost(A, B, q)
bad = [0, 0, 0]; minseen = None
for r in pick:
    F = {int(a): int(b) for a, b in (s.split(':') for s in r['key'].split(','))}
    i, k = map(int, r['pair'].split(',')); v = int(r['pairmin'])
    g = pair_min_general(A, B, F, i, k); bad[0] += (g != v)
    R = np.array([l for l in range(n) if l not in F]); C = np.array([c for c in range(n) if c not in F.values()])
    a, b = F[i], F[k]
    M = np.outer(A[i, R] - A[k, R], B[a, C] - B[b, C]) + np.outer(A[R, i] - A[R, k], B[C, a] - B[C, b])
    p = np.zeros(n, int)
    for l, c in F.items(): p[l] = c
    if len(R):
        rr, cc = linear_sum_assignment(M); p[R[rr]] = C[cc]
    bad[1] += (delta(p, i, k) != v)
    for t in range(220):
        p = np.zeros(n, int)
        for l, c in F.items(): p[l] = c
        p[R] = rng.permutation(C)
        if t >= 200 and len(R) > 1:   # descent on cost by exchanges of free rows
            cur = cost(A, B, p); improved = True
            while improved:
                improved = False
                for u in range(len(R)):
                    for w in range(u + 1, len(R)):
                        q = p.copy(); q[[R[u], R[w]]] = q[[R[w], R[u]]]; cq = cost(A, B, q)
                        if cq < cur: p, cur, improved = q, cq, True
        dl = delta(p, i, k); bad[2] += (dl < v)
        minseen = dl if minseen is None else min(minseen, dl - v)
print(f'flagged nodes checked {len(pick)} ({len(B_)} branched); general LAP != sort value: {bad[0]}; '
      f'minimiser Delta != value: {bad[1]}; sampled completions with Delta < value: {bad[2]} (smallest Delta - value seen: {minseen})')
