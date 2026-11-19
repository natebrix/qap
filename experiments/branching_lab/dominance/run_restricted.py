"""Step 4 on real nodes: for the BRANCHED nodes of a trace (the only ones whose drop saves a subtree),
(a) the pure pair test's minimising completion: its cost above the cutoff (how far outside the useful region the
    LAP's adversary sits), for the best pair;
(b) the restricted test (restricted.py) over the 8 best pairs by pure min and mu in {0..10}: flagged with GLB < UB?
usage: run_restricted.py <instance> <replay.tsv.gz> <UB> [max nodes]"""
import sys, gzip, itertools, numpy as np
from scipy.optimize import linear_sum_assignment
from qapdata import read_dat, cost
from dom import Inst
from restricted import restricted_test
name, path, UB = sys.argv[1], sys.argv[2], float(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
A, B = read_dat(name); inst = Inst(A, B); n = len(A)
rows = []
with gzip.open(path, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        r = dict(zip(hdr, ln.rstrip('\n').split('\t')))
        if r['outcome'] == 'B' and r['key'] != '-' and ',' in r['key']: rows.append(r)
rng = np.random.default_rng(3)
if len(rows) > cap: rows = [rows[j] for j in sorted(rng.choice(len(rows), cap, replace=False))]
print('depth pairmin argmin_cost-UB(%UB) restricted_best mu glb-UB(%UB) flagged')
nf = 0; gaps = []
for r in rows:
    F = {int(a): int(b) for a, b in (s.split(':') for s in r['key'].split(','))}
    mins, pairs = inst.pair_mins(F)
    order = np.argsort(-mins)[:8]
    i, k = pairs[order[0]]; a, b = F[i], F[k]
    R = [l for l in range(n) if l not in F]; C = [c for c in range(n) if c not in F.values()]
    Ri, Ci = np.array(R), np.array(C)
    M = np.outer(A[i, Ri] - A[k, Ri], B[a, Ci] - B[b, Ci]) + np.outer(A[Ri, i] - A[Ri, k], B[Ci, a] - B[Ci, b])
    rr, cc = linear_sum_assignment(M)
    p = np.zeros(n, int)
    for l, c in F.items(): p[l] = c
    p[Ri[rr]] = Ci[cc]
    over = cost(A, B, p) - UB
    v, pr, mu, glb = restricted_test(A, B, F, UB, [(int(pairs[j][0]), int(pairs[j][1])) for j in order])
    fl = v >= 1 and glb < UB; nf += fl; gaps.append(over / UB)
    print(f"{r['depth']:>3} {mins[order[0]]:>10} {100*over/UB:8.2f}% {v:14.0f} {mu} {100*(glb-UB)/UB:8.2f}% {int(fl)}", flush=True)
print(f'nodes {len(rows)}; flagged by the restricted test where GLB alone does not fathom: {nf}; '
      f'median cost of the pure test minimiser above UB: {100*np.median(gaps):.2f}%')
