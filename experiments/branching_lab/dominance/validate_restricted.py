"""Brute force check of restricted.py: whenever the restricted test flags (N, pair, mu), every completion of N with
cost < UB has Delta >= 1 for that pair. Instances with A symmetric (tai10b, tai12b subtrees, nug08, random).
UB is set at opt + slack for several slacks, so that N may hold solutions below UB. Output: validate_restricted.out"""
import itertools, numpy as np
from qapdata import read_dat, read_sln
from restricted import restricted_test
from validate import perms_of, all_costs
rng = np.random.default_rng(7)
tot = dict(nodes=0, flagged=0, flagged_beyond_glb=0, fail=0, nodes_with_sub_ub=0)
for nm in ['tai10b', 'nug08', 'tai10a', 'had10', 'rand_sym8', 'rand_symA_asymB8']:
    if nm.startswith('rand'):
        n = 8; A = rng.integers(0, 4, (n, n)); A = np.triu(A, 1); A = A + A.T
        B = rng.integers(0, 4, (n, n)); 
        if nm == 'rand_sym8': B = B + B.T
    else:
        A, B = read_dat(nm); n = len(A)
    P = perms_of(None, None, {}, n); costs = all_costs(A, B, P); opt = int(costs.min())
    srt = np.sort(costs)
    for t in range(120):
        UB = float(srt[[0, 5, 50, 500][t % 4]] + 0.5)
        d = int(rng.integers(2, n - 1)); p = P[np.argsort(costs)[rng.integers(0, 200)]]
        rows = rng.choice(n, d, replace=False)
        if t % 3 == 0:  # random columns
            F = {int(r): int(c) for r, c in zip(rows, rng.permutation(n)[:d])}
        else:
            F = {int(r): int(p[r]) for r in rows}
        sel = np.all(P[:, list(F)] == np.array(list(F.values()), dtype=np.int8), 1)
        PN, cN = P[sel], costs[sel]
        pairs = [(i, k) for i, k in itertools.combinations(F, 2)]
        v, pr, mu, glb = restricted_test(A, B, F, UB, pairs)
        tot['nodes'] += 1; tot['nodes_with_sub_ub'] += int((cN < UB).any())
        if v >= 1:
            tot['flagged'] += 1; tot['flagged_beyond_glb'] += int(glb < UB)
            i, k = pr; Q = PN.copy(); Q[:, [i, k]] = Q[:, [k, i]]
            dl = cN - all_costs(A, B, Q)
            if ((cN < UB) & (dl < 1)).any(): tot['fail'] += 1; print('FAIL', nm, F, pr, mu, v)
    print(nm, tot, flush=True)
print('TOTAL', tot)
