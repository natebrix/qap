"""Step 4: more than one exchange. For a branched node and a permutation s of a set S of its fixed rows
(p'(r) = p(s(r)) for r in S), Delta = cost(p) - cost(p') is still linear in the free assignment:
  D_l(c) = sum_{r in S} A[r,l](B[p r, c] - B[p' r, c]) + A[l,r](B[c, p r] - B[c, p' r])
  const  = sum over fixed u,v of A[u,v](B[p u, p v] - B[p' u, p' v]) + sum over fixed l not in S of D_l(p l)
min over completions = const + LAP. Tested: every 3-cycle (both orientations) of the 12 fixed rows that appear in
the node's 8 best pairs, and every product of two disjoint exchanges among those pairs.
usage: cycles.py <instance> <replay.tsv.gz> [max nodes]"""
import sys, gzip, itertools, numpy as np
from scipy.optimize import linear_sum_assignment
from qapdata import read_dat
from dom import Inst
name, path = sys.argv[1], sys.argv[2]; cap = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
A, B = read_dat(name); inst = Inst(A, B); n = len(A)
def perm_min(F, s):
    fr = list(F); pp = dict(F); pq = dict(F)
    for r, t in s.items(): pq[r] = F[t]
    R = np.array([l for l in range(n) if l not in F]); C = np.array([c for c in range(n) if c not in F.values()])
    const = sum(int(A[u, v]) * int(B[pp[u], pp[v]] - B[pq[u], pq[v]]) for u in fr for v in fr)
    M = np.zeros((len(R), len(C)))
    for r in s:
        M += np.outer(A[r, R], B[pp[r], C] - B[pq[r], C]) + np.outer(A[R, r], B[C, pp[r]] - B[C, pq[r]])
    if len(R) == 0: return const
    a, b = linear_sum_assignment(M); return const + M[a, b].sum()
rows = []
with gzip.open(path, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for ln in f:
        r = dict(zip(hdr, ln.rstrip('\n').split('\t')))
        if r['outcome'] == 'B' and r['key'].count(',') >= 2: rows.append(r)
rng = np.random.default_rng(5)
if len(rows) > cap: rows = [rows[j] for j in sorted(rng.choice(len(rows), cap, replace=False))]
flag = 0; better = 0
print('depth best_pair best_3cycle best_double')
for r in rows:
    F = {int(a): int(b) for a, b in (s.split(':') for s in r['key'].split(','))}
    mins, pairs = inst.pair_mins(F); order = np.argsort(-mins)[:8]
    top = [(int(pairs[j][0]), int(pairs[j][1])) for j in order]
    S = sorted(set(x for pr in top for x in pr))[:12]
    b3 = max(perm_min(F, {x: y, y: z, z: x}) for x, y, z in itertools.permutations(S, 3) if x < y and x < z)
    dbl = [perm_min(F, {a: b, b: a, c: d, d: c}) for (a, b), (c, d) in itertools.combinations(top, 2) if len({a, b, c, d}) == 4]
    b2 = max(dbl) if dbl else float('nan')
    flag += (b3 >= 1 or b2 >= 1); better += (max(b3, b2 if dbl else -1e18) > mins.max())
    print(r['depth'], int(mins.max()), int(b3), b2, flush=True)
print(f'nodes {len(rows)}; flagged by a 3-cycle or a double exchange: {flag}; best value above the best pair: {better}')
