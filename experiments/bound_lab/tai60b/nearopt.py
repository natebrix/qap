#!/usr/bin/env python3
"""nearopt.py [INST]: how close do the off-path children come? For a list of on-path prefixes (fixes of the best
known permutation p, as the search made them) and a branching row r, every child r:v (v != p[r]) gets an upper
bound on its subtree optimum: 2-swap descent from p with v moved to r (r and the prefix held fixed).
Convention of the node traces: entry r of the permutation = v; cost = sum A[i,k] B[p(i),p(k)], A the first matrix.
Prints, per child, (subtree upper bound - best known); a child can only fathom when the bound's slack at that
depth is below this excess. With argument 'swaps': the distribution of all 2-swap excesses of p."""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/regions')
import gen_regions as g
inst = 'tai60b'
Q = __file__.rsplit('/', 4)[0]
n, A, B = g.read_qaplib(f'{Q}/qaplib/{inst}.dat')
t = open(f'{Q}/QAPSOLN/{inst}.sln').read().split()
best, p0 = int(t[1]), np.array(list(map(int, t[2:2 + n]))) - 1
def cost(p): return int((A * B[np.ix_(p, p)]).sum())
assert cost(p0) == best
def descent(p, fixed):
    p = p.copy(); c = cost(p); free = [i for i in range(n) if i not in fixed]
    improved = True
    while improved:
        improved = False
        for a in range(len(free)):
            for b in range(a + 1, len(free)):
                i, j = free[a], free[b]
                q = p.copy(); q[i], q[j] = q[j], q[i]
                cq = cost(q)
                if cq < c: p, c, improved = q, cq, True
    return p, c
def child(fixed, r, v):
    p = p0.copy(); j = int(np.where(p == v)[0][0]); p[j], p[r] = p[r], p[j]
    return descent(p, set(fixed) | {r})[1] - best
if sys.argv[1] == 'swaps':
    ex = sorted(cost(np.array([p0[j] if k == i else p0[i] if k == j else p0[k] for k in range(n)])) - best
                for i in range(n) for j in range(i + 1, n))
    print('2-swap excess over best known: min', ex[0], 'count below 0.25M/0.5M/1M/2M/3M/5M:',
          [sum(e < x for e in ex) for x in (0.25e6, 0.5e6, 1e6, 2e6, 3e6, 5e6)], 'of', len(ex))
else:
    # nearopt.py ROW fixedrow,fixedrow,... : all children of the on-path node with those rows fixed
    r = int(sys.argv[1]); fixed = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 and sys.argv[2] else []
    out = sorted((child(fixed, r, v), v) for v in range(n) if v != p0[r] and v not in [p0[f] for f in fixed])
    print(f'row {r} (on-path {r}:{p0[r]}), fixed rows {fixed}: excess of the best completion found, per child')
    print(' '.join(f'{r}:{v}={e/1e6:.2f}M' for e, v in out))
    print('children with excess below 0.5M/1M/2M/3M/5M:', [sum(e < x for e, v in out) for x in (0.5e6, 1e6, 2e6, 3e6, 5e6)], 'of', len(out))
