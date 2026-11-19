"""Mirror basics: optimum, its mirror (two rows exchanged), symmetry of the matrices. usage: basics.py"""
import sys, numpy as np
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, read_sln, cost
for name, (i, k) in [('tai35b', (19, 18)), ('tai50b', (31, None)), ('tai60b', (44, 43))]:
    A, B = read_dat(name); opt, p = read_sln(name)
    print(name, 'n', len(p), 'opt', opt, 'cost(sln)', cost(A, B, p), 'A sym', (A == A.T).all(), 'B sym', (B == B.T).all())
    if k is not None:
        q = p.copy(); q[[i, k]] = q[[k, i]]
        print(f'  path fixes {i}:{p[i]},{k}:{p[k]}  mirror cost - opt = {cost(A, B, q) - opt}')
    # smallest exchange gains of the optimum
    g = []
    for a in range(len(p)):
        for b in range(a + 1, len(p)):
            q = p.copy(); q[[a, b]] = q[[b, a]]; g.append((cost(A, B, q) - opt, a, b))
    g.sort(); print('  smallest exchange gains of the optimum:', g[:5])
