"""Magnitudes at tai35b's twin node 19:12,18:26 (pair 19,18): the lifted objective's max |entry| (quadratic part
(a_ik b_jl + a_ki b_lj)/2 over free rows/cols, linear part C/2), the range of D_l(c), the min and max of Delta over
the node's completions (rank-one LAP, A symmetric), Delta at the optimum's mirror. usage: scale.py [INSTANCE FIXES (first two = the pair)] > scale_INSTANCE.out"""
import sys, numpy as np
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, read_sln, cost
from scipy.optimize import linear_sum_assignment
inst = sys.argv[1] if len(sys.argv) > 1 else 'tai35b'
A, B = read_dat(inst); opt, p = read_sln(inst); n = len(A)
A = A.astype(float); B = B.astype(float); Bs = (B + B.T) / 2
fx = dict(tuple(map(int, t.split(':'))) for t in (sys.argv[2] if len(sys.argv) > 2 else '19:12,18:26').split(','))
i, k = list(fx)[:2]; a, b = fx[i], fx[k]
fr = [r for r in range(n) if r not in fx]; fc = [c for c in range(n) if c not in fx.values()]
Aq = A[np.ix_(fr, fr)]; Bq = Bs[np.ix_(fc, fc)]
print('max |A_free| %g  max |B_free (sym)| %g  max quadratic lifted entry ~ %g' % (abs(Aq).max(), abs(Bq).max(), abs(Aq).max() * abs(Bq).max()))
Cl = np.array([[sum(A[l, f] * B[c, fx[f]] + A[f, l] * B[fx[f], c] for f in fx) for c in fc] for l in fr])
print('reduced linear term C: max |C| %g, so max |C/2| %g' % (abs(Cl).max(), abs(Cl).max() / 2))
D = np.array([[(A[i, l] - A[k, l]) * (B[a, c] - B[b, c]) + (A[l, i] - A[l, k]) * (B[c, a] - B[c, b]) for c in fc] for l in fr])
D0 = (A[i, i] - A[k, k]) * (B[a, a] - B[b, b]) + (A[i, k] - A[k, i]) * (B[a, b] - B[b, a])
print('D: max |D_l(c)| %g, D0 %g' % (abs(D).max(), D0))
r, c = linear_sum_assignment(D); dmin = D0 + D[r, c].sum()
r, c = linear_sum_assignment(-D); dmax = D0 + D[r, c].sum()
print('Delta over the twin node: min %g  max %g' % (dmin, dmax))
q = p.copy(); q[[i, k]] = q[[k, i]]
print('the optimum mirrored into the twin: cost - opt %d, Delta %d' % (cost(A.astype(np.int64), B.astype(np.int64), q) - opt,
      cost(A.astype(np.int64), B.astype(np.int64), q) - cost(A.astype(np.int64), B.astype(np.int64), p)))
print('|A[i,l]-A[k,l]| over free l, largest:', sorted(abs(A[i, fr] - A[k, fr]).astype(int))[-8:], ' A[i,k] =', A[i, k])
