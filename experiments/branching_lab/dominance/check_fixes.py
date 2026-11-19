import sys, numpy as np
from qapdata import *
for name in ['tai60b', 'tai50b', 'tai35b', 'nug30', 'tho40', 'tai10b', 'nug12']:
    A, B = read_dat(name); val, p = read_sln(name)
    q = np.argsort(p)
    print(name, 'n', len(p), 'sln', val, 'cost(p)', cost(A, B, p), 'cost(p^-1)', cost(A, B, q),
          'A sym', (A == A.T).all(), 'B sym', (B == B.T).all(), 'diagA0', (np.diag(A) == 0).all(), 'diagB0', (np.diag(B) == 0).all())
