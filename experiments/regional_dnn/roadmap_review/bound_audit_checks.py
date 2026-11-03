"""Small exact-arithmetic audits; no QAP solves, network, or repository edits."""
from fractions import Fraction as F
from itertools import combinations, permutations
import numpy as np


def rational_psd(a):
    """Exact symmetric Gaussian elimination, including zero-pivot row checks."""
    a = [r[:] for r in a]
    positive = zero = 0
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot >= 0, (k, pivot)
        if not pivot:
            assert all(a[k][j] == 0 for j in range(k, len(a)))
            zero += 1
            continue
        positive += 1
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / pivot
                a[j][i] = a[i][j]
    return positive, zero


def higher_cut_witness():
    # A full n=8 QAP-DNN point. Region has 3 locations. Facilities form two
    # groups of 4 whose regional counts are pseudo-deterministically 3/2.
    n, k = 8, 3
    p = [F(3, 8)] * n
    W = [[p[i] if i == j else F(1, 16) if i // 4 == j // 4 else F(9, 64)
          for j in range(n)] for i in range(n)]
    Y = [[F(0) for _ in range(n*n+1)] for _ in range(n*n+1)]
    Y[0][0] = F(1)
    for a in range(n):
        for i in range(n):
            r = 1 + a*n+i
            v = p[i]/k if a < k else (1-p[i])/(n-k)
            Y[0][r] = Y[r][0] = Y[r][r] = v
            for b in range(n):
                for j in range(n):
                    if a == b or i == j:
                        continue
                    s = 1+b*n+j
                    if a < k and b < k:
                        t = W[i][j]/(k*(k-1))
                    elif a >= k and b >= k:
                        t = (1-p[i]-p[j]+W[i][j])/((n-k)*(n-k-1))
                    elif a < k:
                        t = (p[i]-W[i][j])/(k*(n-k))
                    else:
                        t = (p[j]-W[i][j])/(k*(n-k))
                    Y[r][s] = t
    assert all(y >= 0 for row in Y for y in row)
    for col in range(n*n+1):
        for a in range(n):
            assert sum(Y[1+a*n+i][col] for i in range(n)) == Y[0][col]
        for i in range(n):
            assert sum(Y[1+a*n+i][col] for a in range(n)) == Y[0][col]
    ranks = rational_psd(Y)
    tri = [sum(p[t] for t in c)-sum(W[i][j] for i,j in combinations(c,2))-1
           for c in combinations(range(n),3)]
    tri += [W[i][j]+W[i][l]-W[j][l]-p[i] for i,j,l in permutations(range(n),3)]
    assert max(tri) == -F(1,16)
    violation = sum(p[:4])-sum(W[i][j] for i,j in combinations(range(4),2))-1
    assert violation == F(1,8)
    print('Exact full QAP-DNN witness: PSD positive/zero pivots',ranks,
          '; regional triangle max',max(tri),'; conditioned five-point violation',violation)


def child_contraction():
    # Generic objective/dual matrix transfer including the constant and cross
    # terms. Exhaustively check all child permutations for all fixings at n=5.
    rng = np.random.default_rng(20261103)
    n, c = 5, 4
    A = rng.integers(-9, 10, size=(n*n+1,n*n+1))
    M = A+A.T
    checks=0
    for i in range(n):
        for j in range(n):
            fr=[a for a in range(n) if a != i]
            fc=[b for b in range(n) if b != j]
            T=np.zeros((n*n+1,c*c+1),dtype=np.int64)
            T[0,0]=T[1+i*n+j,0]=1
            for a in range(c):
                for b in range(c):
                    T[1+fr[a]*n+fc[b],1+a*c+b]=1
            C=T.T@M@T
            r=1+i*n+j
            assert C[0,0] == M[0,0]+2*M[0,r]+M[r,r]
            for perm in permutations(range(c)):
                y=np.zeros(c*c+1,dtype=np.int64);y[0]=1
                for a,b in enumerate(perm):y[1+a*c+b]=1
                x=T@y
                assert x@M@x == y@C@y
                checks+=1
    print('Exact integer parent-child objective identities:',checks)


if __name__ == '__main__':
    higher_cut_witness()
    child_contraction()
