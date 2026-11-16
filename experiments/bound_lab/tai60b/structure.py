#!/usr/bin/env python3
"""structure.py INSTANCE.dat SLN: check the best known value from the .sln permutation and print the
cluster structure of the locations (first matrix) at several dendrogram levels, with the MDS coordinates
of the cluster centres and, per cluster, which facilities the best known solution puts there."""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/regions')
import gen_regions as g
n, D, F = g.read_qaplib(sys.argv[1])
t = open(sys.argv[2]).read().split()
val, p = int(t[1]), np.array(list(map(int, t[2:2 + n]))) - 1
v1 = int((D * F[np.ix_(p, p)]).sum())          # sum A[i,k] B[p(i),p(k)], A first matrix
pinv = np.argsort(p)
v2 = int((D[np.ix_(p, p)] * F).sum())
print(f'{sys.argv[2].split("/")[-1]}: listed {val}; sum A[i,k]B[p(i),p(k)] = {v1}; sum A[p(i),p(k)]B[i,k] = {v2}')
X, w, stress = g.mds(D, 2)
print('MDS eig', w[:3].round(0), 'stress', round(stress, 4))
for K in map(int, sys.argv[3:]):
    cls = g.top_clusters(D, K)
    print(f'K={K}')
    for s in sorted(cls, key=lambda s: X[s, 0].mean()):
        s = sorted(s)
        diam = D[np.ix_(s, s)].max()
        print(f'  size {len(s):2d} centre ({X[s,0].mean():7.0f},{X[s,1].mean():7.0f}) diam {diam:5.0f}  locs {s}')
