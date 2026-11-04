#!/usr/bin/env python3
"""l1 / hypermetric checks on the expected distance matrix E[d_ij] of a dumped
iterate: metric (triangle) violations, negative type (max eigenvalue of
P E P, P = I - J/m; every l1 metric is of negative type, so a positive
eigenvalue v gives the valid cut v' E v <= 0, sum v = 0), and pentagonal
hypermetric violations (b = (1,1,1,-1,-1)) by random + greedy search.
python3 metric_check.py SPEC YFILE"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
spec, yf = sys.argv[1], sys.argv[2]
nd = gl.Node(spec); Y = nd.load_Y(yf); Y4, _ = nd.blocks(Y)
Er, Ec = nd.expected_dist(Y4)
for nm, E in (('row', Er), ('col', Ec), ('full', Er + Ec)):
    m = len(E)
    tri = (E[:, :, None] - E[:, None, :] - E[None, :, :].transpose(0, 2, 1))  # E_ij - E_ik - E_kj
    i, j, k = np.indices(tri.shape)
    mask = (i != j) & (j != k) & (i != k)
    P = np.eye(m) - np.ones((m, m)) / m
    ev = np.linalg.eigvalsh(P @ E @ P)
    rng = np.random.default_rng(0); best = -1e9
    for _ in range(20000):
        T = rng.choice(m, 5, replace=False); b = np.array([1, 1, 1, -1, -1])
        v = b @ E[np.ix_(T, T)] @ b
        best = max(best, v)
    print(f'{nm:4s}: max metric violation {tri[mask].max():.4f}; max eig of P E P {ev.max():.4e} (min {ev.min():.1f}); '
          f'max pentagonal b\'Eb (<=0 valid) over 20k random 5-sets {best:.4f}')
