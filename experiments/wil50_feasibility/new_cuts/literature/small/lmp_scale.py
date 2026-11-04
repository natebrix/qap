#!/usr/bin/env python3
"""How the DNN's 3-local inconsistency scales with n: on a dumped DNN iterate
(dnntest -dumpY, root of a QAPLIB instance, (n^2+1)^2 doubles), the exact
LMP3 separation value (potentials in [-1,1], as cp_small.lmp_separate) for
triples of side-0 indices and of side-1 indices: the heaviest triples by the
matrix of that side, and random triples.
  python3 -I lmp_scale.py inst.dat PREFIX [--heavy 30] [--rand 30]"""
import argparse, itertools, json, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from relax_small import read_qaplib
from cp_small import lmp_separate

ap = argparse.ArgumentParser(); ap.add_argument('dat'); ap.add_argument('prefix')
ap.add_argument('--heavy', type=int, default=30); ap.add_argument('--rand', type=int, default=30)
ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
n, A, B, C = read_qaplib(a.dat)
N = n * n + 1
Y = np.fromfile(a.prefix + '.Y.bin').reshape(N, N); Y = (Y + Y.T) / 2
rng = np.random.default_rng(a.seed)
res = {}
for side, M in ((0, A), (1, B)):
    T = list(itertools.combinations(range(n), 3))
    w = np.array([abs(M[i, j]) + abs(M[i, k]) + abs(M[j, k]) for i, j, k in T])
    heavy = [T[t] for t in np.argsort(-w)[:a.heavy]]
    rand = [T[t] for t in rng.choice(len(T), min(a.rand, len(T)), replace=False)]
    for kind, lst in (('heavy', heavy), ('rand', rand)):
        out = lmp_separate(Y, n, side, len(lst), triples_list=lst)
        vals = [v for v, _ in out] + [0.0] * (len(lst) - len(out))
        res[f'side{side}_{kind}'] = (float(max(vals)), float(np.median(vals)))
        print(f"n={n} side {side} {kind:5s}: max {max(vals):.4f} median {np.median(vals):.4f}", flush=True)
print("RESULT", json.dumps(res))
