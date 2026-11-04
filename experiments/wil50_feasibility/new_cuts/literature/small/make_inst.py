#!/usr/bin/env python3
"""wil50-like small instances: Manhattan distance on a rows x cols grid (first
matrix, side 0 = locations, as in wil50) and wil50's flow matrix restricted to
a random subset of facilities (second matrix)."""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from relax_small import read_qaplib, write_qaplib
rows, cols, seed, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
n50, A50, B50, _ = read_qaplib('/home/user/qap/qaplib/wil50.dat')
n = rows * cols
rng = np.random.default_rng(seed)
S = np.sort(rng.choice(50, n, replace=False))
A = np.array([[abs(a // cols - b // cols) + abs(a % cols - b % cols) for b in range(n)] for a in range(n)])
B = B50[np.ix_(S, S)]
write_qaplib(out, A, B)
print(out, "facilities", S.tolist())
