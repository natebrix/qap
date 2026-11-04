"""Validity check of cp_small's box rows at the nug12 root: separate the top
200 at the dnntest iterate (wil50/plumb/nug12.Y.bin), evaluate every cut on
100,000 random permutations and on permutations with every count value."""
import sys, numpy as np, itertools
sys.path.insert(0, '.')
from relax_small import read_qaplib, grid_regions
from cp_small import box_separate, box_row
n, A, B, C = read_qaplib('inst/nug12_sw.dat'); N = n * n + 1
Y = np.fromfile('../wil50/plumb/nug12.Y.bin').reshape(N, N); Y = (Y + Y.T) / 2
regions = grid_regions(n, 3, 4)
cuts = box_separate(Y, n, [S for _, S in regions], 200)
rng = np.random.default_rng(3)
perms = np.array([rng.permutation(n) for _ in range(100000)])   # perms[t, a] = side-1 index at side-0 a
worst_all = -np.inf; ntight = 0
sizes = []
for v, (kind, P, qid, b) in cuts:
    coefs, h = box_row(n, list(P), qid, b)
    Q = [j for j in range(n) if (qid >> j) & 1]
    sizes.append((len(P), len(Q), b))
    # value on permutation: sum over pairs of lifted indices present
    keys = np.array(list(coefs.keys())); vals = np.array(list(coefs.values()))
    e1 = keys // N; e2 = keys % N
    a1 = (e1 - 1) // n; j1 = (e1 - 1) % n; a2 = (e2 - 1) // n; j2 = (e2 - 1) % n
    hit = (perms[:, a1] == j1) & (perms[:, a2] == j2)
    val = hit.astype(float) @ vals
    w = val.max() - h
    worst_all = max(worst_all, w); ntight += abs(w) < 1e-9
print("cuts", len(cuts), "violation range %.4f..%.4f" % (cuts[-1][0], cuts[0][0]))
print("worst (value - rhs) over 100k permutations: %.3g; cuts tight on some permutation: %d" % (worst_all, ntight))
from collections import Counter
print("(|P|,|Q|,b) most common:", Counter(sizes).most_common(8))
