#!/usr/bin/env python3
"""Structure of the open instances' data, for the literature notes (README.md, Q2c).
No solver is run.  usage (from the repo root): python3 -I experiments/bound_lab/literature/structure.py

Prints, per instance, the sparsity and tail of each matrix, and for tai50b the cluster
structure of the distance matrix and the split of the best known objective
(QAPSOLN/tai50b.sln) into within-cluster and between-cluster pairs.
"""
import numpy as np


def load(fn):
    v = open(fn).read().split()
    n = int(v[0])
    a = np.array(v[1:1 + 2 * n * n], dtype=float)
    return n, a[:n * n].reshape(n, n), a[n * n:].reshape(n, n)


for nm in ['tho40', 'tai50b', 'wil50', 'sko49', 'sko42']:
    n, A, B = load(f'qaplib/{nm}.dat')
    for lab, M in (('first', A), ('second', B)):
        off = M[~np.eye(n, dtype=bool)]
        nz = np.sort(off[off != 0])[::-1]
        cs = np.cumsum(nz) / nz.sum()
        print(f"{nm:7s} {lab:6s} {'sym ' if np.allclose(M, M.T) else 'asym'} zeros {100 * (off == 0).mean():3.0f}%"
              f" distinct {len(np.unique(off)):4d} max {nz.max():6.0f} median(nz) {np.median(nz):5.0f}"
              f"  top 10% of nonzeros carry {100 * cs[len(nz) // 10]:3.0f}% of the weight, top 1%: {100 * cs[max(len(nz) // 100, 1)]:2.0f}%")

n, D, F = load('qaplib/tai50b.dat')          # first matrix: distances (locations), second: flows


def comps(th):
    lab = -np.ones(n, int)
    c = 0
    for i in range(n):
        if lab[i] >= 0:
            continue
        st = [i]
        lab[i] = c
        while st:
            u = st.pop()
            for w in np.where((D[u] <= th) & (lab < 0))[0]:
                lab[w] = c
                st.append(w)
        c += 1
    return lab


iu = np.triu_indices(n, 1)
print('\ntai50b distances: quantiles 1/5/10/25/50/75/90/99 %:', np.percentile(D[iu], [1, 5, 10, 25, 50, 75, 90, 99]).round(0))
for th in [40, 60, 80, 150, 300]:
    lab = comps(th)
    same = lab[:, None] == lab[None, :]
    w, b = D[iu][same[iu]], D[iu][~same[iu]]
    print(f"single linkage at {th:3d}: {lab.max() + 1:2d} clusters, sizes {sorted(np.bincount(lab).tolist(), reverse=True)},"
          f" within mean {w.mean():.0f} max {w.max():.0f}, between mean {b.mean():.0f} min {b.min():.0f}")
J = np.eye(n) - 1.0 / n
ev = np.linalg.eigvalsh(-0.5 * J @ (D ** 2) @ J)[::-1]
print('classical MDS Gram eigenvalues (top 4):', ev[:4].round(0), ' first axis share of the positive part: %.1f%%' % (100 * ev[0] / ev[ev > 0].sum()))

s = open('QAPSOLN/tai50b.sln').read().replace(',', ' ').split()
val = float(s[1])
p = np.array(s[2:2 + n], int) - 1
T = D * F[np.ix_(p, p)]                       # term (i, j): locations i, j hold facilities p[i], p[j]
assert abs(T.sum() - val) < 0.5, (T.sum(), val)
lab = comps(100)
same = lab[:, None] == lab[None, :]
Fp = F[np.ix_(p, p)]
print(f"best known {val:.0f}: within-cluster pairs carry {100 * T[same].sum() / T.sum():.2f}% of the objective,"
      f" {100 * Fp[same].sum() / Fp.sum():.1f}% of the flow, and are {100 * (same.sum() - n) / (n * n - n):.1f}% of the pairs")
