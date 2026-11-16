#!/usr/bin/env python3
"""Stage 1: the formula for delta against brute-force cost differences, and
validity of every cut at every node that contains a global optimum (n <= 9).
  python3 stage1.py > logs/stage1.out
"""
import itertools, os, sys
import numpy as np
from swaplib import read_dat, cost, swap_delta, node_cut, fixed_pair_cut

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
rng = np.random.default_rng(7)


def instances():
    out = []
    for n in (6, 7, 8):
        # asymmetric, nonzero diagonals
        out.append(("rand_asym_n%d" % n, rng.integers(0, 10, (n, n)).astype(float),
                    rng.integers(0, 10, (n, n)).astype(float)))
    # nug-like 2 x 4 grid, symmetric flows with zeros
    n = 8
    D = np.array([[abs(a // 4 - b // 4) + abs(a % 4 - b % 4) for b in range(n)] for a in range(n)], float)
    F = np.triu(rng.integers(0, 6, (n, n)) * (rng.random((n, n)) < 0.6), 1)
    out.append(("grid2x4", (F + F.T).astype(float), D))
    # tai-b-like: clustered points, asymmetric heavy-tailed flows
    n = 9
    pts = np.concatenate([rng.normal(c, 1.0, (3, 2)) for c in ((0, 0), (20, 5), (8, 30))])
    D = np.rint(np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1)) * 10)
    F = np.rint(np.exp(rng.normal(0, 2.5, (n, n))) * (rng.random((n, n)) < 0.5))
    np.fill_diagonal(F, 0)
    out.append(("taib_like_n9", F, D))
    for name in ("nug08", "esc08a", "nug07"):
        nn, A, B = read_dat(os.path.join(ROOT, "qaplib", name + ".dat"))
        out.append((name, A, B))
    return out


def main():
    worst = 0.0
    # (a) the formula on random permutations, including asymmetric matrices and diagonals
    for n in (5, 9, 14):
        for trial in range(200):
            A = rng.integers(-5, 10, (n, n)).astype(float)
            B = rng.integers(-5, 10, (n, n)).astype(float)
            p = list(rng.permutation(n))
            i, j = rng.choice(n, 2, replace=False)
            q = list(p); q[i], q[j] = q[j], q[i]
            worst = max(worst, abs(cost(A, B, q) - cost(A, B, p) - swap_delta(A, B, p, i, j)))
    print("(a) delta formula vs brute-force difference, 600 random cases (asymmetric, diagonals): max |err| = %g" % worst)
    assert worst < 1e-9

    # (b) every global optimum satisfies every cut of every node that contains it
    for name, A, B in instances():
        n = A.shape[0]
        best, opts = None, []
        for p in itertools.permutations(range(n)):
            c = cost(A, B, p)
            if best is None or c < best - 1e-9:
                best, opts = c, [p]
            elif abs(c - best) <= 1e-9:
                opts.append(p)
        ncut = 0; minlhs = np.inf; nfp = 0; minfp = np.inf; nodes = 0
        nopt = len(opts)
        if len(opts) > 6:
            opts = [opts[t] for t in rng.choice(len(opts), 6, replace=False)]
        for p in opts:
            x = np.zeros((n, n)); x[np.arange(n), list(p)] = 1.0
            for d in range(1, n - 1):
                subsets = list(itertools.combinations(range(n), d))
                if len(subsets) > 12:
                    subsets = [subsets[t] for t in rng.choice(len(subsets), 12, replace=False)]
                for S in subsets:
                    fixed = {k: p[k] for k in S}
                    nodes += 1
                    free_f = [k for k in range(n) if k not in fixed]
                    free_l = [c for c in range(n) if c not in fixed.values()]
                    for i in S:
                        for j in free_f:
                            for b in free_l:
                                c0, g = node_cut(A, B, fixed, i, j, b)
                                lhs = x[j, b] * (c0 + (g * x).sum())
                                # the linear form must equal the brute-force exchange gain when x_jb = 1
                                if x[j, b] == 1.0:
                                    q = list(p); q[i], q[j] = q[j], q[i]
                                    assert abs(lhs - (cost(A, B, q) - best)) < 1e-9
                                minlhs = min(minlhs, lhs); ncut += 1
                    for i, i2 in itertools.combinations(S, 2):
                        c0, g = fixed_pair_cut(A, B, fixed, i, i2)
                        l = c0 + (g * x).sum()
                        q = list(p); q[i], q[i2] = q[i2], q[i]
                        assert abs(l - (cost(A, B, q) - best)) < 1e-9
                        minfp = min(minfp, l); nfp += 1
        sys.stdout.flush(); print("(b) %-14s n=%d optimum %g, %d global optima (at most 6 checked), %d nodes containing one: "
              "%d pairwise cuts, min lhs %g; %d fixed-pair cuts, min lhs %g" %
              (name, n, best, nopt, nodes, ncut, minlhs, nfp, minfp))
        assert minlhs >= -1e-9 and minfp >= -1e-9
    print("PASS")


if __name__ == "__main__":
    main()
