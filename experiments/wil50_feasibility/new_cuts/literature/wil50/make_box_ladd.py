#!/usr/bin/env python3
"""L = mu * sum_k G_k over the top box cuts of box_gate.py's json (each
<G_k, Y> <= h_k with G_k: 2b on (aj),(aj) for a in P, j in Q and -1 on
(aj),(ck) for a != c in P, j != k in Q; h_k = b(b+1)), and the constant
mu * sum h_k for -Ladd L.bin:CONST. Checks every cut on random permutations
of the node (Y = x x') before writing: a cut violated by a permutation aborts.
  python3 -I make_box_ladd.py --json J --mu MU --top K --out L.bin [--minv V]"""
import argparse, json
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--json', required=True); ap.add_argument('--mu', type=float, required=True)
ap.add_argument('--top', type=int, default=50); ap.add_argument('--out', required=True)
ap.add_argument('--minv', type=float, default=0.0); ap.add_argument('--nperm', type=int, default=3000)
a = ap.parse_args()
d = json.load(open(a.json)); m = d['m']; N = m * m + 1
cuts = [c for c in d['top'] if c[0] > a.minv][:a.top]
L = np.zeros((N, N)); const = 0.0; tight = 0
rng = np.random.default_rng(7)
perms = [rng.permutation(m) for _ in range(a.nperm)]   # perm[a] = facility at local location a
for (v, nm, Q, b, P) in cuts:
    idx = np.array([1 + p * m + j for p in P for j in Q])
    G = np.zeros((N, N))
    G[np.ix_(idx, idx)] = -1.0
    for p in P:            # same location, different facilities: gangster zero anyway; keep the formula exact
        for j in Q:
            for k in Q:
                if j != k: G[1 + p * m + j, 1 + p * m + k] = 0.0
    for j in Q:            # different locations, same facility: gangster zero
        for p in P:
            for c in P:
                if p != c: G[1 + p * m + j, 1 + c * m + j] = 0.0
    for e in idx:
        G[e, e] = 2.0 * b
    h = float(b * (b + 1))
    # validity on permutations Y = x x': random ones plus, for every count
    # t = 0..min(|P|,|Q|), a permutation placing exactly t facilities of Q in P
    worst = -np.inf
    Pl = list(P); Ql = list(Q); notP = [x for x in range(m) if x not in P]; notQ = [j for j in range(m) if j not in Q]
    tests = list(perms)
    for t in range(0, min(len(Pl), len(Ql)) + 1):
        if len(Ql) - t > len(notP) or len(Pl) - t > len(notQ):
            continue
        pi = -np.ones(m, dtype=int)
        locs_in = rng.permutation(Pl); locs_out = rng.permutation(notP)
        facs = rng.permutation(Ql)
        for i in range(t): pi[locs_in[i]] = facs[i]
        for i in range(len(Ql) - t): pi[locs_out[i]] = facs[t + i]
        rest_f = list(rng.permutation(notQ)); free_l = [x for x in range(m) if pi[x] < 0]
        for x, f in zip(free_l, rest_f): pi[x] = f
        tests.append(pi)
    for pi in tests:
        idxp = 1 + np.arange(m) * m + pi
        val = G[np.ix_(idxp, idxp)].sum()
        worst = max(worst, val - h)
    assert worst <= 1e-9, ("INVALID cut", nm, Q, b, worst)
    tight += abs(worst) < 1e-9
    L += a.mu * G; const += a.mu * h
L.tofile(a.out)
print(f"{len(cuts)} cuts ({tight} tight on a test permutation; all checked on {a.nperm} random permutations), violations {cuts[-1][0]:.4f}..{cuts[0][0]:.4f}, "
      f"mu {a.mu}, nonzeros {np.count_nonzero(L)}, |L| {np.abs(L).sum():.3f}, const {const:.6f} -> use: -Ladd {a.out}:{const:.9g}")
