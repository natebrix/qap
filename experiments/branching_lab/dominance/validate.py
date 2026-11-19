"""Brute-force validation of the pair and one-fixed exchange-dominance tests.
For each instance: enumerate a solution space (all n! for n <= 10; the completions of a few random
3-fix nodes for n = 12), draw many nodes, and check for each:
  (1) exactness: the formula's min Delta (fast sort form and general LAP form) equals the minimum of
      cost(p) - cost(swap(p)) over every completion p, for every pair (and one-fixed: every (k,i));
  (2) claim: if the test says dominated (min >= 1), min cost over the node > optimum value, and
      every completion's exchange image is strictly cheaper.
Output: validate.out"""
import itertools, sys, numpy as np, math
from qapdata import read_dat, read_sln
from dom import Inst, pair_min_general
rng = np.random.default_rng(20261119)

def all_costs(A, B, P):
    out = np.zeros(len(P), dtype=np.int64)
    for s in range(0, len(P), 200000):
        Q = P[s:s + 200000].astype(np.int64)
        out[s:s + 200000] = (A[None] * B[Q[:, :, None], Q[:, None, :]]).sum((1, 2))
    return out

def perms_of(rows, cols, F, n):
    """all completions of F as an array of full permutations"""
    free_r = [l for l in range(n) if l not in F]
    free_c = [c for c in range(n) if c not in F.values()]
    base = np.full(n, -1, dtype=np.int8)
    for r, c in F.items(): base[r] = c
    pc = np.array(list(itertools.permutations(free_c)), dtype=np.int8) if free_c else np.zeros((1, 0), np.int8)
    P = np.repeat(base[None], len(pc), 0); P[:, free_r] = pc
    return P

stats = {}
def check_node(name, inst, F, P, costs, opt, S):
    A, B = inst.A, inst.B
    n = inst.n
    sel = np.all(P[:, list(F.keys())] == np.array(list(F.values()), dtype=np.int8), 1) if F else np.ones(len(P), bool)
    PN, cN = P[sel], costs[sel]
    mins, pairs = inst.pair_mins(F)
    for (i, k), mv in zip(pairs, mins):
        Q = PN.copy(); Q[:, [i, k]] = Q[:, [k, i]]
        cQ = all_costs(A, B, Q)
        true = int((cN - cQ).min())
        S['pairs'] += 1
        if true != int(mv): S['pair_inexact'] += 1; print('INEXACT pair', name, F, i, k, mv, true)
        g = pair_min_general(A, B, F, int(i), int(k))
        if g != true: S['general_inexact'] += 1; print('INEXACT general', name, F, i, k, g, true)
        if mv >= 1:
            S['pair_dominated'] += 1
            if cN.min() <= opt: S['pair_counter'] += 1; print('COUNTEREXAMPLE', name, F, i, k)
            if not (cQ < cN).all(): S['image_fail'] += 1
    if len(F) >= 1 and len(F) < n:
        best, arg = inst.onefixed(F)
        # brute force: for every (k fixed, i free), min over completions of cost(p) - cost(swap_{i,k}(p))
        tb = None
        for k in F:
            for i in range(n):
                if i in F: continue
                Q = PN.copy(); Q[:, [i, k]] = Q[:, [k, i]]
                v = int((cN - all_costs(A, B, Q)).min())
                if tb is None or v > tb: tb = v
        S['onefixed'] += 1
        if int(round(best)) != tb: S['onefixed_inexact'] += 1; print('INEXACT onefixed', name, F, best, tb)
        if best >= 1:
            S['onefixed_dominated'] += 1
            if cN.min() <= opt: S['onefixed_counter'] += 1; print('COUNTEREXAMPLE onefixed', name, F, arg)
    S['nodes'] += 1
    S['nodes_with_opt'] += int(cN.min() == opt)

def random_nodes(n, base_F, count, P, costs, opt):
    """random nodes extending base_F; half of them biased: fixes taken from a near-optimal permutation
    with one exchange (the mirror situation)."""
    out = []
    good = P[np.argsort(costs)[:50]]
    for t in range(count):
        free = [l for l in range(n) if l not in base_F]
        d = int(rng.integers(1, len(free)))
        rows = list(rng.choice(free, d, replace=False))
        if t % 2 == 0:
            p = good[rng.integers(len(good))].copy()
            if t % 4 == 0:
                u, v = rng.choice(free, 2, replace=False); p[[u, v]] = p[[v, u]]
                rows = list(set(rows) | {int(u), int(v)})
            F = dict(base_F); F.update({int(r): int(p[r]) for r in rows})
        else:
            freec = [c for c in range(n) if c not in base_F.values()]
            cs = rng.permutation(freec)[:len(rows)]
            F = dict(base_F); F.update({int(r): int(c) for r, c in zip(rows, cs)})
        out.append(F)
    return out

def run(name, A, B, opt=None, full=True, nnodes=300):
    inst = Inst(A, B); n = inst.n
    S = dict.fromkeys(['nodes', 'nodes_with_opt', 'pairs', 'pair_inexact', 'general_inexact', 'pair_dominated', 'pair_counter',
                       'image_fail', 'onefixed', 'onefixed_inexact', 'onefixed_dominated', 'onefixed_counter'], 0)
    if full:
        P = perms_of(None, None, {}, n); costs = all_costs(inst.A, inst.B, P)
        o = int(costs.min()); nopt = int((costs == o).sum())
        if opt is not None: assert o == opt, (name, o, opt)
        opt = o
        for F in random_nodes(n, {}, nnodes, P, costs, opt):
            check_node(name, inst, F, P, costs, opt, S)
    else:
        nopt = '?'
        for rep in range(3):
            rows = rng.choice(n, 3, replace=False)
            _, sp = read_sln(name)
            # one base node from the optimum, two random
            cs = sp[rows] if rep == 0 else rng.choice(n, 3, replace=False)
            base = {int(r): int(c) for r, c in zip(rows, cs)}
            P = perms_of(None, None, base, n); costs = all_costs(inst.A, inst.B, P)
            for F in random_nodes(n, base, nnodes // 3, P, costs, opt):
                check_node(name, inst, F, P, costs, opt, S)
    print(f"{name:10s} n={n} Asym={inst.Asym} Bsym={inst.Bsym} diag0={(np.diag(inst.A)==0).all() and (np.diag(inst.B)==0).all()} "
          f"opt={opt} #opt={nopt} " + ' '.join(f'{k}={v}' for k, v in S.items()), flush=True)
    return S

if __name__ == '__main__':
    tot = {}
    def acc(S):
        for k, v in S.items(): tot[k] = tot.get(k, 0) + v
    for nm in ['nug05', 'nug06', 'nug07', 'nug08', 'esc08a', 'esc08b', 'esc08c', 'esc08d', 'esc08e', 'esc08f', 'law07']:
        acc(run(nm, *read_dat(nm), nnodes=300))
    for nm in ['tai10a', 'tai10b', 'had10', 'scr10', 'rou10', 'lipa10a', 'lipa10b', 'car10ga', 'car10pa']:
        acc(run(nm, *read_dat(nm), nnodes=150))
    for nm in ['tai12b', 'tai12a', 'nug12', 'had12', 'scr12', 'chr12a', 'chr12b']:
        A, B = read_dat(nm); opt, _ = read_sln(nm)
        acc(run(nm, A, B, opt=opt, full=False, nnodes=150))
    # random instances with ties: small integers, asymmetric, nonzero diagonals; and symmetric ones
    for t in range(12):
        n = 7 if t < 6 else 8
        A = rng.integers(0, 3, (n, n)); B = rng.integers(0, 3, (n, n))
        if t % 3 == 1: A = A + A.T
        if t % 3 == 2: B = B + B.T; A = np.triu(A) + np.triu(A, 1).T; np.fill_diagonal(A, 0)
        acc(run(f'rand{t}', A, B, nnodes=200))
    print('TOTAL', ' '.join(f'{k}={v}' for k, v in tot.items()))
