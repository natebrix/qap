#!/usr/bin/env python3
"""Cross-check every node the solver drops (trace outcome X) with dom.pair = 1.

  python3 xcheck_drops.py OUT.txt UBMODE inst [inst ...]     (UBMODE opt1 | poor)

For each instance: one bnbdrv run with QAP_NODE_TRACE, then for every X node:
  ref    max over pairs of fixed rows of the min Delta over completions, by
         this file's own general formula (no symmetry assumed, the linear term
         C included, on the original data), then by dominance/dom.py's
         Inst.pair_test where the instance has no C (dom.py has no C); both > 0
  brute  where the node has at most BF_MAX free rows: over every completion,
         (a) the same max-min by enumeration equals ref, and (b) min cost over
         the completions > the optimum (the claim: no optimal solution inside)
Writes one line per instance (and every failure) to OUT.txt.
"""
import itertools, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dominance"))
from dom import Inst
WT = "/home/user/wt-dombuild"
SCR = "/tmp/claude-0/-home-user/92b66c47-7a9f-5ca5-8a65-57e14ac5b075/scratchpad"
BF_MAX = int(os.environ.get("BF_MAX", "8"))
sys.path.insert(0, HERE)
from run_valid import ref_opt


def read_solver_data(name):
    """A, B, C as bnbdrv's readQAPLIB reads them (C from trailing tokens if n*n of them)"""
    t = open(f"{WT}/qaplib/{name}.dat").read().split()
    n = int(t[0]); v = [int(float(x)) for x in t[1:]]
    A = np.array(v[:n * n], dtype=np.int64).reshape(n, n)
    B = np.array(v[n * n:2 * n * n], dtype=np.int64).reshape(n, n)
    C = np.zeros((n, n), dtype=np.int64)
    if len(v) == 3 * n * n:
        C = np.array(v[2 * n * n:], dtype=np.int64).reshape(n, n)
    elif len(v) != 2 * n * n:
        raise SystemExit(f"{name}: trailing tokens, skipped")
    return A, B, C


def cost(A, B, C, p):
    p = np.asarray(p)
    return int((A * B[np.ix_(p, p)]).sum() + C[np.arange(len(p)), p].sum())


def vcost(A, B, C, P):
    """cost of every row of P (a permutation per row)"""
    out = C[np.arange(P.shape[1])[None, :], P].sum(1)
    for u in range(P.shape[1]):
        out = out + (A[u][None, :] * B[P[:, u][:, None], P]).sum(1)
    return out


def pair_min_general(A, B, C, F, i, k):
    from scipy.optimize import linear_sum_assignment
    n = len(A); a, b = F[i], F[k]
    rows = [l for l in range(n) if l not in F]
    used = set(F.values()); cols = [c for c in range(n) if c not in used]
    v = int((A[i, i] - A[k, k]) * (B[a, a] - B[b, b]) + (A[i, k] - A[k, i]) * (B[a, b] - B[b, a])
            + C[i, a] + C[k, b] - C[i, b] - C[k, a])
    for l, c in F.items():
        if l in (i, k): continue
        v += int((A[i, l] - A[k, l]) * (B[a, c] - B[b, c]) + (A[l, i] - A[l, k]) * (B[c, a] - B[c, b]))
    if not rows: return v
    R = np.array(rows); Cc = np.array(cols)
    M = np.outer(A[i, R] - A[k, R], B[a, Cc] - B[b, Cc]) + np.outer(A[R, i] - A[R, k], B[Cc, a] - B[Cc, b])
    r, c = linear_sum_assignment(M)
    return v + int(M[r, c].sum())


def main():
    out, mode = sys.argv[1], sys.argv[2]
    for name in sys.argv[3:]:
        A, B, C = read_solver_data(name)
        n = len(A)
        opt = ref_opt(name, n)[0]
        ub = opt + 1 if mode == "opt1" else 2 * abs(opt) + 100
        tr = f"{SCR}/{name}.xtrace"
        if os.path.exists(tr): os.remove(tr)
        env = dict(os.environ, QAP_NODE_TRACE=tr, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        subprocess.run([f"{WT}/latest/bnbdrv", "-p", f"{HERE}/par/default_dom.par", "-u", str(ub),
                        f"{WT}/qaplib/{name}.dat"], cwd=f"{WT}/latest", env=env, capture_output=True)
        xs = []
        for ln in open(tr):
            t = ln.split()
            if t[7] == "X":
                xs.append(dict(tuple(map(int, s.split(":"))) for s in t[2].split(",")))
        inst = Inst(A, B)
        nodom = nref = nbf = fails = 0
        lines = []
        for F in xs:
            keys = list(F)
            ref = max(pair_min_general(A, B, C, F, keys[u], keys[w])
                      for u in range(len(keys)) for w in range(u + 1, len(keys)))
            ok = ref > 0
            if not C.any():
                pm, _ = inst.pair_test(F)
                ok = ok and pm == ref
            nref += 1
            free_r = [l for l in range(n) if l not in F]
            free_c = [c for c in range(n) if c not in F.values()]
            if ok and len(free_r) <= BF_MAX:
                nbf += 1
                P = np.empty((0, n), dtype=np.int64)
                perms = np.array(list(itertools.permutations(free_c)), dtype=np.int64).reshape(-1, len(free_r))
                P = np.empty((len(perms), n), dtype=np.int64)
                for l, c in F.items(): P[:, l] = c
                P[:, free_r] = perms
                costs = vcost(A, B, C, P)
                bf = None
                for u in range(len(keys)):
                    for w in range(u + 1, len(keys)):
                        i, k = keys[u], keys[w]
                        Q = P.copy(); Q[:, [i, k]] = P[:, [k, i]]
                        m = int((costs - vcost(A, B, C, Q)).min())
                        bf = m if bf is None else max(bf, m)
                ok = ok and bf == ref and costs.min() > opt
                if not ok: lines.append(f"FAIL brute {name} {F} ref={ref} bf={bf} mincost={costs.min()} opt={opt}")
            elif not ok:
                lines.append(f"FAIL ref {name} {F} ref={ref}")
            fails += not ok
        summ = (f"{name} ({mode}, ub {ub}, opt {opt}): {len(xs)} dropped nodes; all > 0 by the general formula"
                f"{' and equal to dom.py' if not C.any() else ' (C present: dom.py not applicable)'}: {nref - fails if not lines else '?'}/{nref}; "
                f"brute-forced {nbf} (<= {BF_MAX} free): max-min equal and min cost > opt; failures {fails}")
        print(summ, flush=True)
        with open(out, "a") as f:
            f.write(summ + "\n")
            for l in lines: f.write(l + "\n")
        os.remove(tr)


if __name__ == "__main__":
    main()
