#!/usr/bin/env python3
"""Stage 2b: the ceiling of the cuts AS LINEAR INEQUALITIES ON Y, exactly.

The exact ceiling of stage 2 (v_ex1: the best exchange-optimal completion) is what the logical
condition gives. A relaxation only sees the linear inequalities x_jb * delta >= 0 on the lifted
matrix, and even a perfect relaxation (the convex hull of the lifted completions xx^T) can
satisfy them with a MIXTURE of completions none of which is exchange-optimal. That ceiling is
an LP over the completions:
    v_lp = min sum_p w_p C(p)   s.t.  sum_p w_p [p(j) = b] delta_{i,j}(p) >= 0  for every cut
                                      (fixed i, free j, free b),  sum w = 1, w >= 0,
with v_node <= v_lp <= v_ex1. No bound from these cuts in any relaxation can pass v_lp.
  python3 lpceil.py > logs/lpceil.out        (m <= 9 free: at most 362,880 columns)
"""
import itertools, os, sys, random, statistics as st
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, coo_matrix
from swaplib import read_dat

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
here = os.path.dirname(os.path.abspath(__file__))


def costs(A, B, P):
    n = A.shape[0]
    c = np.zeros(P.shape[0])
    for k in range(n):
        for l in range(n):
            if A[k, l] != 0:
                c += A[k, l] * B[P[:, k], P[:, l]]
    return c


def node_lp(A, B, fixed):
    n = A.shape[0]
    ff = [k for k in range(n) if k not in fixed]
    fl = [c for c in range(n) if c not in fixed.values()]
    m = len(ff)
    perms = np.array(list(itertools.permutations(fl)), dtype=np.int64)
    P = np.zeros((perms.shape[0], n), dtype=np.int64)
    for k, a in fixed.items():
        P[:, k] = a
    for t, k in enumerate(ff):
        P[:, k] = perms[:, t]
    c = costs(A, B, P)
    N = P.shape[0]
    locidx = {c_: t for t, c_ in enumerate(fl)}
    lmap = np.zeros(n, dtype=np.int64)
    for c_, t in locidx.items():
        lmap[c_] = t
    rows, cols, vals = [], [], []
    ok = np.ones(N, dtype=bool)
    r0 = 0
    for i, a in fixed.items():
        for jt, j in enumerate(ff):
            P2 = P.copy()
            P2[:, i] = P[:, j]
            P2[:, j] = a
            d = costs(A, B, P2) - c
            ok &= d >= 0
            rows.append(r0 + jt * m + lmap[P[:, j]])
            cols.append(np.arange(N))
            vals.append(d)
        r0 += m * m
    G = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r0, N)).tocsr()
    vnode = c.min()
    vex = c[ok].min() if ok.any() else float("inf")
    # only columns below the exchange optimum can matter when it is finite; keep all otherwise
    res = linprog(c, A_ub=-G, b_ub=np.zeros(r0), A_eq=np.ones((1, N)), b_eq=[1.0], bounds=(0, None), method="highs")
    vlp = res.fun if res.status == 0 else float("nan")
    return vnode, vlp, vex, int(ok.sum()), N


def main():
    random.seed(11)
    print("# inst depth m spec v_node v_lp v_ex1 exchange-optimal/completions onpath")
    summary = {}
    jobs = []
    # n = 10: depths 1 and 2 (m = 9, 8): the shallow regime (one or two fixed, most free)
    for inst in ("tai10a", "rou10", "scr10", "had10", "tai10b"):
        path = os.path.join(ROOT, "qaplib", inst + ".dat")
        if not os.path.exists(path):
            continue
        n, A, B = read_dat(path)
        # the optimum by enumeration through the depth-1 nodes of facility 0
        best = (float("inf"), None)
        fl = list(range(n))
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
        c = costs(A, B, perms)
        popt = perms[c.argmin()]
        opt = c.min()
        del perms, c
        order = list(range(n)); random.shuffle(order)
        for line in range(2):
            j = order[line]
            for b in range(n):
                jobs.append((inst, A, B, opt, 1, {j: b}, b == popt[j]))
            i = order[line + 2]
            j2 = order[line + 4]
            for b in range(n):
                if b != popt[i]:
                    jobs.append((inst, A, B, opt, 2, {i: int(popt[i]), j2: b}, b == popt[j2]))
    # n = 12: depths 3 to 5 from the stage 2 nodes (m = 9, 8, 7)
    for inst in ("nug12", "scr12", "tai12a", "rou12", "had12", "tai12b"):
        n, A, B = read_dat(os.path.join(ROOT, "qaplib", inst + ".dat"))
        rows = []; opt = None
        for line in open(os.path.join(here, "logs", "stage2_%s.out" % inst)):
            if line.startswith("#"):
                if "optimum" in line:
                    opt = float(line.split("optimum")[1].split()[0])
                continue
            t = line.split()
            if len(t) != 11 or int(t[1]) not in (3, 4, 5):
                continue
            fixed = {}
            for tok in ([] if t[3] == "-" else t[3].split(",")) + ["%s:%s" % (t[4], t[5])]:
                k, a = tok.split(":"); fixed[int(k)] = int(a)
            rows.append((int(t[1]), fixed, int(t[10])))
        random.shuffle(rows)
        cnt = {3: 0, 4: 0, 5: 0}
        for d, fixed, onp in rows:
            if cnt[d] >= (6 if d == 3 else 12):
                continue
            cnt[d] += 1
            jobs.append((inst, A, B, opt, d, fixed, onp))
    for inst, A, B, opt, d, fixed, onp in jobs:
        vnode, vlp, vex, nok, N = node_lp(A, B, fixed)
        spec = ",".join("%d:%d" % kv for kv in fixed.items())
        print("%s %d %d %s %.0f %.3f %s %d/%d %d" % (inst, d, A.shape[0] - d, spec, vnode, vlp, "inf" if vex == float("inf") else "%.0f" % vex, nok, N, onp))
        sys.stdout.flush()
        if not onp:
            n = A.shape[0]
            key = ("n=%d" % n, d)
            e = summary.setdefault(key, [])
            e.append((vnode, vlp, vex, opt))
    print("\n# off-path nodes: share with v_ex1 > v_node (logical ceiling), share with v_lp > v_node (the linear cuts' ceiling),")
    print("# median of (v_lp - v_node) / opt over all nodes and over the affected ones, median share of the logical rise the LP keeps (finite v_ex1)")
    print("# %-6s %5s %5s %9s %9s %12s %12s %10s" % ("n", "depth", "nodes", "ex>node", "lp>node", "med rise%", "med rise%aff", "lp/ex"))
    for (nn, d), e in sorted(summary.items()):
        k = len(e)
        aex = sum(v[2] > v[0] + 1e-6 for v in e)
        alp = [v for v in e if v[1] > v[0] + 1e-6 * abs(v[0]) + 1e-6]
        rise = [100 * (v[1] - v[0]) / v[3] for v in e]
        risea = [100 * (v[1] - v[0]) / v[3] for v in alp]
        share = [(v[1] - v[0]) / (v[2] - v[0]) for v in e if v[0] + 1e-6 < v[2] < float("inf")]
        print("# %-6s %5d %5d %8.0f%% %8.0f%% %12.3f %12s %10s" % (nn, d, k, 100 * aex / k, 100 * len(alp) / k, st.median(rise),
              "%.3f" % st.median(risea) if risea else "-", "%.2f" % st.median(share) if share else "-"))


if __name__ == "__main__":
    main()
