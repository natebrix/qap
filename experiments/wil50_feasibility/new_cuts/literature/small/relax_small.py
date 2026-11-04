#!/usr/bin/env python3
"""Exact comparison of QAP relaxations on small instances (n <= 12).

Answers: on top of the DNN, how much do (a) the regional triangles, (b) the
3-local marginal polytope on side-1 triples (Astra's "triple consistency",
facility triples for wil50), (c) the same on side-0 triples, (d) both, and
(e) the full projected RLT2 (one cubic variable per assignment triple shared
by both sides) add? And how strong are RLT1 / RLT2 without the PSD cone?

Convention: x_ij = 1 iff side-0 index i (rows of the first QAPLIB matrix A)
goes to side-1 index j (rows of B); objective sum A_ik B_jl x_ij x_kl (+ C).
Lifted index of x_ij is 1 + i*n + j, index 0 the constant. Y = Vh R Vh',
R PSD (facial reduction, ZKRW), Y >= 0, gangster zeros, arrow.

Cubic variables z_{abc}, one per unordered triple of assignments with
distinct side-0 and distinct side-1 indices. For each compatible pair
{a, b}:
  side-0 marginals: for each third side-0 index p:  sum_q z_{ab(pq)} = Y_ab
  side-1 marginals: for each third side-1 index q:  sum_p z_{ab(pq)} = Y_ab
"LMP3 side 1" keeps only the side-1 marginal equations... careful: the
distribution of the side-0 images of a side-1 triple {j_a, j_b, q} is the set
of z with that side-1 triple, and its marginal constraints are the equations
"for third side-1 index q" (sum over p). So:
  lmp1 = side-1 triples consistent  = the 'sum_p' equations,
  lmp0 = side-0 triples consistent  = the 'sum_q' equations,
  rlt2 = both, on the same z.
For wil50 (A = grid distance, B = flow) side 0 = locations, side 1 =
facilities, and Astra's triple cuts are lmp1.

Usage: relax_small.py inst.dat --rows R --cols C [--which dnn,tri,lmp1,...]
"""
import argparse, itertools, json, sys, time
import numpy as np
import scipy.sparse as sp
import cvxpy as cp


def read_qaplib(path):
    tok = open(path).read().split()
    n = int(tok[0]); v = list(map(float, tok[1:]))
    A = np.array(v[:n * n]).reshape(n, n)
    B = np.array(v[n * n:2 * n * n]).reshape(n, n)
    C = None
    if len(v) >= 3 * n * n:
        C = np.array(v[2 * n * n:3 * n * n]).reshape(n, n)
    return n, A, B, C


def write_qaplib(path, A, B):
    n = len(A)
    with open(path, "w") as f:
        f.write(f"{n}\n\n")
        for M in (A, B):
            for r in M:
                f.write(" ".join(str(int(round(x))) for x in r) + "\n")
            f.write("\n")


def lifted_cost(n, A, B, C=None):
    N = n * n + 1
    L = np.zeros((N, N))
    # L[(ij),(kl)] = A_ik B_jl
    L[1:, 1:] = np.einsum("ik,jl->ijkl", A, B).reshape(n * n, n * n)
    if C is not None:
        for i in range(n):
            for j in range(n):
                L[0, 1 + i * n + j] += C[i, j] / 2
                L[1 + i * n + j, 0] += C[i, j] / 2
    return L


def vhat(n):
    """Vh = [[1, 0], [e/n (x) e/n ... ]]: x = e e'/n + V X V'; sparse V."""
    V = np.zeros((n, n - 1))
    for k in range(n - 1):
        V[k, k] = 1.0
        V[n - 1, k] = -1.0
    VV = np.kron(V, V)
    Vh = np.zeros((n * n + 1, (n - 1) ** 2 + 1))
    Vh[0, 0] = 1.0
    Vh[1:, 0] = 1.0 / n
    Vh[1:, 1:] = VV
    return sp.csr_matrix(Vh)


def gangster_mask(n):
    N = n * n + 1
    G = np.zeros((N, N), dtype=bool)
    for i in range(n):
        for j in range(n):
            a = 1 + i * n + j
            for l in range(n):
                if l != j:
                    G[a, 1 + i * n + l] = True
            for k in range(n):
                if k != i:
                    G[a, 1 + k * n + j] = True
    return G


def triples(n):
    """assignment triples {(i1,j1),(i2,j2),(i3,j3)}, i1<i2<i3, j's distinct"""
    T = []
    for i3 in itertools.combinations(range(n), 3):
        for j3 in itertools.permutations(range(n), 3):
            T.append((i3, j3))
    return T


def build_cubic(n, sides):
    """sparse rows: sum over z (coef 1) - Y_ab (coef -1, as y index a*N+b)"""
    N = n * n + 1
    T = triples(n)
    zid = {}
    for t, (ii, jj) in enumerate(T):
        zid[frozenset(zip(ii, jj))] = t
    rows_z, cols_z, rows_y, cols_y = [], [], [], []
    r = 0
    for (i1, i2) in itertools.combinations(range(n), 2):
        for j1 in range(n):
            for j2 in range(n):
                if j1 == j2:
                    continue
                a = 1 + i1 * n + j1; b = 1 + i2 * n + j2
                if 0 in sides:   # third side-0 index p: sum over q
                    for p in range(n):
                        if p in (i1, i2):
                            continue
                        for q in range(n):
                            if q in (j1, j2):
                                continue
                            rows_z.append(r); cols_z.append(zid[frozenset([(i1, j1), (i2, j2), (p, q)])])
                        rows_y.append(r); cols_y.append(a * N + b)
                        r += 1
                if 1 in sides:   # third side-1 index q: sum over p
                    for q in range(n):
                        if q in (j1, j2):
                            continue
                        for p in range(n):
                            if p in (i1, i2):
                                continue
                            rows_z.append(r); cols_z.append(zid[frozenset([(i1, j1), (i2, j2), (p, q)])])
                        rows_y.append(r); cols_y.append(a * N + b)
                        r += 1
    Az = sp.csr_matrix((np.ones(len(rows_z)), (rows_z, cols_z)), shape=(r, len(T)))
    Ay = sp.csr_matrix((-np.ones(len(rows_y)), (rows_y, cols_y)), shape=(r, N * N))
    return Az, Ay, len(T)


def grid_regions(n, rows, cols):
    """threshold regions of a rows x cols grid over side-0 indices i = r*cols + c"""
    assert rows * cols == n
    R = []
    for t in range(1, rows):
        R.append(("row<%d" % t, [r * cols + c for r in range(t) for c in range(cols)]))
    for t in range(1, cols):
        R.append(("col<%d" % t, [r * cols + c for r in range(rows) for c in range(t)]))
    return R


def regional_triangles(n, regions):
    """rows of G (N*N) and h: <G, Y> <= h for every permutation"""
    N = n * n + 1
    data, ri, ci, h = [], [], [], []
    r = 0
    def add(coefs, rhs):
        nonlocal r
        for (idx, v) in coefs.items():
            ri.append(r); ci.append(idx); data.append(v)
        h.append(rhs); r += 1
    for name, S in regions:
        def P(i):  # p_i: sum_{a in S} Y_(ai),(ai)
            return {(1 + a * n + i) * N + (1 + a * n + i): 1.0 for a in S}
        def W(i, j):
            d = {}
            for a in S:
                for b in S:
                    if a == b:
                        continue   # gangster for i != j
                    k = (1 + a * n + i) * N + (1 + b * n + j)
                    d[k] = d.get(k, 0) + 1.0
            return d
        for (i, j, k) in itertools.combinations(range(n), 3):
            p = {x: P(x) for x in (i, j, k)}
            w = {(i, j): W(i, j), (i, k): W(i, k), (j, k): W(j, k)}
            def lin(*terms):
                d = {}
                for s, t in terms:
                    for kk, v in t.items():
                        d[kk] = d.get(kk, 0) + s * v
                return d
            # apex i: W_ij + W_ik - W_jk - p_i <= 0, and the other apexes
            add(lin((1, w[(i, j)]), (1, w[(i, k)]), (-1, w[(j, k)]), (-1, p[i])), 0.0)
            add(lin((1, w[(i, j)]), (1, w[(j, k)]), (-1, w[(i, k)]), (-1, p[j])), 0.0)
            add(lin((1, w[(i, k)]), (1, w[(j, k)]), (-1, w[(i, j)]), (-1, p[k])), 0.0)
            add(lin((1, p[i]), (1, p[j]), (1, p[k]), (-1, w[(i, j)]), (-1, w[(i, k)]), (-1, w[(j, k)])), 1.0)
    G = sp.csr_matrix((data, (ri, ci)), shape=(r, N * N))
    return G, np.array(h)


def solve(n, L, which, regions=None, solver="CLARABEL", verbose=False, psd=True):
    N = n * n + 1
    d = (n - 1) ** 2 + 1
    Vh = vhat(n)
    K = sp.kron(Vh, Vh).tocsr()       # vec(Vh R Vh') = K vec(R), column-major
    R = cp.Variable((d, d), PSD=True) if psd else cp.Variable((d, d), symmetric=True)
    y = cp.Variable(N * N)
    G = gangster_mask(n).reshape(-1)
    cons = [y == K @ cp.vec(R, order="F"), y[0] == 1]
    gz = np.where(G)[0]
    cons.append(y[gz] == 0)
    cons.append(y >= 0)
    # arrow: Y_aa = Y_0a
    aa = np.array([a * N + a for a in range(1, N)]); a0 = np.array([a for a in range(1, N)])
    cons.append(y[aa] == y[a0])
    nz = 0
    if "tri" in which:
        Gt, h = regional_triangles(n, regions)
        cons.append(Gt @ y <= h)
    # each tag gets its own cubic variables: lmp0 and lmp1 together are the two
    # local families decoupled; rlt2 imposes both on one shared z
    for tag, sides in (("lmp0", {0}), ("lmp1", {1}), ("rlt2", {0, 1})):
        if tag in which:
            Az, Ay, nzt = build_cubic(n, sides)
            z = cp.Variable(nzt, nonneg=True)
            cons.append(Az @ z + Ay @ y == 0)
            nz += nzt
    obj = cp.Minimize(L.reshape(-1) @ y)
    prob = cp.Problem(obj, cons)
    t0 = time.time()
    kw = {}
    if solver == "SCS":
        kw = dict(eps_abs=1e-7, eps_rel=1e-7, max_iters=200000)
    prob.solve(solver=solver, verbose=verbose, **kw)
    return prob.value, prob.status, time.time() - t0, y.value


def brute(n, A, B, C=None):
    best = None
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int16)
    # cost = sum_{i,k} A_ik B_{p_i p_k}
    tot = np.zeros(len(perms))
    for i in range(n):
        for k in range(n):
            if A[i, k] != 0:
                tot += A[i, k] * B[perms[:, i], perms[:, k]]
    if C is not None:
        for i in range(n):
            tot += C[i, perms[:, i]]
    k = int(np.argmin(tot))
    return tot[k], perms[k]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--rows", type=int, default=0)
    ap.add_argument("--cols", type=int, default=0)
    ap.add_argument("--which", default="dnn,tri,lmp1,lmp0,lmp01,rlt2,rlt1lp,rlt2lp")
    ap.add_argument("--solver", default="CLARABEL")
    ap.add_argument("--brute", action="store_true")
    ap.add_argument("--opt", type=float, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    n, A, B, C = read_qaplib(args.dat)
    L = lifted_cost(n, A, B, C)
    regions = grid_regions(n, args.rows, args.cols) if args.rows else None
    res = {"dat": args.dat, "n": n}
    if args.brute:
        v, p = brute(n, A, B, C); res["opt"] = v; res["perm"] = p.tolist()
        print("OPT", v, flush=True)
    elif args.opt is not None:
        res["opt"] = args.opt
    variants = {
        "dnn": (["dnn"], True), "tri": (["tri"], True), "lmp1": (["lmp1"], True),
        "lmp0": (["lmp0"], True), "lmp01": (["lmp0", "lmp1"], True), "rlt2": (["rlt2"], True),
        "tri_lmp1": (["tri", "lmp1"], True),
        "rlt1lp": (["dnn"], False), "rlt2lp": (["rlt2"], False),
    }
    for w in args.which.split(","):
        tags, psd = variants[w]
        if "tri" in tags and regions is None:
            continue
        solver = args.solver if psd else "HIGHS"
        try:
            val, st, t, yv = solve(n, L, tags, regions, solver=solver, psd=psd)
        except Exception as e:
            val, st, t = None, "error: %s" % e, 0
        res[w] = val
        print("%-8s %14.4f  %-12s %7.1fs" % (w, val if val is not None else float("nan"), st, t), flush=True)
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print("RESULT", json.dumps(res))


def reduce_instance(A, B, C, fixes):
    """fix side-0 index i -> side-1 index j for (i, j) in fixes; return the
    free problem (A', B', C', const, free0, free1) with objective
    const + sum A'_ik B'_jl x_ij x_kl + sum C'_ij x_ij over the free indices"""
    n = len(A)
    F0 = [i for i, _ in fixes]; F1 = [j for _, j in fixes]
    free0 = [i for i in range(n) if i not in F0]
    free1 = [j for j in range(n) if j not in F1]
    const = 0.0
    for (i, j) in fixes:
        for (k, l) in fixes:
            const += A[i, k] * B[j, l]
        if C is not None:
            const += C[i, j]
    m = len(free0)
    Cp = np.zeros((m, m))
    for a, i in enumerate(free0):
        for b, j in enumerate(free1):
            v = 0.0 if C is None else C[i, j]
            for (k, l) in fixes:
                v += A[i, k] * B[j, l] + A[k, i] * B[l, j]
            Cp[a, b] = v
    return A[np.ix_(free0, free0)], B[np.ix_(free1, free1)], Cp, const, free0, free1


def grid_regions_free(rows, cols, free0):
    """threshold regions over the free side-0 indices (as in production)"""
    out = []
    for t in range(1, rows):
        S = [a for a, i in enumerate(free0) if i // cols < t]
        if 0 < len(S) < len(free0): out.append(("row<%d" % t, S))
    for t in range(1, cols):
        S = [a for a, i in enumerate(free0) if i % cols < t]
        if 0 < len(S) < len(free0): out.append(("col<%d" % t, S))
    return out
