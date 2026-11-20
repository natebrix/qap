"""ajlp.py: the Adams-Johnson level-1 RLT LP of a QAP node, built and solved with HiGHS.

Conventions (as landscape/lsq.py, rechecked in t0.py against QAPSOLN): the instance file holds A then B; a
permutation p has p[i] = the location of facility i; cost = sum_{i,k} A[i,k] B[p[i],p[k]]; a fix 'r:v' means
p[r] = v. A node with fixes F is the QAP on free facilities I and free locations J with
  A' = A[I,I], B' = B[J,J], C[i,j] = A[i,i]B[j,j] + sum_{(r,v) in F} (A[i,r]B[j,v] + A[r,i]B[v,j]),
  s = sum_{(r,v),(r',v') in F} A[r,r']B[v,v'],
as the solver forms it (node value = s + sum_{i!=k} A'_ik B'_{p_i p_k} + sum_i C_{i p_i}).

The LP (one variable per unordered pair of compatible assignments, i.e. Y's symmetry used to halve):
  min sum_{i<k, j!=l} (a_ik b_jl + a_ki b_lj) y_{ij,kl} + sum C_ij x_ij
  s.t. sum_{l!=j} y_{ij,kl} = x_ij        for every (i,j), k != i      (R rows)
       sum_{k!=i} y_{ij,kl} = x_ij        for every (i,j), l != j      (C rows)
       sum_j x_ij = 1, sum_i x_ij = 1                                   (assignment rows)
       0 <= x, y <= 1   (the upper bounds are implied; kept so that any dual gives a valid Lagrangian bound)
"""
import numpy as np, scipy.sparse as sp, time, os
import highspy

Q = '/home/user/qap'


def read_qaplib(path):
    t = open(path).read().split()
    n = int(t[0]); v = np.array(t[1:1 + 2 * n * n], dtype=float)
    return n, v[:n * n].reshape(n, n).copy(), v[n * n:].reshape(n, n).copy()


def inst(name):
    n, A, B = read_qaplib(f'{Q}/qaplib/{name}.dat')
    sln = f'{Q}/QAPSOLN/{name}.sln'
    opt, perm = None, None
    if os.path.exists(sln):
        t = open(sln).read().split(); opt = int(t[1]); perm = np.array([int(x) - 1 for x in t[2:2 + n]])
        if abs(cost(A, B, perm) - opt) > 0.5: perm = np.argsort(perm)   # some .sln files (nug12/15/20) list p^-1
        assert abs(cost(A, B, perm) - opt) < 0.5, name
    return A, B, opt, perm


def parse_fixes(s):
    if s in ('-', ''): return []
    return [tuple(int(u) for u in f.split(':')) for f in s.split(',')]


def reduce(A, B, fixes):
    n = len(A); F = parse_fixes(fixes) if isinstance(fixes, str) else list(fixes)
    fr = {r for r, v in F}; fc = {v for r, v in F}
    I = [i for i in range(n) if i not in fr]; J = [j for j in range(n) if j not in fc]
    Ar, Br = A[np.ix_(I, I)].copy(), B[np.ix_(J, J)].copy()
    C = np.outer(np.diag(A)[I], np.diag(B)[J])
    for r, v in F:
        C += np.outer(A[I, r], B[J, v]) + np.outer(A[r, I], B[v, J])
    s = sum(A[r, r2] * B[v, v2] for r, v in F for r2, v2 in F)
    np.fill_diagonal(Ar, 0); np.fill_diagonal(Br, 0)
    return Ar, Br, C, s, I, J


def cost(A, B, p):
    return float((A * B[np.ix_(p, p)]).sum())


class LP:
    """index maps: x_ij -> i*m+j ; y for i<k, j!=l -> m*m + pk(i,k)*(m*(m-1)) + jl(j,l)"""

    def __init__(self, A, B, C):
        m = len(A); self.m = m
        ii, kk = np.triu_indices(m, 1)            # pairs i<k
        jj, ll = np.where(~np.eye(m, dtype=bool))  # ordered j!=l
        P, Q2 = len(ii), len(jj)
        I = np.repeat(ii, Q2); K = np.repeat(kk, Q2); Jv = np.tile(jj, P); L = np.tile(ll, P)
        self.yI, self.yJ, self.yK, self.yL = I, Jv, K, L
        ny = P * Q2; nx = m * m; self.nx, self.ny = nx, ny
        cy = A[I, K] * B[Jv, L] + A[K, I] * B[L, Jv]
        self.c = np.concatenate([C.ravel(), cy])
        # rows: R[(i,j),k] = (i*m+j)*m + k (k != i: rows with k == i are absent; we keep a dense index and drop)
        #       C[(i,j),l] = m^3 + (i*m+j)*m + l
        yidx = nx + np.arange(ny)
        rR1 = (I * m + Jv) * m + K; rR2 = (K * m + L) * m + I
        rC1 = m ** 3 + (I * m + Jv) * m + L; rC2 = m ** 3 + (K * m + L) * m + Jv
        rows = [rR1, rR2, rC1, rC2]; cols = [yidx] * 4; vals = [np.ones(ny)] * 4
        # -x_ij in R[(ij),k] for k != i and C[(ij),l] for l != j
        xi, xj = np.divmod(np.arange(nx), m)
        for k in range(m):
            msk = xi != k
            rows.append((xi[msk] * m + xj[msk]) * m + k); cols.append(np.arange(nx)[msk]); vals.append(-np.ones(msk.sum()))
            msk = xj != k
            rows.append(m ** 3 + (xi[msk] * m + xj[msk]) * m + k); cols.append(np.arange(nx)[msk]); vals.append(-np.ones(msk.sum()))
        # assignment rows
        base = 2 * m ** 3
        rows.append(base + xi); cols.append(np.arange(nx)); vals.append(np.ones(nx))
        rows.append(base + m + xj); cols.append(np.arange(nx)); vals.append(np.ones(nx))
        r = np.concatenate(rows); cc = np.concatenate(cols); v = np.concatenate(vals)
        # compress away the absent rows (k == i, l == j)
        used = np.unique(r); remap = -np.ones(2 * m ** 3 + 2 * m, dtype=np.int64); remap[used] = np.arange(len(used))
        self.rowid = used
        r = remap[r]
        self.M = sp.csc_matrix((v, (r, cc)), shape=(len(used), nx + ny))
        self.b = np.zeros(len(used)); self.b[used >= base] = 1.0
        self.nrow = len(used)

    def highs(self, opts=None):
        h = highspy.Highs(); h.setOptionValue('output_flag', False)
        for k, v in (opts or {}).items(): h.setOptionValue(k, v)
        lp = highspy.HighsLp()
        n = self.nx + self.ny
        lp.num_col_ = n; lp.num_row_ = self.nrow
        lp.col_cost_ = self.c; lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.ones(n)
        lp.row_lower_ = self.b; lp.row_upper_ = self.b
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = self.M.indptr; lp.a_matrix_.index_ = self.M.indices; lp.a_matrix_.value_ = self.M.data
        lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = self.nrow
        h.passModel(lp)
        return h

    def solve(self, solver='ipm', opts=None):
        o = {'solver': solver}
        if solver == 'ipm': o['run_crossover'] = 'off'
        o.update(opts or {})
        h = self.highs(o)
        t = time.process_time(); w = time.time()
        h.run()
        res = dict(status=h.modelStatusToString(h.getModelStatus()), cpu=time.process_time() - t, wall=time.time() - w)
        sol = h.getSolution()
        res['obj'] = h.getInfo().objective_function_value
        res['x'] = np.array(sol.col_value); res['u'] = np.array(sol.row_dual)
        res['iters'] = (h.getInfo().simplex_iteration_count, h.getInfo().ipm_iteration_count, h.getInfo().pdlp_iteration_count)
        return res

    def lagr(self, u):
        """valid bound from any row dual u, all variables in [0,1]: b'u + sum min(0, d); and d itself"""
        d = self.c - self.M.T @ u
        return float(self.b @ u + np.minimum(d, 0).sum()), d

    def child_bounds(self, u):
        """for every child x_ij = 1: a valid bound on the child's LP (hence its QAP), from the dual u alone.
        The child forbids x_il (l != j), x_kj (k != i) and every y that pairs with one of them or that
        pairs (i,j') / (k,j) entries; all those variables are 0 in the child, so their negative reduced
        costs drop out, and x_ij's own reduced cost enters at 1.  Returns an m x m array."""
        m = self.m
        base, d = self.lagr(u)
        neg = np.minimum(d, 0)
        dx = d[:self.nx].reshape(m, m); negx = neg[:self.nx].reshape(m, m)
        # y variables compatible with x_ij=1 keep their min(0,d); incompatible ones (sharing facility i with
        # location != j, or location j with facility != i) are fixed to 0. Their removal can only raise the
        # bound: sum of -min(0,d) over incompatible y. Computed per (i,j) by accumulation.
        ny = neg[self.nx:]
        I, J, K, L = self.yI, self.yJ, self.yK, self.yL
        # incompatible with (i,j): y_{(i,j'),(k,l)} with j' != j (facility i not at j) unless ... ; simplest:
        # a y_{(a,b),(c,e)} survives in child (i,j) iff each of its two assignments is compatible with (i,j):
        # assignment (a,b) compatible iff (a==i) == (b==j).
        # sum over incompatible = total - sum over compatible; compute compatible sums via inclusion-exclusion
        tot = ny.sum()
        S = np.zeros((m, m))
        # for each y, its two assignments (I,J) and (K,L). y is compatible with child (i,j) iff
        # [(i,j)==(I,J) or (i!=I and j!=J)] and same for (K,L). Count via 4 terms per y:
        #   both "generic": i not in {I,K}, j not in {J,L}
        #   (i,j)=(I,J): then (K,L) generic automatically (K!=I, L!=J): yes since I<K and J!=L
        #   (i,j)=(K,L): likewise
        # generic sum: sum over y with i∉{I,K}, j∉{J,L} = tot - rowsum(i) - colsum(j) + both(i,j)
        rowhit = np.zeros(m); colhit = np.zeros(m); both = np.zeros((m, m))
        np.add.at(rowhit, I, ny); np.add.at(rowhit, K, ny)
        np.add.at(colhit, J, ny); np.add.at(colhit, L, ny)
        for a, b in ((I, J), (I, L), (K, J), (K, L)):
            np.add.at(both, (a, b), ny)
        gen = tot - rowhit[:, None] - colhit[None, :] + both
        own = np.zeros((m, m)); np.add.at(own, (I, J), ny); np.add.at(own, (K, L), ny)
        compat_y = gen + own
        # x part: compatible x are x_ij itself (at 1) and x_kl with k!=i, l!=j
        tx = negx.sum(); rx = negx.sum(1); cx = negx.sum(0)
        compat_x = tx - rx[:, None] - cx[None, :] + negx
        val = self.b @ u + compat_y + compat_x + dx
        return val, base


def solve_full_literal(A, B, C, opts=None):
    """the DNN relaxation with R >= 0 dropped, written as in dnn_bound.tex eq. (dnn) with a free symmetric R:
    Y symmetric (N = m^2+1), Y = P Y P for P = V V^T (equivalently T Y = 0, T the assignment map on y = (1,x)),
    Y_00 = 1, gangster entries 0, 0 <= Y <= 1. Variables: the upper triangle of Y. For checking LP's value."""
    m = len(A); N = m * m + 1
    iu, ju = np.triu_indices(N)
    idx = -np.ones((N, N), dtype=np.int64); idx[iu, ju] = np.arange(len(iu)); idx[ju, iu] = idx[iu, ju]
    Lm = np.zeros((N, N))
    for i in range(m):
        for j in range(m):
            for k in range(m):
                for l in range(m):
                    Lm[1 + i * m + j, 1 + k * m + l] = 0.5 * (A[i, k] * B[j, l] + A[k, i] * B[l, j])
    Lm[0, 1:] = Lm[1:, 0] = 0.5 * C.ravel()
    cvec = np.zeros(len(iu)); np.add.at(cvec, idx[iu, ju], np.where(iu == ju, 1.0, 2.0) * Lm[iu, ju])
    # T: rows: sum_j y_(ij) - y_0 = 0 (i), sum_i y_(ij) - y_0 = 0 (j)
    T = np.zeros((2 * m, N))
    for i in range(m):
        T[i, 0] = -1; T[i, 1 + i * m: 1 + i * m + m] = 1
        T[m + i, 0] = -1; T[m + i, 1 + i::m] = 1
    rr, cc, vv = [], [], []; r = 0
    for t in range(2 * m):
        nz = np.nonzero(T[t])[0]
        for col in range(N):          # (T Y)[t, col] = 0
            for a in nz:
                rr.append(r); cc.append(idx[a, col]); vv.append(T[t, a])
            r += 1
    Mx = sp.csc_matrix((vv, (rr, cc)), shape=(r, len(iu)))
    lo = np.zeros(len(iu)); hi = np.ones(len(iu))
    lo[idx[0, 0]] = hi[idx[0, 0]] = 1.0
    for i in range(m):
        for j in range(m):
            for l in range(m):
                if l != j: hi[idx[1 + i * m + j, 1 + i * m + l]] = 0
            for k in range(m):
                if k != i: hi[idx[1 + i * m + j, 1 + k * m + j]] = 0
    h = highspy.Highs(); h.setOptionValue('output_flag', False)
    lp = highspy.HighsLp(); n = len(iu)
    lp.num_col_ = n; lp.num_row_ = r; lp.col_cost_ = cvec; lp.col_lower_ = lo; lp.col_upper_ = hi
    lp.row_lower_ = np.zeros(r); lp.row_upper_ = np.zeros(r)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = Mx.indptr; lp.a_matrix_.index_ = Mx.indices; lp.a_matrix_.value_ = Mx.data
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = r
    h.passModel(lp); h.run()
    return h.getInfo().objective_function_value, h.modelStatusToString(h.getModelStatus())


def glb(A, B, C):
    """Gilmore-Lawler bound of the reduced QAP: LAP on C_ij + min scalar product of row i of A (k != i) and
    row j of B (l != j), the zero diagonals dropped. Symmetric or not: uses sum_k a_ik b_{j p_k}, rows only."""
    from scipy.optimize import linear_sum_assignment
    m = len(A); G = C.copy()
    for i in range(m):
        a = np.sort(np.delete(A[i], i))
        for j in range(m):
            b = np.sort(np.delete(B[j], j))[::-1]
            G[i, j] += a @ b
    r, c = linear_sum_assignment(G)
    return float(G[r, c].sum())


def pdlp_solve(lp, eps=1e-4, threads=4, time_limit=3600, iters=None, verbosity=0, warm=None):
    """OR-tools PDLP on the LP; returns obj, primal, row dual, its valid Lagrangian bound, time, iterations.
    eps: relative optimality tolerance (eps_optimal_relative = eps, absolute 0)."""
    from ortools.pdlp.python import pdlp
    from ortools.pdlp import solvers_pb2
    qp = pdlp.QuadraticProgram()
    n = lp.nx + lp.ny
    qp.objective_vector = lp.c
    qp.constraint_matrix = sp.csr_matrix(lp.M)
    qp.constraint_lower_bounds = lp.b; qp.constraint_upper_bounds = lp.b
    qp.variable_lower_bounds = np.zeros(n); qp.variable_upper_bounds = np.ones(n)
    prm = solvers_pb2.PrimalDualHybridGradientParams()
    prm.num_threads = threads; prm.verbosity_level = verbosity
    tc = prm.termination_criteria
    tc.simple_optimality_criteria.eps_optimal_relative = eps
    tc.simple_optimality_criteria.eps_optimal_absolute = 0.0
    tc.time_sec_limit = time_limit
    if iters: tc.iteration_limit = iters
    t = time.process_time(); w = time.time()
    if warm is not None:
        ini = pdlp.PrimalAndDualSolution(); ini.primal_solution = warm[0]; ini.dual_solution = warm[1]
        r = pdlp.primal_dual_hybrid_gradient(qp, prm, initial_solution=ini)
    else:
        r = pdlp.primal_dual_hybrid_gradient(qp, prm)
    cpu, wall = time.process_time() - t, time.time() - w
    u = np.array(r.dual_solution); x = np.array(r.primal_solution)
    lag, _ = lp.lagr(u)
    return dict(obj=float(lp.c @ x), u=u, x=x, lagr=lag, cpu=cpu, wall=wall,
                iters=r.solve_log.iteration_count, status=str(r.solve_log.termination_reason))
