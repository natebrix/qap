#!/usr/bin/env python3
"""nearcuts.py INSTANCE [ALPHA] [K] [NCUTS]: Jeff Linderoth's refinement (Nov 26 2026) of the cuts of dualcuts.py.
The dual optimum's eigenvector cuts are "too parallel to the objective"; take instead cutting planes from points
somewhat close to the SDP's optimal Y. Here, in-out separation: with Y* the DNN's optimal matrix and Y_LP the
current LP's, separate the point Y* + ALPHA (Y_LP - Y*) (default 0.25; doubled when that point is already PSD),
add a cut v'Yv >= 0 for each of its K (7) most negative eigenvectors, solve the LP again, until NCUTS (126).
Then, as dualcuts.py, keep the cuts and bound every child of row 0: with the first 63 (the number the dual's
cuts needed) and with all of them. Prints the root bound by round and the children's table."""
import sys, time, os, re, subprocess, numpy as np, scipy.sparse as sp, cvxpy as cp, highspy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../branching_lab/lpbound'))
import ajlp

name = sys.argv[1]; ALPHA = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
K = int(sys.argv[3]) if len(sys.argv) > 3 else 7; NCUTS = int(sys.argv[4]) if len(sys.argv) > 4 else 126
A, B, opt, perm = ajlp.inst(name)
Ar, Br, C, s, I, J = ajlp.reduce(A, B, [])
lp = ajlp.LP(Ar, Br, C); m = lp.m; N = m * m + 1; nx, ny = lp.nx, lp.ny; nz = nx + ny
p1 = 1 + lp.yI * m + lp.yJ; p2 = 1 + lp.yK * m + lp.yL
xs = 1 + np.arange(nx)
r = np.concatenate([xs * N, xs, xs * N + xs, p1 * N + p2, p2 * N + p1])
c = np.concatenate([np.arange(nx)] * 3 + [nx + np.arange(ny)] * 2)
T = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(N * N, nz)); e = np.zeros(N * N); e[0] = 1
here = os.path.dirname(os.path.abspath(__file__))

def Ymat(z): return (T @ z + e).reshape(N, N)

def dnn_primal():
    cache = os.path.join(here, 'out', f'{name}_rootprimal.npz')
    if os.path.exists(cache):
        q = np.load(cache); return float(q['z0']), q['z']
    z = cp.Variable(nz, nonneg=True)
    pr = cp.Problem(cp.Minimize(lp.c @ z), [lp.M @ z == lp.b, cp.reshape(T @ z + e, (N, N), order='F') >> 0])
    pr.solve(solver=cp.CLARABEL)
    np.savez_compressed(cache, z0=pr.value + s, z=z.value)
    return pr.value + s, z.value

DNNTEST = os.environ.get('DNNTEST', '/home/user/wt-blas/latest/dnntest')
def dnn_own(fixes):
    o = subprocess.run([DNNTEST, '-u', '1000000', '-T', '1', '-H', '0', '-i', '20000', '-x', ','.join(f'{i}:{j}' for i, j in fixes),
                        f'{ajlp.Q}/qaplib/{name}.dat'], capture_output=True, text=True).stdout
    return float(re.search(r'RESULT .* lb=(\S+) ', o).group(1))

def lpsolve(fixes, cuts):
    if not cuts:
        h = lp.highs({'solver': 'ipm', 'run_crossover': 'off'})
        for i, j in fixes: h.changeColBounds(i * m + j, 1.0, 1.0)
        h.run(); assert h.getModelStatus() == highspy.HighsModelStatus.kOptimal
        return h.getInfo().objective_function_value + s, np.array(h.getSolution().col_value)
    V = np.array(cuts)
    G = np.hstack([2 * V[:, :1] * V[:, 1:] + V[:, 1:] ** 2, 2 * V[:, p1] * V[:, p2]]); rhs = -V[:, 0] ** 2
    z = cp.Variable(nz, nonneg=True)
    pr = cp.Problem(cp.Minimize(lp.c @ z), [lp.M @ z == lp.b, G @ z >= rhs] + [z[i * m + j] == 1 for i, j in fixes])
    pr.solve(solver=cp.CLARABEL); assert pr.status in ('optimal', 'optimal_inaccurate'), pr.status
    return pr.value + s, z.value

z0, zs = dnn_primal()
l0, zl = lpsolve([], [])
print(f"# {name} m={m} LP {l0:.3f} DNN {z0:.3f} opt {opt}; alpha {ALPHA}, {K} cuts a round, {NCUTS} cuts; "
      f"Y* least eigenvalue {np.linalg.eigvalsh(Ymat(zs))[0]:.2e}", flush=True)
cuts = []; rnd = 0; at63 = None
while len(cuts) < NCUTS:
    a = ALPHA
    while True:
        w, U = np.linalg.eigh(Ymat(zs + a * (zl - zs)))
        nneg = int((w < -1e-7).sum())
        if nneg or a >= 1: break
        a = min(1.0, 2 * a)
    if not nneg: print("the LP's Y is PSD: stop"); break
    for q in range(min(K, nneg, NCUTS - len(cuts))): cuts.append(U[:, q])
    b, zl = lpsolve([], cuts); rnd += 1
    print(f"round {rnd:2d} alpha {a:.2f} point: {nneg:3d} negative, least {w[0]:8.4f}; cuts {len(cuts):3d}  bound {b:9.3f}  "
          f"{100*(b-l0)/(z0-l0):5.1f}% of the gap", flush=True)
sets = [(k, cuts[:k]) for k in sorted({min(63, len(cuts)), len(cuts)})]
print("children of row 0: LP | " + " | ".join(f"LP + {k} near cuts" for k, _ in sets) + " | DNN | share of the gap closed", flush=True)
sh = {k: [] for k, _ in sets}
for j in range(m):
    f = [(0, j)]
    a_, _ = lpsolve(f, []); d = dnn_own(f)
    bs = [lpsolve(f, cs)[0] for _, cs in sets]
    for (k, _), b in zip(sets, bs): sh[k].append((b - a_) / (d - a_))
    print(f"0:{j:<2d} {a_:9.3f} | " + " | ".join(f"{b:9.3f}" for b in bs) + f" | {d:9.3f} | " +
          " ".join(f"{100*(b-a_)/(d-a_):5.1f}%" for b in bs) + ("  fathomed by the cuts" if max(bs) > opt - 2 + 1e-6 else ""), flush=True)
for k, _ in sets:
    print(f"{k} near cuts: mean share closed {100*np.mean(sh[k]):.1f}%, least {100*min(sh[k]):.1f}%  (the dual's 63 cuts: mean 58.7%, least 53.3%)")
