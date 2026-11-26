#!/usr/bin/env python3
"""dualcuts.py INSTANCE [ROW]: Jeff Linderoth's variant (Nov 26 2026; Gunluk, Junger, Linderoth, Lodi, Luedtke,
"Sparse Cuts for the Positive Semidefinite Cone", arXiv 2603.09864, section 4.3's "dense dual SDP cuts"). Take the
eigenvector cuts v'Yv >= 0 once, from the eigenvectors of the dual PSD matrix S* of the node's SDP (here the DNN:
the Adams-Johnson RLT1 LP of ../branching_lab/lpbound/ajlp.py with Y(x, y) PSD); with all of S*'s positive
eigenvectors the LP has the SDP's value at that node. Then branch and keep the cuts: how much of the gap between
the child's LP and the child's DNN bound do the root's cuts still close?

Prints: the root LP, DNN and LP + cuts (by number of cuts kept, largest eigenvalues first); then for every child
of ROW (default 0), and for every child of the next row under the weakest of them: LP, LP + root cuts, DNN, the
share of the gap closed, CPU of each."""
import sys, time, os, re, subprocess, numpy as np, scipy.sparse as sp, cvxpy as cp, highspy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../branching_lab/lpbound'))
import ajlp

name = sys.argv[1]; ROW = int(sys.argv[2]) if len(sys.argv) > 2 else 0
A, B, opt, perm = ajlp.inst(name)
Ar, Br, C, s, I, J = ajlp.reduce(A, B, [])
lp = ajlp.LP(Ar, Br, C); m = lp.m; N = m * m + 1; nx, ny = lp.nx, lp.ny; nz = nx + ny
p1 = 1 + lp.yI * m + lp.yJ; p2 = 1 + lp.yK * m + lp.yL
# vec(Y) = T z + e (column-major)
xs = 1 + np.arange(nx)
r = np.concatenate([xs * N, xs, xs * N + xs, p1 * N + p2, p2 * N + p1])
c = np.concatenate([np.arange(nx)] * 3 + [nx + np.arange(ny)] * 2)
T = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(N * N, nz)); e = np.zeros(N * N); e[0] = 1

def dnn(fixes):
    z = cp.Variable(nz, nonneg=True)
    psd = cp.reshape(T @ z + e, (N, N), order='F') >> 0
    cons = [lp.M @ z == lp.b, psd] + [z[i * m + j] == 1 for i, j in fixes]
    pr = cp.Problem(cp.Minimize(lp.c @ z), cons)
    t = time.process_time(); pr.solve(solver=cp.CLARABEL); cpu = time.process_time() - t
    S = psd.dual_value; S = (S + S.T) / 2
    return pr.value + s, S, cpu

DNNTEST = os.environ.get('DNNTEST', '/home/user/wt-blas/latest/dnntest')
def dnn_own(fixes):
    """the child's DNN bound from latest/dnntest, run to convergence (a cutoff out of reach)"""
    t = time.process_time()
    o = subprocess.run([DNNTEST, '-u', '1000000', '-T', '1', '-H', '0', '-i', '20000', '-x', ','.join(f'{i}:{j}' for i, j in fixes),
                        f'{ajlp.Q}/qaplib/{name}.dat'], capture_output=True, text=True).stdout
    mm = re.search(r'RESULT .* lb=(\S+) .* cpu=(\S+)', o)
    return float(mm.group(1)), None, float(mm.group(2))

def lpsolve(fixes, cuts):
    """the LP, with the cuts as dense rows. HiGHS's IPM stops short of optimal once there are some 40 of them and its
    dual simplex did not finish in 50 minutes, so the LP with cuts goes to Clarabel (an interior-point conic solver)."""
    if not cuts:
        h = lp.highs({'solver': 'ipm', 'run_crossover': 'off'})
        for i, j in fixes: h.changeColBounds(i * m + j, 1.0, 1.0)
        t = time.process_time(); h.run(); cpu = time.process_time() - t
        assert h.getModelStatus() == highspy.HighsModelStatus.kOptimal
        return h.getInfo().objective_function_value + s, cpu
    V = np.array(cuts)                                           # k x N
    G = np.hstack([2 * V[:, :1] * V[:, 1:] + V[:, 1:] ** 2, 2 * V[:, p1] * V[:, p2]]); rhs = -V[:, 0] ** 2
    z = cp.Variable(nz, nonneg=True)
    pr = cp.Problem(cp.Minimize(lp.c @ z), [lp.M @ z == lp.b, G @ z >= rhs] + [z[i * m + j] == 1 for i, j in fixes])
    t = time.process_time(); pr.solve(solver=cp.CLARABEL); cpu = time.process_time() - t
    assert pr.status in ('optimal', 'optimal_inaccurate'), pr.status
    return pr.value + s, cpu

print(f"# {name} m={m} order {N} vars {nz} rows {lp.nrow} opt {opt}", flush=True)
cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', f'{name}_rootdual.npz')
if os.path.exists(cache):
    q = np.load(cache); z0, S, tdnn = float(q['z0']), q['S'], float(q['cpu'])
else:
    z0, S, tdnn = dnn([]); np.savez_compressed(cache, z0=z0, S=S, cpu=tdnn)
l0, tlp = lpsolve([], [])
w, U = np.linalg.eigh(S); w, U = w[::-1], U[:, ::-1]
rank = int((w > 1e-6 * w[0]).sum())
print(f"root: LP {l0:.3f} ({tlp:.1f} s)  DNN {z0:.3f} ({tdnn:.1f} s)  dual S*: {rank} eigenvalues above 1e-6 of the largest, "
      f"{int((w > 1e-3 * w[0]).sum())} above 1e-3, {int((w > 1e-2 * w[0]).sum())} above 1e-2", flush=True)
for k in sorted({10, 20, 40, 80, rank, N}):
    if k > N: continue
    b, t = lpsolve([], [U[:, q] for q in range(k)])
    print(f"root LP + {k:3d} cuts: {b:10.3f}  {100*(b-l0)/(z0-l0):5.1f}% of the gap  {t:6.1f} s", flush=True)
cuts = [U[:, q] for q in range(rank)]

def children(fixes, row, tag):
    out = []
    used = {j for _, j in fixes}
    for j in range(m):
        if j in used: continue
        f = fixes + [(row, j)]
        a, ta = lpsolve(f, []); b, tb = lpsolve(f, cuts); d, S2, td = dnn_own(f)
        sh = 100 * (b - a) / (d - a) if d - a > 1e-6 else float('nan')
        print(f"{tag} {row}:{j:<2d} LP {a:9.3f}  LP+root cuts {b:9.3f}  DNN {d:9.3f}  closed {sh:5.1f}%  "
              f"{'DNN fathoms' if d > opt - 1 + 1e-6 else ''}{'  cuts fathom' if b > opt - 1 + 1e-6 else ''}  cpu {ta:.1f} {tb:.1f} {td:.1f}", flush=True)
        out.append((d, j, a, b))
    return out
ch = children([], ROW, 'depth 1')
g = [(b - a) / (d - a) for d, j, a, b in ch if d - a > 1e-6]
print(f"depth 1: mean share closed {100*np.mean(g):.1f}%, least {100*min(g):.1f}%", flush=True)
d, j, a, b = min(ch)
row2 = (ROW + 1) % m
ch2 = children([(ROW, j)], row2, f'depth 2 under {ROW}:{j}')
g = [(b - a) / (d - a) for d, j, a, b in ch2 if d - a > 1e-6]
print(f"depth 2: mean share closed {100*np.mean(g):.1f}%, least {100*min(g):.1f}%", flush=True)
