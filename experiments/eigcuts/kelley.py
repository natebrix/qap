#!/usr/bin/env python3
"""kelley.py INSTANCE [K] [ROUNDS] [SECONDS]: Steve Wright's question (Nov 25 2026). Start from the LP the DNN
relaxation becomes without its PSD condition (the Adams-Johnson level-1 RLT LP, ../branching_lab/lpbound/ajlp.py),
and add eigenvalue cuts v'Yv >= 0 adaptively: each round, eigendecompose the LP's Y (order m^2+1), add a cut for each
of its K most negative eigenvectors, and re-solve with the dual simplex from the last basis.

  Y_00 = 1, Y_0,(ij) = Y_(ij),(ij) = x_ij, Y_(ij),(kl) = y_{ij,kl} (i != k, j != l), 0 elsewhere, so
  v'Yv = v0^2 + sum (2 v0 v_ij + v_ij^2) x_ij + 2 sum_{i<k, j!=l} v_ij v_kl y_{ij,kl}.

Prints a line per round: cuts so far, the LP value (a valid bound), Y's least eigenvalue and count of negative
ones, the cuts tight at the solution, simplex iterations and CPU of the re-solve."""
import sys, time, numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../branching_lab/lpbound'))
import ajlp, highspy

name = sys.argv[1]; K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 100; budget = float(sys.argv[4]) if len(sys.argv) > 4 else 3600
A, B, opt, perm = ajlp.inst(name)
Ar, Br, C, s, I, J = ajlp.reduce(A, B, [])
lp = ajlp.LP(Ar, Br, C); m = lp.m; N = m * m + 1; nx, ny = lp.nx, lp.ny
h = lp.highs({'solver': 'simplex', 'simplex_strategy': 1})
p1 = 1 + lp.yI * m + lp.yJ; p2 = 1 + lp.yK * m + lp.yL; ar = np.arange(nx + ny, dtype=np.int32)
print(f"# {name} m={m} vars={nx+ny} rows={lp.nrow} opt={opt} K={K}", flush=True)
ncut = 0; t0 = time.process_time(); base_rows = lp.nrow
for rnd in range(rounds + 1):
    t = time.process_time(); h.run(); cpu = time.process_time() - t
    info = h.getInfo(); sol = h.getSolution(); z = np.array(sol.col_value); obj = info.objective_function_value + s
    Y = np.zeros((N, N)); Y[0, 0] = 1; Y[0, 1:] = Y[1:, 0] = z[:nx]; Y[np.arange(1, N), np.arange(1, N)] = z[:nx]
    Y[p1, p2] = z[nx:]; Y[p2, p1] = z[nx:]
    w, U = np.linalg.eigh(Y)
    tight = 0
    if ncut:
        ra = np.array(sol.row_value)[base_rows:]; lo = np.array(h.getLp().row_lower_)[base_rows:]
        tight = int((ra - lo < 1e-7).sum())
    print(f"round {rnd:3d} cuts {ncut:5d} bound {obj:12.3f} lmin {w[0]:9.4f} neg {int((w < -1e-7).sum()):4d} tight {tight:5d} "
          f"simplex_it {info.simplex_iteration_count:7d} cpu {cpu:8.1f} total {time.process_time()-t0:8.1f}", flush=True)
    if w[0] > -1e-6 or time.process_time() - t0 > budget: break
    for q in range(min(K, int((w < -1e-7).sum()))):
        v = U[:, q]
        row = np.concatenate([2 * v[0] * v[1:] + v[1:] ** 2, 2 * v[p1] * v[p2]])
        h.addRow(-v[0] ** 2, highspy.kHighsInf, nx + ny, ar, row); ncut += 1
