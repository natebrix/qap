"""tsmall.py NAME FIXES [DNN_ITERS]: T1/T2 with exact duals on a small node. Parent: the AJ LP by simplex (a vertex
dual), by IPM without crossover (an interior dual) and by PDLP to 1e-4. Every child x_ij = 1: its valid bound from
each dual (ajlp.child_bounds), its own exact LP (IPM), its QPB (bnbdrv, default.par) and its DNN (dnntest -H 0,
DNN_ITERS iterations, default 4,000, no cutoff). Writes out/small_NAME_FIXES.tsv (one row per child) and appends
the summary to out/small_summary.txt."""
import sys, os, numpy as np, ajlp, t0
from scipy.stats import spearmanr
name, fx = sys.argv[1], sys.argv[2]; it = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
t0.name = name
A, B, opt, perm = ajlp.inst(name)
Ar, Br, C, s, I, J = ajlp.reduce(A, B, fx)
lp = ajlp.LP(Ar, Br, C); m = lp.m
duals = {}
rs = lp.solve('simplex'); duals['simplex'] = rs['u']
ri = lp.solve('ipm'); duals['ipm'] = ri['u']
rp = ajlp.pdlp_solve(lp, eps=1e-4, threads=1); duals['pdlp'] = rp['u']
cb = {k: lp.child_bounds(u)[0] + s for k, u in duals.items()}
LPv = rs['obj'] + s
x = rs['x'][:lp.nx].reshape(m, m); xi = ri['x'][:lp.nx].reshape(m, m)
F = ajlp.parse_fixes(fx) if fx != '-' else []
tag = f'{name}_{fx.replace(":", "_").replace(",", "-")}'
rows = []
for i in range(m):
    for j in range(m):
        cf = (fx + ',' if fx != '-' else '') + f'{I[i]}:{J[j]}'
        A2, B2, C2, s2, _, _ = ajlp.reduce(A, B, cf)
        clp = ajlp.LP(A2, B2, C2).solve('ipm')['obj'] + s2 if m > 2 else np.nan
        q = t0.qpb(cf); d = t0.dnn(cf, it)
        rows.append((I[i], J[j], x[i, j], xi[i, j], cb['simplex'][i, j], cb['ipm'][i, j], cb['pdlp'][i, j], clp, q, d[0]))
R = np.array(rows, dtype=float)
with open(f'out/small_{tag}.tsv', 'w') as f:
    f.write('row\tcol\tx_simplex\tx_ipm\tcb_simplex\tcb_ipm\tcb_pdlp\tchild_LP\tchild_QPB\tchild_DNN\n')
    for r in rows: f.write('\t'.join(f'{v:.6g}' if isinstance(v, float) else str(v) for v in r) + '\n')
dn = R[:, 9]
lines = [f'== {name} {fx} m={m}: parent LP {LPv:.2f} (simplex {rs["cpu"]:.1f}s, ipm {ri["cpu"]:.1f}s, pdlp {rp["cpu"]:.1f}s '
         f'dual bound {rp["lagr"] + s:.2f}); parent QPB {t0.qpb(fx):.2f} DNN {t0.dnn(fx, it)[0]:.2f}; opt {opt}',
         f'  LP x support: simplex {(x > 1e-6).sum()} of {m*m} > 0; ipm {(xi > 1e-6).sum()}',
         ]
for k, col in (('simplex', 4), ('ipm', 5), ('pdlp', 6)):
    g = R[:, col] - LPv
    lines.append(f'  {k} dual: child gain median {np.median(g):.2f} p90 {np.percentile(g, 90):.2f} max {g.max():.2f}; '
                 f'gain > 0 at {(g > 1e-6 * max(1, abs(LPv))).sum()} children; Spearman(cb, child DNN) {spearmanr(R[:, col], dn)[0]:.3f}')
lines.append(f'  Spearman with child DNN: child LP {spearmanr(R[:, 7], dn)[0]:.3f}, child QPB {spearmanr(R[:, 8], dn)[0]:.3f}')
lines.append(f'  child DNN - child LP: median {np.median(dn - R[:, 7]):.2f}; child LP - child QPB: median {np.median(R[:, 7] - R[:, 8]):.2f}')
# per line: Spearman within lines (rows and columns), median over lines
for lab, col in (('cb_simplex', 4), ('cb_ipm', 5), ('child_LP', 7), ('child_QPB', 8)):
    rho = []
    for k in range(m):
        for sel in (R[:, 0] == I[k], R[:, 1] == J[k]):
            rho.append(spearmanr(R[sel, col], dn[sel])[0])
    lines.append(f'  within-line Spearman with child DNN, median over {len(rho)} lines: {lab} {np.nanmedian(rho):.3f}')
print('\n'.join(lines), flush=True)
open('out/small_summary.txt', 'a').write('\n'.join(lines) + '\n')
