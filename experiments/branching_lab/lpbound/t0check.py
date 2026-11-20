"""t0check.py: (1) conventions: the .sln's cost; (2) the compact AJ LP equals the literal 'DNN without R >= 0'
LP on small nodes; (3) child bounds from the dual are valid (<= the child's own LP) on a small node."""
import sys, numpy as np, itertools, ajlp
for name in ['nug12', 'had12', 'tai12b']:
    A, B, opt, p = ajlp.inst(name); print(name, 'sln cost', ajlp.cost(A, B, p), 'opt', opt)
rng = np.random.default_rng(1)
for name, fixes in [('nug12', '0:3,5:7,2:1'), ('had12', '1:2,3:3,7:9'), ('tai12b', '2:5,8:1')]:
    A, B, opt, p = ajlp.inst(name)
    Ar, Br, C, s, I, J = ajlp.reduce(A, B, fixes)
    lp = ajlp.LP(Ar, Br, C); r = lp.solve('simplex')
    lit = ajlp.solve_full_literal(Ar, Br, C)
    # brute-force node optimum
    m = len(I); best = min(s + ajlp.cost(Ar, Br, list(q)) + C[range(m), list(q)].sum() for q in itertools.permutations(range(m)))
    print(name, fixes, 'm', m, 'AJ LP', round(r['obj'] + s, 4), 'literal', round(lit[0] + s, 4), lit[1], 'node opt', best, 'GLB', round(ajlp.glb(Ar, Br, C) + s, 2))
    cb, base = lp.child_bounds(r['u'])
    worst = 1e9
    for i in range(m):
        for j in range(m):
            Fc = ajlp.parse_fixes(fixes) + [(I[i], J[j])]
            A2, B2, C2, s2, _, _ = ajlp.reduce(A, B, Fc)
            v = ajlp.LP(A2, B2, C2).solve('simplex')['obj'] + s2
            worst = min(worst, v - (cb[i, j] + s))
    print('  lagrangian base', round(base + s, 4), ' min over children of (child LP - dual child bound)', round(worst, 6))
