"""Ceiling of ALL fixed-pair exchange cuts together at a node: any convex relaxation of the node with every cut
Delta_ik <= 0 (i, k fixed rows) added, or any Lagrangian bound min cost + sum mu_ik Delta_ik, certifies at most
    max_{mu >= 0} min_{p in node} cost(p) + sum mu_ik Delta_ik(p)  <=  the same max with p restricted to a set S.
S is grown by column generation: solve the LP max_{mu, t} t s.t. t <= cost(s) + sum mu Delta(s) for s in S
(mu <= MUMAX), run lagts (simulated annealing) on cost + sum mu Delta at the LP's mu, add what it finds, repeat.
Every LP value is a valid upper bound on what the cuts can certify; the last is the tightest found.
usage: ceiling_multi.py INSTANCE FIXES CUTOFF ROUNDS ITERS RESTARTS [MUMAX=1e4]"""
import sys, subprocess, itertools, numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, read_sln, cost
inst, fixes, cutoff, rounds, iters, rs = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
MUMAX = float(sys.argv[7]) if len(sys.argv) > 7 else 1e4
A, B = read_dat(inst); opt, p = read_sln(inst)
fx = dict(tuple(map(int, t.split(':'))) for t in fixes.split(','))
pairs = list(itertools.combinations(sorted(fx), 2))
S = {}
def add(q):
    q = np.asarray(q); key = tuple(q)
    if key in S: return
    c = cost(A, B, q); ds = []
    for i, k in pairs:
        q2 = q.copy(); q2[[i, k]] = q2[[k, i]]; ds.append(c - cost(A, B, q2))
    S[key] = (c, np.array(ds, float))
for i, k in pairs:   # the optimum's images under one exchange of two fixed rows that land in the node
    q = p.copy(); q[[i, k]] = q[[k, i]]
    if all(q[r] == c for r, c in fx.items()): add(q)
def sa(mu, seed):
    spec = ','.join(f'{i}:{k}:{m:.6g}' for (i, k), m in zip(pairs, mu) if m > 0) or '-'
    out = subprocess.run(['./lagts', f'/home/user/wt-cognac/qaplib/{inst}.dat', fixes, spec, iters, rs, str(seed)], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        if ln.startswith('sol'): add([int(x) for x in ln.split(' p ')[1].split(',')])
sa(np.zeros(len(pairs)), 1)
print(f'{inst} node {fixes} (depth {len(fx)}, {len(pairs)} pairs), cutoff {cutoff:.2f}')
for r in range(rounds):
    C = np.array([v[0] for v in S.values()], float); D = np.array([v[1] for v in S.values()])
    P = len(pairs)
    # variables x = (t, mu_1..mu_P): maximise t; t - D mu <= C
    res = linprog(np.r_[-1.0, np.zeros(P)], A_ub=np.c_[np.ones(len(C)), -D], b_ub=C,
                  bounds=[(None, None)] + [(0, MUMAX)] * P, method='highs')
    t, mu = res.x[0], res.x[1:]
    before = len(S); sa(mu, 100 + r)
    C2 = np.array([v[0] for v in S.values()], float); D2 = np.array([v[1] for v in S.values()])
    found = (C2 + D2 @ mu).min()
    act = ', '.join(f'({i},{k}) {m:.3g}' for (i, k), m in zip(pairs, mu) if m > 1e-9)
    print(f'round {r}: LP ceiling cutoff {t - cutoff:+,.0f}; SA at its mu finds cutoff {found - cutoff:+,.0f}; |S| {before} -> {len(S)}; mu: {act}')
    if found >= t - 1: break
