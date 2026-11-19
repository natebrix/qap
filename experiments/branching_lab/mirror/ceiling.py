"""Ceiling of the exchange cut at a node: any convex relaxation (DNN included) of the node with the cut Delta_ik <= 0
added, or any Lagrangian bound min cost + mu Delta, certifies at most
    max_mu min_{p in node} cost(p) + mu Delta(p)   (= min cost over conv(node) with Delta <= 0, LP duality)
<=  max_mu min_{p in S} cost(p) + mu Delta(p)       for any set S of completions.
S: the optimum's mirror (cost and Delta computed here) and every restart's best of lagts (simulated annealing on
cost + mu Delta) over a ladder of mu. usage: ceiling.py INSTANCE FIXES I K CUTOFF ITERS RESTARTS MU1,MU2,...
prints S's points, the ceiling (a valid upper bound on what the cut can reach), the mu that attains it."""
import sys, subprocess, numpy as np
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, read_sln, cost
inst, fixes, i, k, cutoff, iters, rs, mus = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6], sys.argv[7], sys.argv[8]
A, B = read_dat(inst); opt, p = read_sln(inst)
fx = dict(tuple(map(int, t.split(':'))) for t in fixes.split(','))
S = []
def add(q, tag):
    q = np.asarray(q); q2 = q.copy(); q2[[i, k]] = q2[[k, i]]; c = cost(A, B, q); S.append((c, c - cost(A, B, q2), tag))
# the optimum's image in the node, if the node is its exchange
q = p.copy(); q[[i, k]] = q[[k, i]]
if all(q[r] == c for r, c in fx.items()): add(q, 'mirror of the optimum')
for mu in mus.split(','):
    out = subprocess.run(['./lagts', f'/home/user/wt-cognac/qaplib/{inst}.dat', fixes, f'{i}:{k}:{mu}', iters, rs, '1'], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        if ln.startswith('sol'):
            add([int(x) for x in ln.split(' p ')[1].split(',')], f'SA mu {mu}')
S = sorted(set(S))
print(f'{inst} node {fixes}, pair ({i},{k}), cutoff {cutoff:.2f}; {len(S)} completions (cost - cutoff, Delta):')
for c, d, t in S[:12]: print(f'   {c - cutoff:>14,.0f} {d:>12,.0f}   {t}')
neg = [s for s in S if s[1] <= 0]
print(f'   ... cheapest completion with Delta <= 0: ' + (f'{neg[0][0] - cutoff:+,.0f} (Delta {neg[0][1]:,})' if neg else 'none found'))
C = np.array([s[0] for s in S], float); D = np.array([s[1] for s in S], float)
grid = np.concatenate([[0], np.geomspace(1e-4, 1e4, 4001)])
H = np.array([(C + m * D).min() for m in grid])
j = H.argmax()
print(f'ceiling: max over mu of min over these completions of cost + mu Delta = cutoff {H[j] - cutoff:+,.0f} at mu {grid[j]:.4g}')
for m in [0, 1, 2, 4, 7, 8, 16, 32, 64, 128, 256, 1024]:
    print(f'   mu {m:>5}: min over S of cost + mu Delta - cutoff = {(C + m * D).min() - cutoff:+,.0f}')
