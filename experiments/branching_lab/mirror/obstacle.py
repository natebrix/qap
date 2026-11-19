"""The completions that hold the exchange cuts' ceiling down at a node: lagts at the ceiling's mu on pair (I,K),
each found completion's cost, its Delta for every fixed pair (max shown), and its best improving exchange of any
two rows (fixed-fixed, fixed-free or free-free). usage: obstacle.py INSTANCE FIXES I K MU CUTOFF"""
import sys, subprocess, itertools, numpy as np
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, cost
inst, fixes, I, K, mu, cutoff = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], float(sys.argv[6])
A, B = read_dat(inst); n = len(A); fx = dict(tuple(map(int, t.split(':'))) for t in fixes.split(','))
out = subprocess.run(['./lagts', f'/home/user/wt-cognac/qaplib/{inst}.dat', fixes, f'{I}:{K}:{mu}', '3000000', '6', '5'], capture_output=True, text=True).stdout
seen = set()
for ln in sorted([l for l in out.splitlines() if l.startswith('sol')], key=lambda l: float(l.split()[2])):
    q = np.array([int(x) for x in ln.split(' p ')[1].split(',')]); 
    if tuple(q) in seen: continue
    seen.add(tuple(q)); c = cost(A, B, q)
    best = (0, None)
    fixed_pair_max = -1e18
    for i, k in itertools.combinations(range(n), 2):
        q2 = q.copy(); q2[[i, k]] = q2[[k, i]]; g = cost(A, B, q2) - c   # < 0: the exchange improves q
        if i in fx and k in fx: fixed_pair_max = max(fixed_pair_max, -g)
        if g < best[0]: best = (g, (i, k))
    kind = '-' if best[1] is None else ('fixed-fixed' if best[1][0] in fx and best[1][1] in fx else ('free-free' if best[1][0] not in fx and best[1][1] not in fx else 'fixed-free'))
    print(f'cost - cutoff {c - cutoff:+,.0f}  max over fixed pairs of Delta {fixed_pair_max:+,.0f}  best exchange {best[1]} gains {-best[0]:,} ({kind})')
