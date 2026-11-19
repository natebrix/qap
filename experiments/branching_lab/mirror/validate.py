"""Brute-force validation of the QAP_MIRROR hook (mirror_hook.h) on small instances.
For every search of validate.sh's list: bnbdrv.hook with QAP_MIRROR (several pairs, one mu) and QAP_MIRROR_CHECK=1.
Checks, per instance run:
  1. every 'mirror-check' line: the reduced objective equals cost + sum mu Delta on every completion (m <= 8)
     or 500 random ones, and the D0 + sum D_l formula equals Delta (max differences must be 0);
  2. every 'mirror:' node with m <= 9: the certified values (QPB, DNN first solve, cut phase) never exceed
     the brute-force minimum over the node's completions of cost + sum mu Delta;
  3. that minimum is <= the minimum cost over the completions with every applied Delta <= 0 (the claim);
  4. the search ends at the optimum (brute force over all n! permutations, or QAPSOLN), from optimum + 1.
usage: validate.py  (reads runs/val_*/err.txt, log.txt and cmd.txt; writes to stdout)"""
import sys, glob, os, re, itertools, numpy as np
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_dat, read_sln

H = os.path.dirname(os.path.abspath(__file__))

def all_costs(A, B, P):
    """P: k x n array of permutations -> costs"""
    return np.einsum('uv,kuv->k', A, B[P[:, :, None], P[:, None, :]])

def node_eval(A, B, fixes, pairs):
    n = len(A); fr = [r for r in range(n) if r not in fixes]; used = set(fixes.values())
    fc = [c for c in range(n) if c not in used]
    perms = np.array(list(itertools.permutations(fc)), dtype=np.int64)
    P = np.empty((len(perms), n), dtype=np.int64)
    for r, c in fixes.items(): P[:, r] = c
    P[:, fr] = perms
    cost = all_costs(A, B, P).astype(np.float64)
    f = cost.copy(); ok = np.ones(len(P), bool)
    for i, k, mu in pairs:
        Q = P.copy(); Q[:, [i, k]] = Q[:, [k, i]]
        d = cost - all_costs(A, B, Q)
        f += mu * d; ok &= d <= 0
    return f.min(), (cost[ok].min() if ok.any() else np.inf), cost.min()

tot = dict(runs=0, checks=0, check_bad=0, nodes=0, viol=0, claim_bad=0, wrong_end=0)
for d in sorted(glob.glob(f'{H}/runs/val_*')):
    cmd = open(f'{d}/cmd.txt').read().split('\n')[0]
    inst = re.search(r'/([a-z0-9]+)\.dat', cmd).group(1)
    A, B = read_dat(inst); n = len(A)
    try: opt = read_sln(inst)[0]
    except Exception: opt = min(int((A * B[np.ix_(p, p)]).sum()) for p in map(list, itertools.permutations(range(n))))
    err = open(f'{d}/err.txt').read().splitlines(); log = open(f'{d}/log.txt').read()
    tot['runs'] += 1
    ends = re.findall(r'New objective value: ([0-9.e+]+)', log)
    end = float(ends[-1]) if ends else None
    if end is None or abs(end - opt) > 0.5: tot['wrong_end'] += 1
    worst = 0; nodes = 0; viol = 0; claim = 0; maxexcess = -np.inf; tight = 0
    for ln in err:
        if ln.startswith('mirror-check'):
            t = ln.split(); a, b = float(t[13]), float(t[17]); tot['checks'] += 1
            if abs(a) > 1e-6 or abs(b) > 1e-6: tot['check_bad'] += 1   # double rounding of ~4e7 values: ~1e-9
        if ln.startswith('mirror: '):
            head, *rest = ln.split(' | ')
            fx = {int(s.split(':')[0]): int(s.split(':')[1]) for s in head.split()[4:]}
            if n - len(fx) > 8: continue
            pairs = []
            for part in rest:
                if part.startswith('pair'):
                    t = part.split(); i, k = map(int, t[1].split(',')); pairs.append((i, k, float(t[3])))
            vals = re.search(r'qpb (\S+) dnn (\S+) cut (\S+)', ln)
            vs = [float(x) for x in vals.groups() if x not in ('nan', '-nan')]
            fmin, cmin_ok, cmin = node_eval(A, B, fx, pairs)
            nodes += 1
            for v in vs:
                maxexcess = max(maxexcess, v - fmin)
                if v > fmin + 1e-6: viol += 1
            if fmin > cmin_ok + 1e-6: claim += 1
            if max(vs, default=-np.inf) > cmin + 1e-6: tight += 1
    tot['nodes'] += nodes; tot['viol'] += viol; tot['claim_bad'] += claim
    print(f'{os.path.basename(d):28s} {inst:7s} opt {opt:>10} end {end} | nodes checked {nodes:4d}  certified > min(cost+muDelta): {viol}'
          f'  max(certified - min) {maxexcess:.4f}  claim failures {claim}  nodes whose bound passed their own optimum {tight}')
print(tot)
