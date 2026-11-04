#!/usr/bin/env python3
"""Validity check of the objective-part (Lagrangian decomposition) bound on
small grid instances, against enumeration:

   bound = DNN(node of A' = s*D - k*D_K) / s + (k/s) * OPT_K(node)

where D is the Manhattan distance of an R x C grid, D_K the sum of the cut
matrices of a threshold subset K, OPT_K(node) the exact (enumerated) optimum
of the K-part alone at the node. Every certified value must be <= the node's
optimum. Also checks the plumbing of gridlib.Node (row/col split, linear
terms) on these instances.

  python3 validate_objpart.py [-N 10] [--seed 1]
"""
import argparse, itertools, os, random, re, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl

DNNTEST = gl.ROOT + '/latest/dnntest'


def node_opt(D, F, fixes):
    """min over permutations (location a -> facility p[a]) with fixes, of
    sum_{a,b} D[a,b] F[p[a],p[b]]"""
    n = len(D); fl = dict(fixes); used_f = set(fl.values())
    free_loc = [a for a in range(n) if a not in fl]; free_fac = [i for i in range(n) if i not in used_f]
    best = None
    for perm in itertools.permutations(free_fac):
        p = dict(fl); p.update(zip(free_loc, perm))
        pv = np.array([p[a] for a in range(n)])
        v = (D * F[np.ix_(pv, pv)]).sum()
        best = v if best is None or v < best else best
    return best


def thr(R, C, name):
    n = R * C; t = int(name[4:])
    ind = np.array([((a // C) if name.startswith('row') else (a % C)) < t for a in range(n)])
    return np.abs(ind[:, None].astype(float) - ind[None, :].astype(float))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('-N', type=int, default=10); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed); tmp = tempfile.mkdtemp(dir=gl.ROOT + '/experiments/wil50_feasibility/new_cuts/grid/runs')
    checks = bad = 0; gains = []
    for inst in range(a.N):
        R, C = rng.choice([(2, 4), (2, 3), (3, 3)])
        n = R * C
        Dr, Dc = gl.grid_dist(R, C); D = Dr + Dc
        F = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                F[i, j] = F[j, i] = rng.randint(0, 9)
        names = [f'row<{t}' for t in range(1, R)] + [f'col<{t}' for t in range(1, C)]
        K = rng.sample(names, rng.randint(1, len(names)))
        DK = sum(thr(R, C, nm) for nm in K)
        nfix = rng.choice([0, 1, 2])
        locs = rng.sample(range(n), nfix); facs = rng.sample(range(n), nfix)
        fixes = list(zip(locs, facs))
        spec = ','.join(f'{x}:{y}' for x, y in fixes) or '-'
        opt = node_opt(D, F, fixes); optK = node_opt(DK, F, fixes)
        # plumbing: Node's row/col split of the optimal value via a brute-force permutation
        vals = {}
        for s_, k_ in ((1, 0), (2, 1), (4, 1), (4, 3)):
            Ap = s_ * D - k_ * DK
            f = f'{tmp}/v{inst}_{s_}_{k_}.dat'; gl.write_qaplib(f, Ap, F)
            out = subprocess.run(['taskset', '-c', '2', DNNTEST, '-x', spec, '-c', '1e12', '-T', '1', '-i', '3000', '-H', '0', f],
                                 capture_output=True, text=True).stdout
            m = re.search(r'RESULT .* lb=(\S+)', out)
            if not m:
                print('no result', out[-300:]); continue
            lb = float(m.group(1)) / s_ + (k_ / s_) * optK
            checks += 1
            if lb > opt + 1e-6:
                bad += 1; print(f'  INVALID: inst {inst} {R}x{C} K={K} node {spec} s={s_} k={k_}: bound {lb:.4f} > opt {opt}')
            vals[(s_, k_)] = lb
        print(f'inst {inst}: {R}x{C} K={"+".join(K)} node {spec} opt {opt:.0f} optK {optK:.0f} DNN {vals.get((1,0),0):.3f} '
              f'alpha .5 {vals.get((2,1),0):.3f} .25 {vals.get((4,1),0):.3f} .75 {vals.get((4,3),0):.3f}', flush=True)
    print(f'{checks} certified values, {bad} above the node optimum')


if __name__ == '__main__':
    main()
