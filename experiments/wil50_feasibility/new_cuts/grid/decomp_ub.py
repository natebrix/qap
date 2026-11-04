#!/usr/bin/env python3
"""Row / column decomposition of a wil50 node: heuristic optima (simulated
annealing, bin/qapsa) of the full node, the row-only problem (distance
|r_a - r_b|, i.e. facilities into the 5 rows with the node's free capacities)
and the column-only problem (|c_a - c_b|), plus the row and column parts of
the full problem's best permutation.

  python3 decomp_ub.py SPEC [--runs 20] [--sweeps 200] [--out prefix]

OPT_row + OPT_col <= OPT_node always (the node's permutation is feasible for
both); the sum of the two SA values is an upper bound on OPT_row + OPT_col.
If it is below the node's DNN bound, the decoupled bound is useless and a
Lagrangian decomposition can only gain through inequalities of the row and
column hulls that the DNN does not imply."""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl

HERE = os.path.dirname(os.path.abspath(__file__))


def write_sa(path, D, F, C):
    m = len(D)
    with open(path, 'w') as f:
        f.write(f'{m}\n')
        for M in (D, F, C):
            for r in M:
                f.write(' '.join(f'{x:.6g}' for x in r) + '\n')


def run_sa(path, runs, sweeps, seed):
    out = subprocess.run(['taskset', '-c', '2', HERE + '/bin/qapsa', path, str(runs), str(sweeps), str(seed)],
                         capture_output=True, text=True, check=True).stdout
    best = float(out.split()[1]); perm = list(map(int, out.split('PERM')[1].split()))
    return best, perm


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('--runs', type=int, default=20)
    ap.add_argument('--sweeps', type=int, default=200); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tmp', default=HERE + '/runs')
    a = ap.parse_args()
    nd = gl.Node(a.spec)
    tag = a.spec.replace(':', '_').replace(',', '-') or 'root'
    res = {}
    for nm, D, C, c0 in (('full', nd.Dr + nd.Dc, nd.lin_r + nd.lin_c, nd.const_r + nd.const_c),
                         ('row', nd.Dr, nd.lin_r, nd.const_r), ('col', nd.Dc, nd.lin_c, nd.const_c)):
        p = f'{a.tmp}/sa_{tag}_{nm}.txt'
        write_sa(p, D, nd.F, C)
        best, perm = run_sa(p, a.runs, a.sweeps, a.seed)
        res[nm] = (best + c0, perm)
        extra = ''
        if nm == 'full':
            extra = f"  (row part {nd.perm_cost(perm, 'row'):.0f}, col part {nd.perm_cost(perm, 'col'):.0f})"
        print(f'{nm:5s} SA best {best + c0:.0f}{extra}', flush=True)
    s = res['row'][0] + res['col'][0]
    print(f'row + col SA optima {s:.0f}; full SA {res["full"][0]:.0f}; coupling loss >= {res["full"][0] - s:.0f}')


if __name__ == '__main__':
    main()
