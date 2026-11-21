#!/usr/bin/env python3
"""budget.py: what the DNN lookahead costs at the top of each tree, as a share of the run. No solver runs: the
inputs are read off the estimates and runs already in the repository (sources in README.md).

A firing at a node with m free facilities bounds the children of K rows and K columns: 2*K*m - K*K solves (a
cell on a chosen row and a chosen column is solved once), each of at most CAP warm ADMM iterations at
dimension m - 1. Where symmetry is left only rows are scored (K*m solves); the full count is used here, and
every child is charged all CAP iterations (at tho30's depth 1 they averaged 110-137 of 150), so the costs
are upper bounds.   budget.py [K=4] [CAP=150]"""
import sys
K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 150
# instance: n, the run in core-hours (low, high) on MACHINE, and per depth of firing: (branching nodes, seconds
# per ADMM iteration of a CHILD of such a node on MACHINE, how the count was got)
I = {
 'tho40 (settings c)': dict(n=40, run=(5182, 6255), machine='c7a core-hours, estimate of Nov 19', depths={
     1: (12, 0.161 * 1.29, 'every node'),
     2: (403, 0.135 * 1.29, 'every node (settings b; c: 404 nodes, the split not recorded)'),
     3: (round(13459 * 2292 / 4737), 0.135 * 1.29 * (36 / 37) ** 5.2, 'estimate x a sample of open nodes')}),
 'wil50': dict(n=50, run=(7500, 13200), machine='cloud Xeon core-hours, fit of Nov 11', depths={
     1: (15, 0.37, 'every node'),
     2: (304, 0.32, 'every node (cut.zone 0.005; fewer under 0.0075)'),
     3: (round(14192 * 103 / 400), 0.32 * (46 / 47) ** 5.2, 'frontier x a sample of 400, no cut phase')}),
 'sko49': dict(n=49, run=(37000, 41000), machine='cloud Xeon core-hours, estimate Cw of Nov 13', depths={
     1: (10, 0.383, 'every node'),
     2: (313, 0.347, 'every node (arm C)'),
     3: (round(13682 * 198 / 250), 0.307, 'frontier x a sample of 250, no cut phase')}),
 # the two solved trees, as back-tests: the depth is the node's own, the root included (it has marginals there)
 'tai50b (solved: tai50b-r9)': dict(n=50, run=(116.5, 116.5), machine='c7a core-hours, the run', depths={
     0: (1, 0.70, 'the tree'), 1: (10, 0.70 * (48 / 49) ** 5.2, 'the tree'), 2: (12, 0.70 * (47 / 49) ** 5.2, 'the tree')}),
 'tai60b (solved: tai60b-261117)': dict(n=60, run=(3597, 3597), machine='m7a core-hours, the run', depths={
     0: (1, 1.72 * (59 / 60) ** 5.2, 'the tree'), 1: (5, 1.72 * (58 / 60) ** 5.2, 'the tree'),
     2: (4, 1.72 * (57 / 60) ** 5.2, 'the tree'), 3: (12, 1.72 * (56 / 60) ** 5.2, 'the tree')}),
}
print(f'K = {K} rows and {K} columns, {CAP} iterations a child\n')
print('| instance | firing depth | branching nodes | core-hours a firing | core-hours, this depth | cumulative | share of the run, cumulative |')
print('|---|---|---|---|---|---|---|')
for name, d in I.items():
    cum = 0.0
    for dep in sorted(d['depths']):
        cnt, sit, how = d['depths'][dep]
        m = d['n'] - dep
        one = (2 * K * m - K * K) * CAP * sit / 3600.0
        cum += cnt * one
        lo, hi = d['run']
        share = f'{100 * cum / hi:.1f}%' if lo == hi else f'{100 * cum / hi:.1f}-{100 * cum / lo:.1f}%'
        print(f'| {name} | {dep} | {cnt:,} | {one:.1f} | {cnt * one:,.0f} | {cum:,.0f} | {share} |')
    print(f'| | | | | | | of {d["run"][0]:,.0f}' + (f'-{d["run"][1]:,.0f}' if d['run'][0] != d['run'][1] else '') + f' {d["machine"]} |')
