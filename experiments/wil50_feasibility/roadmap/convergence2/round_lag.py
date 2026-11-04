#!/usr/bin/env python3
"""The cut phase's per-round lag (Lagrangian primal - certified lb at the end of
each Polyak round) and the quiet round's gain, from the budget study's logs:
arm B restarts every round from Y = e0e0' (as the search does today), Y12 and
Y12s restore Y as well as Z (dnncuts -warmY).
usage: round_lag.py <repo root> [extra logs...]
"""
import os, re, sys

root = sys.argv[1]
bs = os.path.join(root, 'experiments/same_node_continuation/budget_study/results/runs')


def parse(p):
    out = []
    for l in open(p):
        m = re.match(r'round\s+(\d+)\s+lb (\S+) .*?primal (\S+)\s+status', l)
        q = re.match(r'quiet round\s+lb (\S+) .*?primal (\S+)\s+status', l)
        if m:
            out.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
        elif q:
            out.append(('q', float(q.group(1)), float(q.group(2))))
    return out


def report(name, r):
    rounds = [x for x in r[1:] if x[0] != 'q']
    lags = [pr - lb for k, lb, pr in rounds]
    if not lags:
        return
    med = sorted(lags)[len(lags) // 2]
    qq = [x for x in r if x[0] == 'q']
    qg = f"{qq[0][1] - max(lb for k, lb, pr in rounds):6.1f}" if qq else '     -'
    print(f"{name:14} rounds {len(lags):2d}  median lag {med:6.1f}  quiet gain over best round {qg}   "
          f"lags {' '.join(f'{x:.0f}' for x in lags)}")


for n in ['n14', 'n6', 'n9', 'n13', 'd2', 'n10']:
    for arm in ['B', 'Y12', 'Y12s']:
        report(f'{n} {arm}', parse(os.path.join(bs, f'{n}_{arm}.log')))
for p in sys.argv[2:]:
    report(os.path.basename(p), parse(p))
