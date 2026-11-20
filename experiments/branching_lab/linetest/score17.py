#!/usr/bin/env python3
"""score17.py: the depth-17 line test. For each of the three depth-17 parents under c47 (8:49, 8:57, the path node
22:31): the taken line's children (production's rule, from runs/down/d18/nodes) and the forced arm's
(runs/arm17/NAME/nodes.txt): D / C / B (other), the sum and the largest of the branching children's shortfalls, CPU;
with the signal's predicted survivors (out/signal_TAG_k300.tsv, excess <= the parent's own slack). Writes out/score17.txt."""
import glob, os, collections
H = os.path.dirname(os.path.abspath(__file__)); CUT = 608215054.05
P = '44:45,43:58,1:33,6:24,7:19,0:32,10:39,'
A = P + '46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,'
PARENTS = [('8:49', A + '8:49', 'p849', 'p849_r51'), ('8:57', A + '8:57', 'p857', 'p857_r51'),
           ('22:31', P + '47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6,22:31', 'p2231', 'p2231_r13')]
def trace(paths):
    d = {}
    for p in paths:
        if os.path.exists(p):
            for ln in open(p):
                f = ln.split()
                if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
def lb(f):
    for i in (5, 4, 3):
        if f[i] != '-': return float(f[i])
down = trace(glob.glob(H + '/runs/down/d18/nodes/*.txt') + glob.glob(H + '/runs/down/d17/nodes/*.txt'))
def kids(T, par):
    return [f for k, f in T.items() if k.startswith(par + ',') and k.count(',') == par.count(',') + 1]
def row(L):
    c = collections.Counter(f[7] for f in L); sf = [CUT - lb(f) for f in L if f[7] == 'B']
    return f"{len(L)} children: D {c['D']} C {c['C']} B {c['B']} other {len(L) - c['D'] - c['C'] - c['B']}; B shortfall sum {sum(sf)/1e3:,.0f}k max {max(sf, default=0)/1e3:,.0f}k; CPU {sum(float(f[11]) for f in L):,.0f} s", c['B']
def predicted(tag, line):
    p = f'{H}/out/signal_{tag}_k300.tsv'
    if not os.path.exists(p): return '?'
    return sum(1 for l in open(p).read().split('\n')[1:] if l and l.split('\t')[0] == line and float(l.split('\t')[3]) <= 1)
tot_t = tot_a = 0
for name, par, tag, arm in PARENTS:
    pf = down[par]; taken = pf[9]
    A_ = trace([f'{H}/runs/arm17/{arm}/nodes.txt']); af = A_.get(par)
    rt, bt = row(kids(down, par)); tot_t += bt
    print(f'== parent {name}: lb {lb(pf):,.0f} (slack {CUT - lb(pf):,.0f})')
    print(f'  taken {taken} (predicted {predicted(tag, taken)}): {rt}')
    if af:
        ra, ba = row(kids(A_, par)); tot_a += ba
        print(f'  arm {af[9]} (predicted {predicted(tag, af[9])}): parent bound {lb(af):,.0f} (run {lb(pf):,.0f}); {ra}')
    else: print('  arm: not bounded yet')
print(f'total branching children: taken {tot_t}, signal arms {tot_a}')
