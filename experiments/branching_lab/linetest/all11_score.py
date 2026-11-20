#!/usr/bin/env python3
"""all11_score.py: rule 9 (the descent, runs/down) against rule 2 (runs/scan/s16_*, s17_*; the three depth-17 parents
of the line test from runs/arm17/*_rule2) at every branched depth-16 and depth-17 node of the c47 descent. Per node:
branching children, shortfall sum, CPU, under each; totals per depth and split by on / off the best known solution's
path. Writes out/all11.txt."""
import glob, os, collections
H = os.path.dirname(os.path.abspath(__file__)); CUT = 608215054.05
t = open('/home/user/qap/QAPSOLN/tai60b.sln').read().split(); p0 = [int(x) - 1 for x in t[2:62]]
def trace(paths):
    d = {}
    for p in paths:
        if os.path.exists(p):
            for l in open(p):
                f = l.split()
                if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
lb = lambda f: float(f[5] if f[5] != '-' else f[4] if f[4] != '-' else f[3])
DOWN = trace(glob.glob(H + '/runs/down/d*/nodes/*.txt'))
R2 = {}
for d in glob.glob(H + '/runs/scan/s1*') + glob.glob(H + '/runs/arm17/*_rule2'):
    T = trace([d + '/nodes.txt'])
    if not T or not any(l.startswith('exit') for l in open(d + '/cmd.txt')): continue
    dep = min(int(f[1]) for f in T.values()); k = [f for f in T.values() if int(f[1]) == dep][0][2]; R2[k] = T
def kids(T, k): return [f for kk, f in T.items() if kk.startswith(k + ',') and kk.count(',') == k.count(',') + 1]
def stat(L): return sum(f[7] == 'B' for f in L), sum(CUT - lb(f) for f in L if f[7] == 'B'), sum(float(f[11]) for f in L)
out = []
for D in (17, 16):
    nodes = sorted(k for k, f in DOWN.items() if int(f[1]) == D and f[7] == 'B')
    out += [f'depth {D}: {len(nodes)} branched nodes in the descent (the whole level)', '',
            '| node | path | rule 9: line, B, shortfall sum, CPU s | rule 2: line, B, shortfall sum, CPU s |', '|---|---|---|---|']
    tot = collections.defaultdict(lambda: [0, 0.0, 0.0, 0, 0.0, 0.0, 0])
    for k in nodes:
        on = all(p0[int(a)] == int(b) for a, b in (x.split(':') for x in k.split(',')))
        b9 = stat(kids(DOWN, k)); l9 = DOWN[k][9]
        if k in R2:
            T = R2[k]; b2 = stat(kids(T, k)); l2 = T[k][9]; c2 = f'{l2}: **{b2[0]}**, {b2[1]/1e3:,.0f}k, {b2[2]:,.0f}'
        else: b2 = None; c2 = 'not done'
        out.append(f"| {k.split(',', 7)[-1]} | {'yes' if on else 'no'} | {l9}: **{b9[0]}**, {b9[1]/1e3:,.0f}k, {b9[2]:,.0f} | {c2} |")
        for g in ('all', 'path' if on else 'off'):
            x = tot[g]
            if b2:
                x[0] += b9[0]; x[1] += b9[1]; x[2] += b9[2]; x[3] += b2[0]; x[4] += b2[1]; x[5] += b2[2]; x[6] += 1
    for g in ('all', 'path', 'off'):
        x = tot[g]
        out.append(f'| **{g}** ({x[6]} nodes with both) | | rule 9: **{x[0]}**, {x[1]/1e3:,.0f}k, {x[2]:,.0f} | rule 2: **{x[3]}**, {x[4]/1e3:,.0f}k, {x[5]:,.0f} |')
    out.append('')
open(H + '/out/all11.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
