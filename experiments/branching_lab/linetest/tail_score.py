#!/usr/bin/env python3
"""tail_score.py: the tail rule's arms (runs/scan/tail_*) and the already-bounded tail lines, beside rule 9's line (the
descent) and the best alternative bounded at the node: branching children, shortfall sum, CPU. Writes out/tail_score.txt."""
import glob, os
H = os.path.dirname(os.path.abspath(__file__)); CUT = 608215054.05
def trace(p):
    d = {}
    if os.path.exists(p):
        for l in open(p):
            f = l.split()
            if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
lb = lambda f: float(f[5] if f[5] != '-' else f[4] if f[4] != '-' else f[3])
def st(T, k):
    K = [f for kk, f in T.items() if kk.startswith(k + ',') and kk.count(',') == k.count(',') + 1]
    return sum(f[7] == 'B' for f in K), sum(CUT - lb(f) for f in K if f[7] == 'B'), sum(float(f[11]) for f in K)
DOWN = {}
for p in [H + '/runs/sc47/nodes.txt'] + glob.glob(H + '/runs/d9_*/nodes.txt') + glob.glob(H + '/runs/down/d*/nodes/*.txt'): DOWN.update(trace(p))
ARMS = {}
for d in glob.glob(H + '/runs/scan/*') + glob.glob(H + '/runs/arm17/*'):
    T = trace(d + '/nodes.txt')
    if not T or not any(l.startswith('exit') for l in open(d + '/cmd.txt')): continue
    dep = min(int(f[1]) for f in T.values()); f = [f for f in T.values() if int(f[1]) == dep][0]
    ARMS.setdefault(f[2], {})[f[9]] = st(T, f[2])
rows = [l.split('|') for l in open(H + '/out/tail_rule.txt').read().split('\n')[2:] if l.strip()]
out = ['| node | depth | rule 9: line, B, shortfall sum, CPU s | tail: line, B, shortfall sum, CPU s | best alternative bounded: line, B | tail - rule 9 |', '|---|---|---|---|---|---|']
fail = {'f12', 'p849', 'p857', 'p2231'}
for r in rows:
    name = r[1].strip(); dep = r[2].strip(); tl = r[4].strip().split()[0]
    k = [f for f in trace(H + f'/runs/dump/{name}/nodes.txt').values()]; k = min(k, key=lambda f: int(f[1]))[2]
    r9 = DOWN[k][9]; s9 = st(DOWN, k)
    alts = dict(ARMS.get(k, {}))
    st_t = s9 if tl == r9 else alts.get(tl)
    oth = {l: v for l, v in alts.items() if l not in (tl, r9)}
    ba = min(oth.items(), key=lambda kv: kv[1][0]) if oth else None
    if st_t is None and name not in fail and not os.path.exists(H + f'/runs/scan/tail_{name}'): continue
    tt = f"{tl}: **{st_t[0]}**, {st_t[1]/1e3:,.0f}k, {st_t[2]:,.0f}" if st_t else f"{tl}: running"
    out.append(f"| {name}{' (failure)' if name in fail else ''} | {dep} | {r9}: **{s9[0]}**, {s9[1]/1e3:,.0f}k, {s9[2]:,.0f} | {tt} | "
               f"{(ba[0] + ': ' + str(ba[1][0])) if ba else '-'} | {(st_t[0] - s9[0]) if st_t else '?'} |")
open(H + '/out/tail_score.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
