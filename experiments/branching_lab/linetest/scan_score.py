#!/usr/bin/env python3
"""scan_score.py: the rule 2 / rule 9 scan (README). Per scanned node: depth, on the best known solution's path, its QPB
bound as a fraction below the cutoff, rule 9's line (from the c47 descent, runs/down) and rule 2's (runs/scan/NAME),
each with its children's D / C / B (other), branching shortfall sum and CPU; rule 2's count (QAP_MAXF_DEBUG: how many
children QPB + U put at the cutoff, for its chosen line, and the best of the other side). Writes out/scan.txt."""
import glob, os, collections, re
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
DOWN = trace(glob.glob(H + '/runs/down/d*/nodes/*.txt') + glob.glob(H + '/runs/d9_*/nodes.txt'))
def kids(T, k): return [f for kk, f in T.items() if kk.startswith(k + ',') and kk.count(',') == k.count(',') + 1]
def row(L):
    c = collections.Counter(f[7] for f in L); sf = sum(CUT - lb(f) for f in L if f[7] == 'B')
    o = len(L) - c['D'] - c['C'] - c['B']
    return f"{c['D']} / {c['C']} / **{c['B']}**" + (f" ({o} Q)" if o else ''), sf, sum(float(f[11]) for f in L), c['B']
out = ['| node (depth) | path | QPB below cutoff | rule 9: line, D / C / B, shortfall sum, CPU s | rule 2: line, D / C / B, shortfall sum, CPU s | rule 2 count (chosen side / other side) |', '|---|---|---|---|---|---|']
for d in sorted(glob.glob(H + '/runs/scan/*')):
    T = trace([d + '/nodes.txt'])
    if not T: continue
    dep = min(int(f[1]) for f in T.values()); node = [f for f in T.values() if int(f[1]) == dep][0]; k = node[2]
    r9 = DOWN[k]; on = all(p0[int(a)] == int(b) for a, b in (x.split(':') for x in k.split(',')))
    a9 = row(kids(DOWN, k)); a2 = row(kids(T, k)); done = any(l.startswith('exit') for l in open(d + '/cmd.txt'))
    cnt = '?'
    for l in open(d + '/err.txt'):
        m = re.match(r'maxfathoms: depth (\d+) best_row (-?\d+) \(count (\d+) usum (\S+)\) best_col (-?\d+) \(count (\d+) usum (\S+)\)', l)
        if m and int(m.group(1)) == dep:
            rc, cc = int(m.group(3)), int(m.group(6)); ch = node[9][0]
            cnt = f'{rc} / {cc}' if ch == 'r' else f'{cc} / {rc}'; break
    q = (CUT - float(node[3])) / CUT
    out.append(f"| {os.path.basename(d)}: {k.split(',',7)[-1]} ({dep}) | {'yes' if on else 'no'} | {q:.3f} | {r9[9]}: {a9[0]}, {a9[1]/1e3:,.0f}k, {a9[2]:,.0f} | "
               f"{node[9]}: {a2[0]}, {a2[1]/1e3:,.0f}k, {a2[2]:,.0f}{'' if done else ' (running)'} | {cnt} |")
open(H + '/out/scan.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
