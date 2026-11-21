#!/usr/bin/env python3
"""tail_rule.py: the "tail" rule (the coordinator's, fixed Nov 21 00:42): choose the free line with the fewest x-bar
entries above TH (0.001); ties by the larger sum of squares (rule 9's score), then rows before columns. At each dumped
node: its line at TH = 0.0001, 0.001, 0.01, how many lines tie at the minimum, whether the line is already bounded
(rule 9's in the descent, or an arm in runs/scan, runs/arm17, runs/alt) and with how many branching children.
Writes out/tail_rule.txt and runs/tail_jobs.txt (the lines to bound)."""
import glob, os, sys
H = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0]]
def parse(p):
    X = []; rows = cols = None
    for l in open(p):
        if l.startswith('ls-rows:'): rows = [int(x) for x in l.split()[1:]]
        elif l.startswith('ls-cols:'): cols = [int(x) for x in l.split()[1:]]
        elif l.startswith('ls-X '): X.append([float(x) for x in l.split(':', 1)[1].split()])
    return X, rows, cols
def trace(p):
    d = {}
    if os.path.exists(p):
        for l in open(p):
            f = l.split()
            if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
DOWN = {}
for p in [H + '/runs/sc47/nodes.txt'] + glob.glob(H + '/runs/d9_*/nodes.txt') + glob.glob(H + '/runs/down/d*/nodes/*.txt'): DOWN.update(trace(p))
BOUND = {}
def nb(T, k): return sum(1 for kk, f in T.items() if kk.startswith(k + ',') and kk.count(',') == k.count(',') + 1 and f[7] == 'B')
for d in glob.glob(H + '/runs/scan/*') + glob.glob(H + '/runs/arm17/*') + glob.glob(H + '/runs/alt/*'):
    T = trace(d + '/nodes.txt')
    if not T or not any(l.startswith('exit') for l in open(d + '/cmd.txt')): continue
    dep = min(int(f[1]) for f in T.values()); f = [f for f in T.values() if int(f[1]) == dep][0]
    BOUND.setdefault(f[2], {})[f[9]] = nb(T, f[2])
jobs = {l.split()[0]: l.split()[1:] for l in open(H + '/runs/dump_jobs.txt')}
jobs['f12'] = [H + '/runs/rebuild/b335/d11.out.state', '12', '11:29']
def tail(X, rows, cols, th):
    m = len(X); L = []
    for i, r in enumerate(rows): e = X[i]; L.append((sum(v > th for v in e), -sum(v * v for v in e), 0, f'r{r}'))
    for j, c in enumerate(cols): e = [X[i][j] for i in range(m)]; L.append((sum(v > th for v in e), -sum(v * v for v in e), 1, f'c{c}'))
    L.sort(); return L[0][3], L[0][0], sum(1 for x in L if x[0] == L[0][0])
out = ['| node | depth | rule 9: line (B) | tail @0.001: line (entries; lines tied at the min) | @0.0001 | @0.01 | tail line bounded? (B) |', '|---|---|---|---|---|---|---|']
need = []
for d in sorted(glob.glob(H + '/runs/dump/*'), key=lambda d: (int([l for l in open(d + '/err.txt') if l.startswith('ls-node')][0].split()[2]), os.path.basename(d))):
    X, rows, cols = parse(d + '/err.txt')
    if not X: continue
    name = os.path.basename(d); T = trace(d + '/nodes.txt'); f = min(T.values(), key=lambda f: int(f[1])); k = f[2]; dep = int(f[1])
    r9 = DOWN[k][9]; r9b = nb(DOWN, k)
    t3 = tail(X, rows, cols, 1e-3); t4 = tail(X, rows, cols, 1e-4); t2 = tail(X, rows, cols, 1e-2)
    known = dict(BOUND.get(k, {})); known[r9] = r9b
    b = known.get(t3[0])
    out.append(f"| {name} | {dep} | {r9} ({r9b}) | {t3[0]} ({t3[1]}; {t3[2]}) | {t4[0]} ({t4[1]}; {t4[2]}) | {t2[0]} ({t2[1]}; {t2[2]}) | {'yes (' + str(b) + ')' if b is not None else 'no'} |")
    if b is None: need.append((name, dep, t3[0]))
open(H + '/out/tail_rule.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
with open(H + '/runs/tail_jobs.txt', 'w') as w:
    for name, dep, ln in need: w.write(f"{name} {dep} {ln} {' '.join(jobs[name])}\n")
print('not bounded:', need)
