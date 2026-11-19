#!/usr/bin/env python3
"""d9score.py: one level deeper. Per depth-8 survivor of the strong c47 arm (runs/d9_*/nodes.txt), and per
depth-8 branching child of production's c54 (tai60b-261117): the line taken at depth 8, the depth-9 children's
D / C / B counts, the branching ones' shortfalls (cutoff - final bound), CPU. Writes out/d9score.txt."""
import glob, gzip, os, collections
H = os.path.dirname(os.path.abspath(__file__))
NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39,'
fv = lambda s: None if s == '-' else float(s)
def lb(f): return fv(f[5]) if f[5] != '-' else (fv(f[4]) if f[4] != '-' else fv(f[3]))
def report(name, T, lines):
    by = collections.defaultdict(list)
    for f in T.values():
        if int(f[1]) == 9 and f[2].startswith(NODE): by[f[2][len(NODE):].split(',')[0]].append(f)
    tot = collections.Counter(); tcpu = 0
    print(f'== {name}')
    print('depth-8 node  line  depth-9 bounded  D  C  B  CPU s   branching shortfalls')
    for k in sorted(by, key=lambda k: -len(by[k])):
        L = by[k]; c = collections.Counter(f[7] for f in L); cpu = sum(float(f[11]) for f in L)
        sf = sorted((float(f[6]) - lb(f) for f in L if f[7] == 'B'), reverse=True)
        tot.update(c); tcpu += cpu
        print(f'{k:8s} {lines.get(k, "?"):>5s}  {len(L):6d}  {c["D"]:3d} {c["C"]:3d} {c["B"]:3d}  {cpu:8,.0f}   ' + ' '.join(f'{s/1e3:,.0f}k' for s in sf))
    print(f'total: {sum(tot.values())} bounded, D {tot["D"]} C {tot["C"]} B {tot["B"]}, CPU {tcpu:,.0f} s')
def trace(paths):
    d = {}
    for p in paths:
        op = gzip.open if p.endswith('.gz') else open
        for ln in op(p, 'rt'):
            f = ln.split()
            if len(f) >= 12 and f[0] == 'node' and f[2].startswith(NODE): d[f[2]] = f
    return d
sc = trace([H + '/runs/sc47/nodes.txt'])
lines = {f[2][len(NODE):]: f[9] for f in sc.values() if int(f[1]) == 8 and f[7] == 'B'}
report('strong c47, depth 9 (runs/d9_*)', trace(sorted(glob.glob(H + '/runs/d9_*/nodes.txt'))), lines)
P = trace(['/home/user/wt-cognac/infra/phase1/runs/tai60b-261117/nodes.txt.gz'])
plines = {f[2][len(NODE):]: f[9] for f in P.values() if int(f[1]) == 8 and f[7] == 'B'}
report('production c54, depth 9 (tai60b-261117)', P, plines)
