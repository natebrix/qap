#!/usr/bin/env python3
"""pathcheck.py: the best known solution's path (p0, QAPSOLN) through the c47 runs: at every depth, the node whose
fixes all agree with p0 (sc47 arm, runs/d9_*, runs/down/d*), its outcome and bound. Each must branch (B) or be a
leaf (E) with bound <= 608,215,054 (the cutoff; the incumbent was best known + 1). Writes out/pathcheck.txt."""
import glob, os
H = os.path.dirname(os.path.abspath(__file__))
t = open('/home/user/qap/QAPSOLN/tai60b.sln').read().split(); p0 = [int(x) - 1 for x in t[2:62]]
files = [H + '/runs/sc47/nodes.txt'] + sorted(glob.glob(H + '/runs/d9_*/nodes.txt')) + sorted(glob.glob(H + '/runs/down/d*/nodes/*.txt'))
on = {}
for p in files:
    for ln in open(p):
        f = ln.split()
        if len(f) < 12 or f[0] != 'node': continue
        if all(p0[int(a)] == int(b) for a, b in (x.split(':') for x in f[2].split(','))): on[f[2]] = f
ok = True
for k, f in sorted(on.items(), key=lambda kv: int(kv[1][1])):
    b = float(f[5] if f[5] != '-' else f[4] if f[4] != '-' else f[3])
    good = f[7] in 'BE' and b <= 608215054.05
    ok &= good
    print(f'depth {f[1]:>2} {f[2].split(",")[-1]:>6} {f[7]} bound {b:,.0f} line {f[9]} {"ok" if good else "FAIL"}')
print('path check:', 'ok' if ok else 'FAILED', f'({len(on)} path nodes, depths {min(int(f[1]) for f in on.values())}-{max(int(f[1]) for f in on.values())})')
