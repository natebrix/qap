#!/usr/bin/env python3
"""clusters_xbar.py: at each dumped node (runs/dump/*), the marginal mass x-bar of rule 9's line, the best alternative
line and the "support" line, by facility cluster (tai60b's six clusters from the census,
experiments/branching_lab/failures/census/out/clusters.json on claude/cool-babbage-ilnogj: 0-13, 14-25, 26-42+45-47,
43-44, 48-51, 52-59). Only columns have facility children; a row's children are locations, which have no cluster
structure (census), so rows are reported by their own facility's cluster only. Per column: mass per cluster, the
number of free members of each cluster and how many carry x > 1/m and x > 0.001; and for the realised children, which
clusters branched. Writes out/clusters_xbar.txt."""
import glob, os, sys
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
CL = [[0,1,2,3,4,5,6,7,8,9,10,11,12,13],[14,15,16,17,18,19,20,21,22,23,24,25],[26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,45,46,47],[43,44],[48,49,50,51],[52,53,54,55,56,57,58,59]]
cof = {f: k for k, c in enumerate(CL) for f in c}
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
import glob as G
DOWN = {}
for p in ['runs/sc47/nodes.txt'] + G.glob('runs/d9_*/nodes.txt') + G.glob('runs/down/d*/nodes/*.txt'): DOWN.update(trace(p))
ARMS = {}
for d in G.glob('runs/scan/*') + G.glob('runs/arm17/*') + G.glob('runs/alt/*'):
    T = trace(d + '/nodes.txt')
    if T:
        dep = min(int(f[1]) for f in T.values()); f = [f for f in T.values() if int(f[1]) == dep][0]; ARMS.setdefault(f[2], {})[f[9]] = T
out = []
for d in sorted(G.glob('runs/dump/*')):
    if not os.path.exists(d + '/done.txt'): continue
    X, rows, cols = parse(d + '/err.txt')
    if not X: continue
    m = len(X); node = [f for f in trace(d + '/nodes.txt').values()]; node = min(node, key=lambda f: int(f[1])); k = node[2]
    lines = {node[9]: DOWN}
    for ln, T in ARMS.get(k, {}).items(): lines[ln] = T
    out.append(f"== {os.path.basename(d)} (depth {node[1]}, 1/m = {1/m:.4f})")
    for ln, T in lines.items():
        bkids = [f[2].rsplit(',', 1)[1] for f in T.values() if f[2].startswith(k + ',') and f[2].count(',') == k.count(',') + 1 and f[7] == 'B']
        if ln[0] == 'c':
            j = cols.index(int(ln[1:])); ent = {rows[i]: X[i][j] for i in range(m)}
            mass = [sum(v for f, v in ent.items() if cof[f] == c) for c in range(6)]
            free = [sum(1 for f in ent if cof[f] == c) for c in range(6)]
            big = [sum(1 for f, v in ent.items() if cof[f] == c and v > 1 / m) for c in range(6)]
            nz = [sum(1 for f, v in ent.items() if cof[f] == c and v > 1e-3) for c in range(6)]
            bc = [0] * 6
            for b in bkids: bc[cof[int(b.split(':')[0])]] += 1
            out.append(f"  {ln}: {len(bkids)} branching children, by cluster {bc}")
            out.append(f"      mass by cluster {[round(x, 3) for x in mass]}; free members {free}; entries > 1/m {big}; entries > 0.001 {nz}")
        else:
            i = rows.index(int(ln[1:])); e = X[i]
            out.append(f"  {ln} (row; facility {ln[1:]} in cluster {cof[int(ln[1:])]}): {len(bkids)} branching children; entries > 1/m {sum(v > 1/m for v in e)}, > 0.001 {sum(v > 1e-3 for v in e)}")
open(H + '/out/clusters_xbar.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
