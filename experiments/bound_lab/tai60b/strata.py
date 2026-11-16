#!/usr/bin/env python3
"""strata.py PARENTRUN DEPTH: split the children (depth DEPTH) of the on-path node of runs/PARENTRUN (its
out.state; the root: rootM, root.state) into the on-path child, a 'near' stratum and a 'far' stratum, for the
dives. Every child gets the excess of the best completion a 2-swap descent finds for it (model.py's rule); near =
excess below NEARF (default 2) times the first-solve shortfall of the on-path node at DEPTH-1... no: at DEPTH
if the chain has it, else the deepest measured one. Writes runs/strata/offD_near.state, offD_far.state and
runs/strata/offD.txt (one line per child: stratum, fixes, excess), and prints the counts."""
import sys, os, re, glob, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/../regions')
import gen_regions as g
Q = H.rsplit('/', 3)[0]
n, A, B = g.read_qaplib(f'{Q}/qaplib/tai60b.dat'); Bs = (B + B.T) / 2.0
t = open(f'{Q}/QAPSOLN/tai60b.sln').read().split(); best, p0 = int(t[1]), np.array(list(map(int, t[2:2 + n]))) - 1
def cost(p): return float((A * Bs[np.ix_(p, p)]).sum())
def deltas(p):
    Bp = Bs[np.ix_(p, p)]; M = A @ Bp.T; s = (A * Bp).sum(1)
    return 2 * (M + M.T - s[:, None] - s[None, :] + 2 * A * Bp)
def descent(p, c, free):
    while True:
        D = deltas(p)[np.ix_(free, free)]
        k = int(np.argmin(D)); a, b = divmod(k, len(free))
        if D[a, b] > -0.5: return p, c
        i, j = free[a], free[b]; p[i], p[j] = p[j], p[i]; c += D[a, b]
def excess(fix):
    """fix: {row: value}; the best completion found from p0 with those entries forced"""
    p = p0.copy()
    for r, v in fix.items():
        if p[r] != v:
            j = int(np.where(p == v)[0][0]); p[j], p[r] = p[r], p[j]
    free = [i for i in range(n) if i not in fix]
    return descent(p, cost(p), free)[1] - best
def first_short():
    """on-path first-solve shortfall by depth, from the chain's traces"""
    out = {}
    for f in [x for x in glob.glob(H + '/runs/d*path/dnn.txt') if 'under' not in x]:
        m = re.match(r'pid \d+ depth (\d+) fixes.*?\| first lb (\S+)', open(f).read())
        if m: out[int(m[1])] = best - float(m[2])
    return out
if __name__ == '__main__':
    run, D = sys.argv[1], int(sys.argv[2]); nearf = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    src = f'{H}/runs/{run}/' + ('root.state' if D == 1 else 'out.state')
    fs = first_short(); f = fs.get(D, fs[max(fs)])
    head, nodes, tail = [], [], []
    for ln in open(src):
        if ln.startswith('nodes '): continue
        if ln.startswith('node '): nodes.append(ln)
        elif nodes: tail.append(ln)
        else: head.append(ln)
    os.makedirs(H + '/runs/strata', exist_ok=True)
    S = {'near': [], 'far': [], 'path': []}; rows = []
    for ln in nodes:
        tk = ln.split(); assert int(tk[1]) == D
        fix = {i: int(tk[5 + i]) for i in range(n) if int(tk[5 + i]) >= 0}
        if all(p0[r] == v for r, v in fix.items()): k, ex = 'path', 0.0
        else:
            ex = excess(fix); k = 'near' if ex < nearf * f else 'far'
        S[k].append(ln); rows.append((ex, k, ','.join(f'{r}:{v}' for r, v in fix.items() if p0[r] != v)))
    for k in ('near', 'far'):
        with open(f'{H}/runs/strata/off{D}_{k}.state', 'w') as w:
            w.writelines(head); w.write(f'nodes {len(S[k])}\n'); w.writelines(S[k]); w.writelines(tail)
    with open(f'{H}/runs/strata/off{D}.txt', 'w') as w:
        w.write(f'# depth {D} children of runs/{run}: near = excess < {nearf} x {f:.0f} (the on-path first-solve shortfall)\n')
        for ex, k, fx in sorted(rows): w.write(f'{k} {fx or "(on-path)"} {ex:.0f}\n')
    print(f'depth {D}: {len(nodes)} children, near {len(S["near"])} (excess < {nearf * f / 1e6:.2f}M), far {len(S["far"])}, on-path {len(S["path"])}')
