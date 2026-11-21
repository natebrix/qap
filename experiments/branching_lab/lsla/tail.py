#!/usr/bin/env python3
"""tail.py: the coordinator's "tail" predictor (Nov 21) at the nodes whose lines were bounded here. From each node's
DNN marginals x-bar after its own bound (runs/xdump/NODE/xbar.txt, xdump.sh), per line: entries above 0.001, 0.0001
and 0.01, the sum of squares, and the measured branching children (B, from this README's bounded arms). Rule
"tail": fewest entries above 0.001, ties by the larger sum of squares, then rows before columns.
Writes out/tail.md."""
import os
H = os.path.dirname(os.path.abspath(__file__))
B = {'d7': {'c59': 3, 'r11': 4, 'c47': 5, 'c10': 5, 'c54': 9},
     'd8': {'c59': 3, 'r3': 4, 'c2': 7},
     'd9': {'c59': 3, 'r8': 4, 'c53': 4, 'r3': 4, 'c30': 12},
     'r50': {'r31': 10, 'c26': 28}}
NAME = {'d7': 'tai60b 10:39 (depth 7)', 'd8': 'tai60b 14:54 (depth 8)', 'd9': 'tai60b 24:2 (depth 9)',
        'r50': 'tai50b root'}
if os.path.exists(f'{H}/runs/root50_c25/cmd.txt') and 'exit' in open(f'{H}/runs/root50_c25/cmd.txt').read():
    B['r50']['c25'] = sum(1 for l in open(f'{H}/runs/root50_c25/nodes.txt') if l.split()[:2] == ['node', '1'] and l.split()[7] == 'B')
def load(p):
    rows = cols = None; X = []; hdr = ''
    for l in open(p):
        if l.startswith('xbar'): hdr = l.strip()
        elif l.startswith('rows'): rows = [int(x) for x in l.split()[1:]]
        elif l.startswith('cols'): cols = [int(x) for x in l.split()[1:]]
        elif l.startswith('x '): X.append([float(x) for x in l.split()[1:]])
    return hdr, rows, cols, X
def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]: j += 1
        for k in range(i, j + 1): r[s[k]] = (i + j) / 2 + 1
        i = j + 1
    return r
def spearman(a, b):
    ra, rb = rank(a), rank(b); n = len(a); ma = sum(ra) / n; mb = sum(rb) / n
    c = sum((x - ma) * (y - mb) for x, y in zip(ra, rb)); va = sum((x - ma) ** 2 for x in ra); vb = sum((y - mb) ** 2 for y in rb)
    return c / (va * vb) ** 0.5 if va and vb else float('nan')
out = ['# The "tail" predictor at the bounded nodes (tail.py)', '']
allx, ally = [], []
for N in ['d7', 'd8', 'd9', 'r50']:
    p = f'{H}/runs/xdump/{N}/xbar.txt'
    if not os.path.exists(p): out += [f'## {NAME[N]}: no dump yet', '']; continue
    hdr, rows, cols, X = load(p); m = len(X)
    L = {}
    for a, r in enumerate(rows):
        v = X[a]; L[f'r{r}'] = v
    for b, c in enumerate(cols):
        L[f'c{c}'] = [X[a][b] for a in range(m)]
    S = {k: (sum(x > 1e-3 for x in v), sum(x > 1e-4 for x in v), sum(x > 1e-2 for x in v), sum(x * x for x in v)) for k, v in L.items()}
    order = sorted(S, key=lambda k: (S[k][0], -S[k][3], 0 if k[0] == 'r' else 1))
    pos = {k: i + 1 for i, k in enumerate(order)}
    mn = S[order[0]][0]; ties = sum(1 for k in S if S[k][0] == mn)
    out += [f'## {NAME[N]}', '', f'`{hdr}`', '',
            f'"tail" picks **{order[0]}** ({mn} entries above 0.001; {ties} line(s) tie at that minimum; '
            f'bounded: {"yes, B = %d" % B[N][order[0]] if order[0] in B[N] else "no"}). Next: ' +
            ', '.join(f'{k} ({S[k][0]}, {S[k][3]:.3f})' for k in order[1:6]) + '.', '',
            '| line | B measured | > 0.001 | > 0.0001 | > 0.01 | sum x^2 | tail rank |', '|---|---|---|---|---|---|---|']
    xs, ys = [], []
    for k in sorted(B[N], key=lambda k: B[N][k]):
        s = S[k]; out.append(f'| {k} | {B[N][k]} | {s[0]} | {s[1]} | {s[2]} | {s[3]:.4f} | {pos[k]} of {len(S)} |')
        xs.append(s[0]); ys.append(B[N][k])
    allx += xs; ally += ys
    out += ['', f'Spearman (entries above 0.001 against B), this node: {spearman(xs, ys):+.2f} over {len(xs)} lines.', '']
if len(allx) > 2: out += [f'Pooled over the nodes dumped: Spearman {spearman(allx, ally):+.2f} over {len(allx)} lines (pooling mixes nodes).', '']
os.makedirs(f'{H}/out', exist_ok=True); open(f'{H}/out/tail.md', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
