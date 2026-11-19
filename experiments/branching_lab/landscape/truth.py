#!/usr/bin/env python3
"""truth.py SET [KICKS]: every child of every candidate line at the truth nodes gets its excess (best completion
found - best known), for searches of increasing effort: one descent from each start (ex0), then an ILS of 30, 100,
KICKS (default 300) kicks from the best of those (ex30, ex100, exK). Starts: the best known permutation and its
images under the instance's symmetries (tho: the 4 automorphisms of the grid; tai60b: p0 and the mirror), each
adjusted to the child's fixes.
SET: tho30 (3 depth-1 nodes, 8 lines each, ../tho30_solve/line_truth/results), tho40 (3 depth-3 nodes, 8 lines,
../tho40_line_truth/results; instance tho40t), tai60b (43:58 under 44:45, the line study's 4 lines with per-child
outcomes from its arms), tai60ball (the same node, every free line).
Writes truth/SET.tsv: node, line, child fixes, the truth (subtree nodes and whether it branched, or the arm's
outcome), ex0, ex30, ex100, exK."""
import sys, os, json, glob, time, numpy as np, lsq, autom

SET = sys.argv[1]; KICKS = int(sys.argv[2]) if len(sys.argv) > 2 else 300
rng = np.random.default_rng(11)
EXP = f'{lsq.WT}/experiments'

def kids(node_fix, line, n):
    used_r = {r for r, _ in node_fix}; used_c = {c for _, c in node_fix}
    side, k = line[0], int(line.lstrip('rowcl'))
    if side == 'r': return [(k, v) for v in range(n) if v not in used_c]
    return [(r, k) for r in range(n) if r not in used_r]

def excess(I, starts, fix, marks):
    best = min((I.descend(I.adjust(s, fix), fix) for s in starts), key=lambda x: x[0])
    c, p = best; free = np.array([i for i in range(I.n) if i not in fix]); out = {0: c}
    for it in range(1, max(marks) + 1):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix)
        if cq < c: c, p = cq, q
        if it in marks: out[it] = c
    return [out[m] - I.best for m in [0] + marks]

marks = sorted({30, 100, KICKS})
rows = []
if SET in ('tho30', 'tho40'):
    if SET == 'tho30':
        t = open(lsq.H + '/truth/tho30.perm').read().split()   # QAPSOLN/tho30.sln has no permutation: one found by ts.c
        I = lsq.Inst('tho30', perm=[int(x) for x in t[1:]], best=int(t[0])); files = sorted(glob.glob(f'{EXP}/tho30_solve/line_truth/results/*.json'))
    else:
        t = open(f'{lsq.Q}/QAPSOLN/tho40.sln').read().split(); p = np.array([int(x) - 1 for x in t[2:42]])
        I = lsq.Inst('tho40t', dat=f'{lsq.Q}/qaplib/tho40t.dat', perm=np.argsort(p), best=int(t[1]))
        files = sorted(glob.glob(f'{EXP}/tho40_line_truth/results/*.json'))
    SA = autom.automorphisms(I.A); SB = autom.automorphisms(I.B)
    starts = {tuple(sb[I.p0[sa]]) for sa in SA for sb in SB}       # p(i) -> sb(p(sa(i)))
    starts = [np.array(s, dtype=np.int32) for s in starts]
    assert all(abs(I.cost(s) - I.best) < 0.5 for s in starts)
    for f in files:
        d = json.load(open(f)); node = d['node']; nf = lsq.parse_fixes(node)
        for line, L in d['lines'].items():
            for ch, rec in L['children'].items():
                r, c = map(int, ch.split(':')); fix = dict(nf); fix[r] = c
                tr = rec['truth']
                rows.append([node, line, ch, tr['nodes'], int(tr['nodes'] > 1), rec['score'].get('lb')] + excess(I, starts, fix, marks))
            print(SET, node, line, flush=True)
else:
    I = lsq.Inst('tai60b')
    mir = lsq.mirror(I)
    starts = [I.p0, mir]; node = '44:45,43:58'; nf = lsq.parse_fixes(node)
    if SET == 'tai60b':
        for arm in ['r1', 'r45', 'c26', 'r8']:
            T = {t['fixes']: t for t in lsq.read_trace(f'{EXP}/tai60b_lines/runs/{arm}/nodes.txt') if t['depth'] == 3}
            for (r, c) in kids(nf, arm, I.n):
                fx = f'{node},{r}:{c}'; t = T[fx]; fix = dict(nf); fix[r] = c
                rows.append([node, arm, f'{r}:{c}', t['out'], int(t['out'] == 'B'), t['dnn']] + excess(I, starts, fix, marks))
            print(SET, arm, flush=True)
    else:
        done = {}
        for side in 'rc':
            for k in range(I.n):
                if (side == 'r' and k in (44, 43)) or (side == 'c' and k in (45, 58)): continue
                for (r, c) in kids(nf, f'{side}{k}', I.n):
                    if (r, c) not in done:
                        fix = dict(nf); fix[r] = c; done[(r, c)] = excess(I, starts, fix, marks)
                    rows.append([node, f'{side}{k}', f'{r}:{c}', '-', '-', '-'] + done[(r, c)])
            print(SET, side, len(done), flush=True)
os.makedirs(lsq.H + '/truth', exist_ok=True)
with open(f'{lsq.H}/truth/{SET}.tsv', 'w') as w:
    w.write('\t'.join(['node', 'line', 'child', 'truth', 'branched', 'lb'] + [f'ex{m}' for m in [0] + marks]) + '\n')
    for r in rows: w.write('\t'.join(str(x) for x in r) + '\n')
print(SET, 'done', len(rows))
