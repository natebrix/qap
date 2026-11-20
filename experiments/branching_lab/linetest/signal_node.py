#!/usr/bin/env python3
"""signal_node.py TAG FIXES LB KICKS [LINES]: the landscape signal at one node (../landscape README, section 4): every
child of every free line (or of LINES, comma-separated) gets its excess, the best of 2-exchange descents from p0, the
mirror and the ILS-start solution adjusted to the child's fixes, then an ILS of KICKS kicks (6 random exchanges + a
descent; seed 0). s = cutoff - LB (the node's own shortfall after its cut phase). Per line: cnt = children with
excess <= s; off = sum over the others with excess > 3,000 of P(x) max(x, 0.05)^-1.07, x = excess / s
(model_chain.json), the README's tie-break. Writes out/signal_TAG_kKICKS.tsv (children) and _lines.tsv (lines ranked)."""
import sys, os, json, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/../landscape')
import lsq
tag, fixes, lbv, K = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
I = lsq.Inst('tai60b'); s = 608215054.05 - lbv
M = json.load(open(H + '/../landscape/out/model_chain.json'))
nf = dict(lsq.parse_fixes(fixes))
fr = [r for r in range(I.n) if r not in nf]; fc = [c for c in range(I.n) if c not in nf.values()]
t = open(H + '/../landscape/truth/tai60b_ils.perm').read().split()
starts = [I.p0, lsq.mirror(I), np.array(t[1:], dtype=np.int32)]
lines = sys.argv[5].split(',') if len(sys.argv) > 5 else [f'r{r}' for r in fr] + [f'c{c}' for c in fc]
memo = {}
def excess(r, v):
    if (r, v) in memo: return memo[(r, v)]
    fix = dict(nf); fix[r] = v; free = np.array([i for i in range(I.n) if i not in fix])
    c, p = min((I.descend(I.adjust(st, fix), fix) for st in starts), key=lambda x: x[0])
    rng = np.random.default_rng(0)
    for _ in range(K):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix)
        if cq < c: c, p = cq, q
    memo[(r, v)] = c - I.best; return memo[(r, v)]
os.makedirs(H + '/out', exist_ok=True)
w = open(f'{H}/out/signal_{tag}_k{K}.tsv', 'w'); w.write('line\tchild\texcess\tx\n')
res = []
for ln in lines:
    k = int(ln[1:]); ch = [(k, v) for v in fc] if ln[0] == 'r' else [(r, k) for r in fr]
    ex = [excess(r, v) for r, v in ch]
    for (r, v), e in zip(ch, ex): w.write(f'{ln}\t{r}:{v}\t{e:.0f}\t{e / s:.3f}\n')
    x = np.array(ex) / s
    cnt = int((x <= 1).sum()); big = x[(np.array(ex) > 3000) & (x > 1)]
    off = float(lsq.model_cost(M, big).sum()) if len(big) else 0.0
    res.append((cnt, off, ln, sorted(x)[:6]))
w.close()
res.sort()
with open(f'{H}/out/signal_{tag}_k{K}_lines.tsv', 'w') as o:
    o.write(f'# {tag} {fixes} lb {lbv:.1f} slack {s:.0f} kicks {K}\nrank\tline\tcnt\toff\tsmallest_x\n')
    for i, (cnt, off, ln, xs) in enumerate(res):
        o.write(f'{i+1}\t{ln}\t{cnt}\t{off:.3f}\t' + ' '.join(f'{v:.2f}' for v in xs) + '\n')
print(tag, 'slack', round(s), 'top lines:', [(ln, cnt, round(off, 2)) for cnt, off, ln, _ in res[:5]])
