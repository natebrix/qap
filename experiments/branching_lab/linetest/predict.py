#!/usr/bin/env python3
"""predict.py: the landscape signal's prediction, child by child, at tai60b's depth-7 chain node 10:39
(fixes 44:45,43:58,1:33,6:24,7:19,0:32,10:39) for the lines c54, c47, r31: each child's excess (best completion
found - best known), by the search of ../landscape/chain.py (one 2-exchange descent from each start adjusted to
the child's fixes, then an ILS of 30 kicks of 6 random exchanges), under two sets of starts: 'p0' (the best known
solution and its mirror, the landscape README's main run) and 'ils' (truth/tai60b_ils.perm, the solution found
without knowing p0). Writes out/predict.tsv: start, line, child, ex0, exK, and x = exK / s for production's
slack s = 405,921 (the node's cutoff - its cut-phase bound in tai60b-261117).
The ILS draws differ from chain.py's (which searched every free line in another order), so a child's exK can
differ from the landscape run's by what 30 kicks find; ex0 (descents) is deterministic."""
import sys, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/../landscape')
import lsq
I = lsq.Inst('tai60b')
NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39'
S_PROD = 608215054.05 - 607809132.8147
nf = dict(lsq.parse_fixes(NODE))
fr = [r for r in range(I.n) if r not in nf]; fc = [c for c in range(I.n) if c not in nf.values()]
lines = sys.argv[1:] or ['c54', 'c47', 'r31']
starts = {'p0': [I.p0, lsq.mirror(I)]}
t = open(H + '/../landscape/truth/tai60b_ils.perm').read().split()
starts['ils'] = [np.array(t[1:], dtype=np.int32)]
os.makedirs(H + '/out', exist_ok=True)
w = open(H + '/out/predict.tsv', 'w'); w.write('start\tline\tchild\tex0\texK\tx_prod\n')
for sname, st in starts.items():
    rng = np.random.default_rng(5)
    for ln in lines:
        k = int(ln[1:])
        ch = [(k, v) for v in fc] if ln[0] == 'r' else [(r, k) for r in fr]
        for r, v in ch:
            fix = dict(nf); fix[r] = v
            c, p = min((I.descend(I.adjust(s, fix), fix) for s in st), key=lambda x: x[0]); c0 = c
            free = np.array([i for i in range(I.n) if i not in fix])
            for _ in range(30):
                q = p.copy()
                for _ in range(6):
                    i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
                cq, q = I.descend(q, fix)
                if cq < c: c, p = cq, q
            w.write(f'{sname}\t{ln}\t{r}:{v}\t{c0 - I.best:.0f}\t{c - I.best:.0f}\t{(c - I.best) / S_PROD:.3f}\n')
        w.flush()
w.close()
import collections
d = collections.defaultdict(list)
for l in open(H + '/out/predict.tsv').read().split('\n')[1:]:
    if l: s, ln, ch, e0, ek, x = l.split('\t'); d[(s, ln)].append(float(x))
for (s, ln), xs in d.items():
    xs.sort(); print(s, ln, 'predicted survivors (x <= 1):', sum(x <= 1 for x in xs), ' 1 < x <= 1.5:', sum(1 < x <= 1.5 for x in xs), ' smallest x:', ' '.join(f'{x:.2f}' for x in xs[:12]))
