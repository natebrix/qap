#!/usr/bin/env python3
"""effort_all.py: every child of c54, c47, r31 at 10:39, an ILS of 3,000 kicks (one seed, 0) from the best of the
p0 / mirror / ILS-start descents; the excess after 30, 300 and 3,000 kicks. Writes out/effort_all.tsv."""
import sys, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/../landscape')
import lsq
I = lsq.Inst('tai60b'); NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39'
nf = dict(lsq.parse_fixes(NODE))
fr = [r for r in range(I.n) if r not in nf]; fc = [c for c in range(I.n) if c not in nf.values()]
t = open(H + '/../landscape/truth/tai60b_ils.perm').read().split()
starts = [I.p0, lsq.mirror(I), np.array(t[1:], dtype=np.int32)]
w = open(H + '/out/effort_all.tsv', 'w'); w.write('line\tchild\tex30\tex300\tex3000\n')
for ln in sys.argv[1:] or ['c47', 'c54', 'r31']:
    k = int(ln[1:]); ch = [(k, v) for v in fc] if ln[0] == 'r' else [(r, k) for r in fr]
    for r, v in ch:
        fix = dict(nf); fix[r] = v; free = np.array([i for i in range(I.n) if i not in fix])
        c, p = min((I.descend(I.adjust(s, fix), fix) for s in starts), key=lambda x: x[0])
        rng = np.random.default_rng(0); rec = {}
        for kk in range(1, 3001):
            q = p.copy()
            for _ in range(6):
                i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
            cq, q = I.descend(q, fix)
            if cq < c: c, p = cq, q
            if kk in (30, 300, 3000): rec[kk] = c - I.best
        w.write(f'{ln}\t{r}:{v}\t{rec[30]:.0f}\t{rec[300]:.0f}\t{rec[3000]:.0f}\n'); w.flush()
