#!/usr/bin/env python3
"""effort.py: are the misses of the strong arms a search failure? For chosen children of 10:39 (c47's and c54's
borderline ones), an ILS of K kicks (6 random exchanges + a descent) from the best of the p0 / mirror / ILS-start
descents, 3 seeds; prints the excess after 30, 300, 3,000 kicks. Writes out/effort.tsv."""
import sys, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/../landscape')
import lsq
I = lsq.Inst('tai60b'); NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39'; S = 608215054.05 - 607809132.8147
nf = dict(lsq.parse_fixes(NODE))
t = open(H + '/../landscape/truth/tai60b_ils.perm').read().split()
starts = [I.p0, lsq.mirror(I), np.array(t[1:], dtype=np.int32)]
kids = sys.argv[1:] or ['45:47', '29:47', '31:47', '28:47', '26:47', '40:47', '35:47', '46:47', '21:54', '15:54', '25:54']
w = open(H + '/out/effort.tsv', 'w'); w.write('child\tseed\tex30\tex300\tex3000\n')
for ch in kids:
    r, v = map(int, ch.split(':')); fix = dict(nf); fix[r] = v
    free = np.array([i for i in range(I.n) if i not in fix])
    c0, p0 = min((I.descend(I.adjust(s, fix), fix) for s in starts), key=lambda x: x[0])
    for seed in range(3):
        rng = np.random.default_rng(seed); c, p = c0, p0.copy(); rec = {}
        for k in range(1, 3001):
            q = p.copy()
            for _ in range(6):
                i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
            cq, q = I.descend(q, fix)
            if cq < c: c, p = cq, q
            if k in (30, 300, 3000): rec[k] = c - I.best
        w.write(f'{ch}\t{seed}\t{rec[30]:.0f}\t{rec[300]:.0f}\t{rec[3000]:.0f}\n'); w.flush()
    print(ch, 'x at 30/300/3000 kicks (best of 3 seeds):', *(f'{min(float(l.split()[2+i]) for l in open(H + "/out/effort.tsv") if l.startswith(ch + chr(9)))/S:.2f}' for i in range(3)), flush=True)
