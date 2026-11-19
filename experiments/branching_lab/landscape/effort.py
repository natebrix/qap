#!/usr/bin/env python3
"""effort.py RUN [ND] [KICKS] [TAG]: does more search effort change the excess? For every branched (B) and cut-phase (C)
child of a branched node and ND (default 400) random first-solve (D) children, an iterated local search: start
from the best of the p0-adjusted (and, tai60b, mirror-adjusted) descents, then KICKS (default 300) times: a random
kick of 6 exchanges among the free rows of the current best, a descent, keep it if better. Records the excess
after 0, 30, 100, 300 kicks (ex_ils0, ex_ils30, ...). Also the no-knowledge version from a random start
(ex_rils*), the same kicks. Writes out/effortTAG_RUN.tsv (fixes and those columns)."""
import sys, time, numpy as np, pandas as pd, lsq

run = sys.argv[1]; ND = int(sys.argv[2]) if len(sys.argv) > 2 else 400; KICKS = int(sys.argv[3]) if len(sys.argv) > 3 else 300
TAG = sys.argv[4] if len(sys.argv) > 4 else ''
inst = run.split('-')[0].split('_')[0]
I = lsq.Inst(inst)
k = pd.read_pickle(f'{lsq.H}/out/kid_{run}.pkl')
rng = np.random.default_rng(7)
NBC = int(sys.argv[5]) if len(sys.argv) > 5 else 0   # >0: a random sample of that many B/C children (tai60b), weighted in effort_score
bc = k[k.out.isin(['B', 'C', 'E'])]
if NBC and len(bc) > NBC: bc = bc.sample(NBC, random_state=7)
sel = pd.concat([bc, k[k.out == 'D'].sample(min(ND, int((k.out == 'D').sum())), random_state=7)])
starts = [I.p0]
if inst == 'tai60b':
    mir = lsq.mirror(I)
    starts.append(mir)
marks = [m for m in (0, 30, 100, 300, 1000, 3000) if m <= KICKS]

def ils(p, c, fix, free):
    out = {}
    if 0 in marks: out[0] = c
    for it in range(1, KICKS + 1):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix)
        if cq < c: c, p = cq, q
        if it in marks: out[it] = c
    return out

rows = []; t0 = time.time()
for n_, (_, r) in enumerate(sel.iterrows()):
    fix = dict(lsq.parse_fixes(r.fixes)); free = np.array([i for i in range(I.n) if i not in fix])
    best = min((I.descend(I.adjust(s, fix), fix) for s in starts), key=lambda x: x[0])
    o = ils(best[1], best[0], fix, free)
    rp = I.adjust(rng.permutation(I.n).astype(np.int32), fix); c, p = I.descend(rp, fix)
    o2 = ils(p, c, fix, free)
    rows.append([r.fixes] + [o[m] - I.best for m in marks] + [o2[m] - I.best for m in marks])
    if n_ % 200 == 0: print(run, n_, len(sel), f'{time.time() - t0:.0f}s', flush=True)
pd.DataFrame(rows, columns=['fixes'] + [f'ex_ils{m}' for m in marks] + [f'ex_rils{m}' for m in marks]).to_csv(
    f'{lsq.H}/out/effort{TAG}_{run}.tsv', sep='\t', index=False)
print(run, 'done', len(sel), f'{time.time() - t0:.0f}s')
