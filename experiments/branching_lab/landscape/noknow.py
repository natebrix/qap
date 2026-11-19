#!/usr/bin/env python3
"""noknow.py INST SEEDS KICKS: a search that does not know the best known solution: ILS (kick 6 exchanges, descent,
keep if better) from random starts, SEEDS runs of KICKS kicks. Writes truth/INST_ils.perm (cost, then the
permutation) for the best found; noknow_chain.py uses it in place of p0."""
import sys, time, numpy as np, lsq
inst, S, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
I = lsq.Inst(inst); best = None; t0 = time.time()
for s in range(S):
    rng = np.random.default_rng(100 + s); c, p = I.descend(rng.permutation(I.n).astype(np.int32), {})
    for it in range(K):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(I.n, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, {})
        if cq < c: c, p = cq, q
    print(inst, 'seed', s, 'excess', c - I.best, f'{time.time() - t0:.0f}s', flush=True)
    if best is None or c < best[0]: best = (c, p)
open(f'{lsq.H}/truth/{inst}_ils.perm', 'w').write(f'{int(best[0])} ' + ' '.join(map(str, best[1])) + '\n')
