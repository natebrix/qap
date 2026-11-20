#!/usr/bin/env python3
"""hardsearch.py: is the signal's x for a child an overestimate? For tai50b children r:v of the root, an ILS
that does not start from the incumbent: SEEDS random starts (the fix forced), KICKS kicks of 6 exchanges among
the free rows + descent, keep if better. Prints the best excess over the best known as a multiple of the root's
slack, beside the signal's x (root_children.tsv)."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/qap/experiments/branching_lab/landscape")
import lsq
S, K = int(sys.argv[1]), int(sys.argv[2]); kids = sys.argv[3:]
I = lsq.Inst("tai50b"); slack = 3372320.9
sig = {l.split("\t")[1]: l.split("\t") for l in open("root_children.tsv") if not l.startswith("line")}
for kid in kids:
    r, v = map(int, kid.split(":")); fix = {r: v}; free = [i for i in range(I.n) if i != r]; best = None; t0 = time.time()
    for s in range(S):
        rng = np.random.default_rng(1000 + s)
        c, p = I.descend(I.adjust(rng.permutation(I.n).astype(np.int32), fix), fix)
        for it in range(K):
            q = p.copy()
            for _ in range(6):
                i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
            cq, q = I.descend(q, fix)
            if cq < c: c, p = cq, q
        assert p[r] == v
        if best is None or c < best: best = c
    row = sig[kid]
    print(f"{kid}\tsignal x {row[3]}\tILS x {(best - I.best) / slack:.2f}\t{row[4]}\tshortfall {row[5]}\t{time.time() - t0:.0f}s", flush=True)
