"""timing.py: CPU of the landscape signal for every child of every free line at one node (m^2 distinct children),
descents only (one per start) and with an ILS of K kicks, at tai50b's and tai60b's depth-1 chain node."""
import time, numpy as np, lsq
for inst, nf in [('tai50b', {31: 7}), ('tai60b', {44: 45})]:
    I = lsq.Inst(inst); rng = np.random.default_rng(1)
    fr = [r for r in range(I.n) if r not in nf]; fc = [c for c in range(I.n) if c not in nf.values()]
    kids = [(r, c) for r in fr for c in fc]
    for K in (0, 30, 100):
        t = time.process_time(); sub = kids if K == 0 else kids[::10]
        for r, c in sub:
            fix = dict(nf); fix[r] = c
            cc, p = I.descend(I.adjust(I.p0, fix), fix); free = np.array([i for i in range(I.n) if i not in fix])
            for _ in range(K):
                q = p.copy()
                for _ in range(6):
                    i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
                cq, q = I.descend(q, fix)
                if cq < cc: cc, p = cq, q
        dt = (time.process_time() - t) * len(kids) / len(sub)
        print(f'{inst} m={len(fr)}: {len(kids)} children, one start, {K} kicks: {dt:.1f} CPU s for the node (all lines)', flush=True)
