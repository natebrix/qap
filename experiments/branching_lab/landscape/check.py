"""check.py: the conventions (lsq.py's docstring) against the traces, and ls.c against a slow recomputation."""
import numpy as np, lsq
for inst, run in [('tai60b', 'tai60b-261117'), ('tai50b', 'tai50b-261116'), ('tai50b', 'tai50b-r9'),
                  ('tai35b', None), ('tai40b', None)]:
    I = lsq.Inst(inst)
    path = f'{lsq.WT}/infra/phase1/runs/{run}/nodes.txt.gz' if run else \
        f'{lsq.WT}/experiments/line_signals/runs/{inst}_base/nodes.txt.gz'
    T = lsq.read_trace(path)
    inv = np.argsort(I.p0)
    on = [t for t in T if all(I.p0[r] == v for r, v in t['fx'])]
    onT = [t for t in T if all(inv[r] == v for r, v in t['fx'])]   # the other reading: r a location, v a facility
    print(f'{run or inst}: {len(T)} nodes; consistent with p0 as p[r]=v: {len(on)} nodes, outcomes',
          sorted(set(t['out'] for t in on)), f'deepest {max(t["depth"] for t in on)}; as the inverse: {len(onT)}',
          sorted(set(t['out'] for t in onT)))
# ls.c: descents end at true local minima with the right cost
I = lsq.Inst('tai60b'); rng = np.random.default_rng(0)
for trial in range(5):
    fix = {int(r): int(I.p0[r]) for r in rng.choice(I.n, 10, replace=False)}
    p = I.adjust(rng.permutation(I.n), fix)
    c, q = I.descend(p, fix)
    assert abs(c - I.cost(q)) < 1e-3 * abs(c), (c, I.cost(q))
    free = [i for i in range(I.n) if i not in fix]
    best = min(I.cost(np.where(np.arange(I.n) == i, q[j], np.where(np.arange(I.n) == j, q[i], q))) for i in free for j in free if i < j)
    assert best >= c - 0.5 and all(q[r] == v for r, v in fix.items())
    print('descent ok', c, best - c)
