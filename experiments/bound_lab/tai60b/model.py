#!/usr/bin/env python3
"""model.py ORDER SLACK [CAP]: a model of the tai60b tree from the near-optimal solutions around the best known
permutation p (QAPSOLN/tai60b.sln). [e] throughout: it is an estimate, not a search.

A node (a prefix of fixes along a fixed branching order of rows) is taken to SURVIVE (branch) when the best
completion found for it (2-swap best-improvement descent from its parent's completion, the fixes held) costs less
than best known + s(d), where s(d) is the cut phase's slack at depth d: the shortfall of the on-path node at the
end of its phase, measured for the depths the chain reached (SLACK = comma list for d = 1, 2, ...), and continued
below by the ratio RATIO = last argument of SLACK after '/', e.g. 2.72,1.58,1.07,0.73/0.7 (in millions).
ORDER: 'path' = the rows the search branched on down the on-path chain (PATHROWS below), then the remaining rows
cluster by cluster (smallest cluster first); 'rand:SEED' = PATHROWS then random.
Counts branch nodes per depth (the on-path node included); a level with more than CAP branch nodes (default
200000) is sampled down to CAP and weighted.
The descent finds an upper bound on a subtree's optimum, so 'survives' errs towards too few survivors where the
descent is poor and towards too many where the true bound is tighter off the path than on it."""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/regions')
import gen_regions as g
Q = __file__.rsplit('/', 4)[0]
n, A, B = g.read_qaplib(f'{Q}/qaplib/tai60b.dat')
Bs = (B + B.T) / 2.0
t = open(f'{Q}/QAPSOLN/tai60b.sln').read().split()
best, p0 = int(t[1]), np.array(list(map(int, t[2:2 + n]))) - 1
PATHROWS = [int(x) for x in open(__file__.rsplit('/', 1)[0] + '/pathrows.txt').read().split()]
def cost(p): return float((A * Bs[np.ix_(p, p)]).sum())
assert abs(cost(p0) - best) < 0.5
def deltas(p):
    Bp = Bs[np.ix_(p, p)]
    M = A @ Bp.T; s = (A * Bp).sum(1)
    T = M + M.T - s[:, None] - s[None, :]
    return 2 * (T + 2 * A * Bp)
# self-test of the delta formula
_p = p0.copy(); _D = deltas(_p); _q = _p.copy(); _q[3], _q[17] = _q[17], _q[3]
assert abs(cost(_q) - cost(_p) - _D[3, 17]) < 1e-3, (cost(_q) - cost(_p), _D[3, 17])
def descent(p, c, free):
    while True:
        D = deltas(p)[np.ix_(free, free)]
        k = int(np.argmin(D)); a, b = divmod(k, len(free))
        if D[a, b] > -0.5: return p, c
        i, j = free[a], free[b]; p[i], p[j] = p[j], p[i]; c += D[a, b]
order_spec, slack_spec = sys.argv[1], sys.argv[2]
cap = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
sl, ratio = slack_spec.split('/'); sl = [float(x) * 1e6 for x in sl.split(',')]; ratio = float(ratio)
def s(d): return sl[d - 1] if d <= len(sl) else sl[-1] * ratio ** (d - len(sl))
cls = sorted(g.top_clusters(A, 6), key=len)
rest = [r for c in cls for r in sorted(c) if r not in PATHROWS]
if order_spec.startswith('rand'):
    np.random.default_rng(int(order_spec.split(':')[1])).shuffle(rest)
order = PATHROWS + rest
branch = [0] * (n + 1); kids = [0] * (n + 1); phasef = [0] * (n + 1); total = 0
FS = 1.6   # first-solve slack / phase slack (on-path nodes, depths 1-6: 1.3-2.2)
branch[0] = 1; kids[1] = n
# level by level; above W branch nodes at a level a uniform sample of W is expanded and weighted (cap = W)
level = [(p0.copy(), float(best))]; wt = 1.0; capped = False
rng = np.random.default_rng(7)
for d in range(1, n - 1):
    if not level: break
    if len(level) > cap:
        idx = rng.choice(len(level), cap, replace=False); wt *= len(level) / cap; level = [level[i] for i in idx]; capped = True
    r = order[d - 1]; free = order[d:]; sd = s(d); nxt = []
    for q, c in level:
        D = deltas(q)
        for j in free + [r]:
            if j == r: qq, cc = q, c
            else:
                if c + D[r, j] - best > 4 * sd + 1e6: continue     # screen: far children are not descended
                qq = q.copy(); qq[r], qq[j] = qq[j], qq[r]
                qq, cc = descent(qq, c + D[r, j], free)
            ex = cc - best
            if sd <= ex < FS * sd: phasef[d] += wt
            if ex < sd:
                branch[d] += wt; kids[d + 1] += wt * (n - d); nxt.append((qq.copy(), cc))
    level = nxt
branch = [int(round(x)) for x in branch]; kids = [int(round(x)) for x in kids]; phasef = [int(round(x)) for x in phasef]
print(f'order {order_spec} slack {slack_spec}: nodes by depth')
print('depth  s(d)      bounded  branch  fathom in a phase  iterations  CPU-h')
# cost [e]: a branching node 110 + 800 iterations, a phase-fathomed one 110 + 400, a first-solve fathom 45;
# 1.3 s per iteration at m = 59, scaled by (m/59)^4.5 (measured 0.93 at m = 56, 0.72 at m = 54, 0.46 at m = 48, 0.43 at m = 46)
IT = 0.0; CPU = 0.0
for d in range(0, n - 1):
    nb = kids[d] if d else 1
    it = branch[d] * 910 + phasef[d] * 510 + max(0, nb - branch[d] - phasef[d]) * 45
    if d == 0: it = 3470
    cpu = it * 1.3 * ((n - d) / 59.0) ** 4.5 / 3600
    IT += it; CPU += cpu
    if branch[d] > 1 or d < 10: print(f'{d:3d}  {s(d) if d else float("nan"):10.0f}  {nb:7d} {branch[d]:7d}  {phasef[d]:7d}  {it:10.0f}  {cpu:7.2f}')
print(f'total: branch nodes {sum(branch)}{" (levels above the cap were sampled)" if capped else ""}, nodes bounded {sum(kids) + 1}, iterations {IT:.3g}, core-hours on this box {CPU:.1f}')
