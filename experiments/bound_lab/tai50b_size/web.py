"""A model of the tree's "web" from the landscape [e]:  python3 web.py [MAXOPEN]  -> logs/web.txt
A node that holds a solution delta above the best known value cannot fathom while the bound's gap on its path
exceeds delta. Model: a child at depth d is OPEN (branches) iff the best solution found inside it (the parent's
best solution with the child's assignment swapped in, then 2-swap descent over the free rows) is within g(d) of the
best known value, g(d) = the shortfall measured on the spine at depth d after the phase (runs/spine_D); every
node branches on the line the spine branched on at its depth (the twin probe did: r33, r32, r45, c39).
Output: per depth the number of nodes and of open nodes, for g and for g/2 and 2g. Depths 1-8 use measured g;
depth 9 assumes g(9) = g(8). The local search can miss a better solution, which would add open nodes."""
import sys, numpy as np
n = 50
d = open('/home/user/qap/qaplib/tai50b.dat').read().split()
A = np.array(d[1:1+n*n], dtype=np.float64).reshape(n, n); B = np.array(d[1+n*n:1+2*n*n], dtype=np.float64).reshape(n, n)
BK = np.array([5,23,14,12,21,4,29,40,13,8,16,19,25,18,46,11,15,37,17,30,47,20,33,26,36,48,22,10,3,31,2,7,0,32,6,42,24,28,44,27,39,35,43,45,1,34,49,38,9,41])
val = lambda p: float((A * B[np.ix_(p, p)]).sum())
V0 = val(BK)
def deltas(p):
    """D[u,v] = val(p with rows u and v swapped) - val(p), all pairs"""
    Bp = B[np.ix_(p, p)]
    M1 = A @ Bp.T; M2 = A.T @ Bp
    d1 = np.diag(M1); d2 = np.diag(M2)
    S = M1 + M1.T - d1[:, None] - d1[None, :] + M2 + M2.T - d2[:, None] - d2[None, :]
    dA = np.diag(A); dB = np.diag(Bp)
    # remove the k in {u, v} terms of both sums, then add the 2x2 block's true change
    k1 = (dA[:, None] - A.T) * (Bp.T - dB[:, None]) + (A - dA[None, :]) * (dB[None, :] - Bp)
    k2 = (dA[:, None] - A) * (Bp - dB[:, None]) + (A.T - dA[None, :]) * (dB[None, :] - Bp.T)
    blk = dA[:, None] * dB[None, :] + dA[None, :] * dB[:, None] + A * Bp.T + A.T * Bp - (dA * dB)[:, None] - (dA * dB)[None, :] - A * Bp - A.T * Bp.T
    return S - k1 - k2 + blk
def descent(p, free):
    p = p.copy(); v = val(p); idx = np.array(sorted(free))
    while True:
        D = deltas(p)[np.ix_(idx, idx)]
        a = np.unravel_index(np.argmin(D), D.shape)
        if D[a] >= -0.5: return p, v
        u, w = idx[a[0]], idx[a[1]]; p[u], p[w] = p[w], p[u]; v += D[a]
if __name__ == "__main__":
    rng = np.random.default_rng(0); p = rng.permutation(n); D = deltas(p)
    for _ in range(20):
        u, w = rng.choice(n, 2, replace=False); q = p.copy(); q[u], q[w] = q[w], q[u]
        assert abs(D[u, w] - (val(q) - val(p))) < 1e-3 * abs(val(p)) * 1e-6 + 1.0, (D[u, w], val(q) - val(p))
    LINES = ["r31", "r33", "r32", "r45", "c39", "r30", "c31", "c20", "c10"]          # the spine's branching line at depths 0..8
    G = [None, 3222232, 2652428, 2185798, 1029797, 626571, 199906, 112268, 107349, 107349]   # spine shortfall at depths 1..8; g(9) = g(8) assumed
    MAXOPEN = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    out = []
    for scale in (1.0, 0.5, 2.0):
        level = [(dict(), BK.copy())]          # (fixed rows -> location, best permutation found inside)
        out.append("\ngap scale %.1f:\n| depth | nodes | open (model) | of which on the spine's twins' side (off the spine) | smallest / median delta of the open off-spine nodes |\n|---|---|---|---|---|" % scale)
        for dep in range(0, 9):
            nxt = []; total = 0
            for fixed, pbest in level:
                line = LINES[dep]; free_rows = [r for r in range(n) if r not in fixed]; used = set(fixed.values())
                if line[0] == "r":
                    r = int(line[1:]); kids = [(r, l) for l in range(n) if l not in used] if r not in fixed else []
                else:
                    l = int(line[1:]); kids = [(r, l) for r in free_rows] if l not in used else []
                if not kids:        # the spine's line is already fixed here: use the first free row
                    r = free_rows[0]; kids = [(r, l) for l in range(n) if l not in used]
                total += len(kids)
                for r, l in kids:
                    q = pbest.copy(); j = int(np.where(q == l)[0][0])
                    if j in fixed and j != r: continue          # the location is held by a fixed row (cannot happen: l is free)
                    q[r], q[j] = q[j], q[r]
                    f2 = dict(fixed); f2[r] = l
                    q, v = descent(q, [x for x in range(n) if x not in f2])
                    if v - V0 <= scale * G[dep + 1]: nxt.append((f2, q, v - V0))
            off = [x[2] for x in nxt if not all(BK[r] == l for r, l in x[0].items())]
            out.append("| %d | %d | %d | %d | %s |" % (dep + 1, total, len(nxt), len(off), ("%.0f / %.0f" % (min(off), float(np.median(off)))) if off else "-"))
            print(out[-1], flush=True)
            level = [(f, q) for f, q, _ in nxt]
            if len(level) > MAXOPEN: out.append("(stopped: more than %d open nodes)" % MAXOPEN); break
    open("logs/web.txt", "w").write("\n".join(out) + "\n")
