"""Near-optimal solutions around the best known one (QAPSOLN/tai50b.sln), value = sum A[i,k] B[p(i),p(k)]:
  python3 landscape.py            2-swap neighbours by excess over 458,821,517; per depth-1 child 31:x the best
                                  solution found with p(31) = x (swap/3-cycle from the best known, then 2-swap descent
                                  with p(31) held) -> logs/landscape.txt, logs/perm_31_X.txt for the probes
A node whose best completion is delta above the best known value cannot fathom while the bound's gap exceeds delta."""
import sys, numpy as np
n = 50
d = open('/home/user/qap/qaplib/tai50b.dat').read().split()
A = np.array(d[1:1+n*n], dtype=np.int64).reshape(n, n); B = np.array(d[1+n*n:1+2*n*n], dtype=np.int64).reshape(n, n)
BK = np.array([5,23,14,12,21,4,29,40,13,8,16,19,25,18,46,11,15,37,17,30,47,20,33,26,36,48,22,10,3,31,2,7,0,32,6,42,24,28,44,27,39,35,43,45,1,34,49,38,9,41])
val = lambda p: int((A * B[np.ix_(p, p)]).sum())
V0 = val(BK); assert V0 == 458821517
def swapped(p, i, j): q = p.copy(); q[i], q[j] = q[j], q[i]; return q
sw = sorted((val(swapped(BK, i, j)) - V0, i, j) for i in range(n) for j in range(i + 1, n))
print("2-swap neighbours of the best known solution: %d; none better: %s" % (len(sw), sw[0][0] >= 0))
for t in (1e5, 3e5, 1e6, 2e6, 3e6, 5e6): print("  within %9.0f: %d" % (t, sum(1 for s in sw if s[0] <= t)))
print("  ten nearest:", [(s[0], "%d<->%d" % (s[1], s[2])) for s in sw[:10]])
def descent(p, hold):
    v = val(p)
    while True:
        best = (0, None)
        for i in range(n):
            for j in range(i + 1, n):
                if i in hold or j in hold: continue
                dv = val(swapped(p, i, j)) - v
                if dv < best[0]: best = (dv, (i, j))
        if best[1] is None: return p, v
        p = swapped(p, *best[1]); v += best[0]
print("\nbest solution found with p(31) = x (the depth-1 children of the root, which branches on row 31):")
rows = []
for x in range(n):
    if x == BK[31]: rows.append((0, x, BK)); continue
    j = int(np.where(BK == x)[0][0]); p, v = descent(swapped(BK, 31, j), {31}); rows.append((v - V0, x, p))
for dv, x, p in sorted(rows, key=lambda r: r[0]):
    print("  31:%-2d  +%d" % (x, dv)); np.savetxt("logs/perm_31_%d.txt" % x, p[None], fmt="%d")
