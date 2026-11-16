#!/usr/bin/env python3
"""cutset.py IN.state OUT.state K SEED [F:L ...]: keep K open nodes of a bnbdrv state file drawn uniformly at
random (seed SEED), plus the nodes whose permutation has entry F equal to L for each F:L given (always kept,
not counted in K). Prints the fixes of each kept node (random ones first, then the forced ones)."""
import random, sys
src, dst, K, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
forced = [tuple(map(int, a.split(':'))) for a in sys.argv[5:]]
head, nodes = [], []
with open(src) as fh:
    for ln in fh:
        if ln.startswith('nodes '): continue
        (nodes if ln.startswith('node ') else head).append(ln)
def fixes(ln, n):
    t = ln.split()
    return [(i, int(t[5 + i])) for i in range(n) if int(t[5 + i]) >= 0]
n = int(next(l for l in head if l.startswith('n ')).split()[1])
idx = list(range(len(nodes)))
fidx = [i for i in idx if any(f in fixes(nodes[i], n) for f in forced)]
rest = [i for i in idx if i not in fidx]
random.Random(seed).shuffle(rest)
keep = rest[:K] + fidx
with open(src) as fh, open(dst, 'w') as w:
    k = -1
    for ln in fh:
        if ln.startswith('nodes '):
            w.write(f'nodes {len(keep)}\n'); continue
        if ln.startswith('node '):
            k += 1
            if k in keep: w.write(ln)
            continue
        w.write(ln)
print(len(nodes), 'open nodes; kept', len(keep))
for i in keep:
    print(('forced ' if i in fidx else 'random ') + ','.join(f'{a}:{b}' for a, b in fixes(nodes[i], n)))
