#!/usr/bin/env python3
"""Gate for the Jünger-Kaibel 1-box (count) inequalities with LARGE facility
sets on a dumped wil50 iterate.

For a location set P (side 0) and a facility set Q (side 1), the count
q = #{facilities of Q placed in P} is an integer, so for every integer b
  E[(q - b)(q - b - 1)] >= 0  <=>  2b sum_{a in P, j in Q} x_aj
                                   - sum_{a != c in P, j != k in Q} Y_(aj),(ck) <= b(b+1),
violated iff Var(q) < f(1 - f) with f = frac(E q), b = floor(E q). Step 7
(5-point count inequalities) is the |Q| = 5 case on the threshold regions;
the triangles are |Q| = 3. Jünger-Kaibel (MP 2001) found |Q| up to n/2,
b in {2, 3} the most effective at the LP level.

Search: for each region P (13 threshold regions, plus all grid rectangles
[r0,r1) x [c0,c1) if --rects), local search over Q (add/drop/swap moves from
several random starts, best-improvement on the violation), sizes 3..qmax.
Reports per region the max violation by |Q| and the top cuts.

  python3 -I box_gate.py PREFIX --fixes SPEC [--rects] [--starts 40] [--qmax 24]
"""
import argparse, json, time
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('prefix'); ap.add_argument('--fixes', default='-'); ap.add_argument('--n', type=int, default=50)
ap.add_argument('--grid', default='5:10'); ap.add_argument('--rects', action='store_true')
ap.add_argument('--starts', type=int, default=40); ap.add_argument('--qmax', type=int, default=24)
ap.add_argument('--seed', type=int, default=1); ap.add_argument('--out', default=None)
ap.add_argument('--top', type=int, default=0, help='dump the top cuts (P, Q, b, violation) to --out')
a = ap.parse_args()
n = a.n
fixes = [] if a.fixes == '-' else [tuple(map(int, t.split(':'))) for t in a.fixes.split(',')]
free_loc = [i for i in range(n) if i not in {x for x, _ in fixes}]
free_fac = [j for j in range(n) if j not in {y for _, y in fixes}]
m = len(free_loc); N = m * m + 1
Y = np.fromfile(a.prefix + '.Y.bin').reshape(N, N)
Y = (Y + Y.T) / 2
rows, cols = map(int, a.grid.split(':'))
Yb = Y[1:, 1:].reshape(m, m, m, m)       # Yb[a, j, c, k] = Y_(aj),(ck)
X = Y[0, 1:].reshape(m, m)

regions = []
for t in range(1, rows):
    regions.append((f'row<{t}', [k for k, loc in enumerate(free_loc) if loc // cols < t]))
for t in range(1, cols):
    regions.append((f'col<{t}', [k for k, loc in enumerate(free_loc) if loc % cols < t]))
if a.rects:
    for r0 in range(rows):
        for r1 in range(r0 + 1, rows + 1):
            for c0 in range(cols):
                for c1 in range(c0 + 1, cols + 1):
                    if (r0 == 0 and r1 == rows) or (c0 == 0 and c1 == cols):
                        continue      # the threshold regions and their complements
                    S = [k for k, loc in enumerate(free_loc)
                         if r0 <= loc // cols < r1 and c0 <= loc % cols < c1]
                    if 2 <= len(S) <= m - 2:
                        regions.append((f'rect{r0}{r1}x{c0}{c1}', S))
regions = [(nm, S) for nm, S in regions if 0 < len(S) < m]


def moments(S):
    S = np.array(S)
    p = X[S, :].sum(0)
    Wf = Yb[np.ix_(S, np.arange(m), S, np.arange(m))].sum(axis=(0, 2))
    W = Wf - np.diag(np.diag(Wf))
    return p, W


def viol(mask, p, W):
    mu = p @ mask
    Eq2 = mu + mask @ W @ mask
    b = np.floor(mu)
    return -(Eq2 - (2 * b + 1) * mu + b * (b + 1)), b


rng = np.random.default_rng(a.seed)
res = {'regions': {}, 'top': []}
allcuts = []
t0 = time.time()


def g(mask, p, W, b):
    """violation of the b-inequality: 2b p.m - m'Wm - b(b+1)  (= -E[(q-b)(q-b-1)])"""
    return 2 * b * (p @ mask) - mask @ W @ mask - b * (b + 1)


def local(mask, p, W, b):
    v = g(mask, p, W, b)
    while True:
        # marginal gains of flipping each element (W symmetric, zero diagonal)
        Wm = W @ mask
        gain_add = 2 * b * p - 2 * Wm            # adding j (mask_j = 0)
        gain_drop = -2 * b * p + 2 * Wm          # dropping j (mask_j = 1); W_jj = 0
        gains = np.where(mask > 0, gain_drop, gain_add)
        sz = mask.sum()
        gains = np.where((mask > 0) & (sz <= 3), -np.inf, gains)
        gains = np.where((mask == 0) & (sz >= a.qmax), -np.inf, gains)
        j = int(np.argmax(gains))
        if gains[j] <= 1e-10:
            return mask, v
        mask = mask.copy(); mask[j] = 1 - mask[j]
        v = g(mask, p, W, b)


for nm, S in regions:
    p, W = moments(S)
    best_by_b = {}
    for b in range(1, 8):
        for st in range(a.starts):
            if st == 0:
                mask = np.zeros(m)
                order = np.argsort(-p)
                mask[order[:b + 2]] = 1
            else:
                k0 = int(rng.integers(b + 1, min(a.qmax, 3 * b + 6) + 1))
                w = p * (1 - p) + 1e-3
                mask = np.zeros(m); mask[rng.choice(m, k0, replace=False, p=w / w.sum())] = 1
            mask, v = local(mask, p, W, b)
            if v > best_by_b.get(b, (-np.inf,))[0]:
                best_by_b[b] = (float(v), [int(j) for j in np.where(mask)[0]])
            allcuts.append((float(v), nm, [int(j) for j in np.where(mask)[0]], b, [int(x) for x in S]))
    mx = max(x[0] for x in best_by_b.values())
    res['regions'][nm] = {'size': len(S), 'max': mx,
                          'by_b': {b: (best_by_b[b][0], len(best_by_b[b][1])) for b in best_by_b}}
    print(f"{nm:12s} |P|={len(S):2d} max {mx:.4f}  by b (viol,|Q|): " +
          " ".join(f"{b}:({best_by_b[b][0]:.3f},{len(best_by_b[b][1])})" for b in sorted(best_by_b)), flush=True)
allcuts.sort(key=lambda r: -r[0])
seen = set(); top = []
for c in allcuts:
    key = (c[1], tuple(c[2]), c[3])
    if key in seen: continue
    seen.add(key); top.append(c)
res['top'] = top[:max(a.top, 50)]
res['m'] = m; res['fixes'] = a.fixes
print("top violations:", [round(c[0], 4) for c in top[:20]])
print("time %.1fs" % (time.time() - t0))
if a.out:
    json.dump(res, open(a.out, 'w'), indent=1)
print("RESULT max %.4f" % (top[0][0] if top else 0))
