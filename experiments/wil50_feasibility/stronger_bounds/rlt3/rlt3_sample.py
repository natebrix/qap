#!/usr/bin/env python3
"""Violation sample for projected RLT-k consistency on a dumped DNN iterate.

For k facilities (k = 3: the triple consistency of step 8, projected RLT2;
k = 4: quadruple consistency, projected RLT3), the k(k-1)/2 pairwise blocks of
Y must be the marginals of one joint placement of the k facilities. To keep
the separation LP small, locations outside the s heaviest for these
facilities are merged into one "elsewhere" cell; an inequality valid on the
merged placements is valid on Y (it is a sum of Y entries), so the violation
found is a lower bound on the true one. Potentials G in [-1, 1] as Astra's
separator; after rounding to 1/1024 the right-hand side is the exact maximum
of sum G over all merged placements, so every reported violation is that of
a valid inequality, whatever the LP did.

  python3 rlt3_sample.py PREFIX --fixes SPEC --k 4 --count 100 [--seed 1] [--n 50] [--from-triples '../triples/n13.triplet_*.json']
PREFIX.Y.bin as dumped by dnntest -dumpY. Prints one line per sampled set
(facilities, violation, LP seconds; for k = 4 also the best sub-triple's
violation and the excess over it) and a summary."""
import argparse, itertools, json, random, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

def node_m(spec, n):
    fixes = [] if spec == '-' else [tuple(map(int, t.split(':'))) for t in spec.split(',')]
    return n - len(fixes)

def all_placements(k, m):
    """All injective k-tuples of locations, as an int16 array (m^k candidates masked)."""
    grids = np.indices((m,) * k, dtype=np.int16).reshape(k, -1)
    ok = np.ones(grids.shape[1], dtype=bool)
    for a in range(k):
        for c in range(a + 1, k):
            ok &= grids[a] != grids[c]
    return grids[:, ok].T

_pl = {}
def placements_for(k, m):
    if (k, m) not in _pl: _pl[(k, m)] = all_placements(k, m)
    return _pl[(k, m)]

def violation(Y, m, fac, seconds, max_rounds=60, add=3000, init=20000, seed=0):
    """Exact separation of a projected RLT-(k-1) consistency inequality <G, Y> <= 0 for
    the k facilities fac, by row generation over the injective placements: the LP
    maximizes <G, Y> over G in [-1, 1] with sum of G <= 0 on the placements seen so
    far; the most violated placements (a vectorized scan of all of them) are added
    until none is. The reported violation is for the 1/1024-rounded G against the
    exact maximum over all placements, so it is that of a valid inequality."""
    k = len(fac); pairs = list(itertools.combinations(range(k), 2)); nb = len(pairs)
    ids = [1 + np.arange(m) * m + f for f in fac]
    blocks = [Y[np.ix_(ids[a], ids[c])] for a, c in pairs]
    cvec = np.concatenate([B.ravel() for B in blocks])
    bounds = [(-1.0, 1.0)] * (nb * m * m)
    for b in range(nb):                                   # same-location entries appear in no placement
        for a in range(m): bounds[b * m * m + a * m + a] = (0.0, 0.0)
    pl = placements_for(k, m)
    colidx = np.stack([b * m * m + pl[:, a].astype(np.int64) * m + pl[:, c] for b, (a, c) in enumerate(pairs)], 1)  # P x nb
    rng = np.random.default_rng(seed)
    rows = rng.choice(len(pl), size=min(init, len(pl)), replace=False) if len(pl) > init else np.arange(len(pl))
    t = time.monotonic(); G = None; n_lp = 0
    while True:
        ci = colidx[rows]
        A = coo_matrix((np.ones(ci.size), (np.repeat(np.arange(len(rows)), nb), ci.ravel())), shape=(len(rows), nb * m * m)).tocsr()
        lp = linprog(-cvec, A_ub=A, b_ub=np.zeros(len(rows)), bounds=bounds, method='highs', options={'time_limit': seconds})
        n_lp += 1
        if lp.status != 0: return float('nan'), time.monotonic() - t, n_lp
        G = lp.x
        sums = G[colidx].sum(1)                            # every placement
        bad = np.flatnonzero(sums > 1e-7)
        if len(bad) == 0 or len(rows) == len(pl) or n_lp >= max_rounds: break
        worst = bad[np.argsort(-sums[bad])[:add]]
        rows = np.union1d(rows, worst)
    G = np.round(G * 1024) / 1024
    rhs = G[colidx].sum(1).max()
    return float(cvec @ G - rhs), time.monotonic() - t, n_lp

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefix'); ap.add_argument('--fixes', default='-'); ap.add_argument('--n', type=int, default=50)
    ap.add_argument('--k', type=int, default=4); ap.add_argument('--count', type=int, default=100)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--from-triples', default=None, help='glob of step-8 triplet json files: extend their most violated triples by a random 4th facility')
    ap.add_argument('--seconds', type=float, default=120); ap.add_argument('--json', default=None)
    a = ap.parse_args()
    m = node_m(a.fixes, a.n); N = m * m + 1
    Y = np.fromfile(a.prefix + '.Y.bin').reshape(N, N)
    rng = random.Random(a.seed); seen = set(); recs = []
    known = {}                                            # step-8 triple violations, by local facilities
    seeds = []
    if a.from_triples:
        import glob
        for j in sorted(glob.glob(a.from_triples)):
            d = json.load(open(j))
            if d.get('success'): known[tuple(sorted(d['facilities_local']))] = d['rounded_cut_violation']
        seeds = [t for t, v in sorted(known.items(), key=lambda kv: -kv[1])]
    while len(recs) < a.count:
        if seeds and a.k == 4 and len(recs) < a.count // 2:   # first half: extend the most violated known triples
            t = seeds[len(recs) % len(seeds)]
            fac = tuple(sorted(t + (rng.choice([f for f in range(m) if f not in t]),)))
        else:
            fac = tuple(sorted(rng.sample(range(m), a.k)))
        if fac in seen: continue
        seen.add(fac)
        v, dt, nl = violation(Y, m, fac, a.seconds)
        r = dict(facilities=fac, violation=v, seconds=round(dt, 1), lps=nl)
        if a.k == 4:
            subs = [known.get(t) if t in known else violation(Y, m, t, a.seconds)[0] for t in itertools.combinations(fac, 3)]
            r['best_triple'] = max(subs); r['excess'] = v - max(subs)
        recs.append(r)
        extra = f" best3 {r['best_triple']:.4f} excess {r['excess']:+.4f}" if a.k == 4 else ""
        print(f"{','.join(map(str, fac))} viol {v:.4f} {nl} LPs {dt:.0f}s{extra}", flush=True)
    vs = np.array([r['violation'] for r in recs]); vs = vs[~np.isnan(vs)]
    pct = lambda z, p: float(np.percentile(z, p)) if len(z) else float('nan')
    print(f"SUMMARY k={a.k} node={a.prefix} m={m} sets={len(recs)} violated(>1e-3)={int((vs > 1e-3).sum())} "
          f"median={pct(vs,50):.4f} p90={pct(vs,90):.4f} max={vs.max() if len(vs) else float('nan'):.4f} sum={vs.sum():.3f}")
    if a.k == 4:
        ex = np.array([r['excess'] for r in recs]); b3 = np.array([r['best_triple'] for r in recs])
        print(f"SUMMARY k=4 excess over best sub-triple: median={pct(ex,50):+.4f} p90={pct(ex,90):+.4f} max={ex.max():+.4f}; "
              f"best sub-triple median={pct(b3,50):.4f}")
    if a.json: json.dump(dict(args=vars(a), m=m, records=recs), open(a.json, 'w'), indent=1)

if __name__ == '__main__':
    main()
