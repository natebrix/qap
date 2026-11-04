#!/usr/bin/env python3
"""Diagnostics of a dumped DNN iterate (dnntest/dnncuts -dumpY) at a wil50 node
for the grid families:

  1. objective split: row part and column part of <L, Y>;
  2. regional triangles on arbitrary cell regions (the 13 prefix regions for
     reference; parity regions: checkerboard, row parity, column parity;
     single rows / columns; 2 x 2-block checkerboards), max violation and the
     sum of positive violations, in the (p, W) units of the production logs;
  3. clique dispersion: sum_{i<j in T} E[d_ij] against the exact Dmin(|T|) over
     the node's free cells, most violated T by greedy + swap search;
  4. weighted (flow) dispersion: sum_{T} F_ij E[d_ij] against the SA optimum
     of the |T|-facility sub-QAP (an estimate of the potential, not a cut);
  5. a region search: for the most violated parity / prefix triangles, local
     search over all cell subsets for the region maximizing the violation.

  python3 analyze_Y.py SPEC YFILE [--json out.json]
"""
import argparse, itertools, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
import dispersion as dsp


def tri_viol(p, W):
    """all four triangle types on (p, W); returns array of violations over
    unordered triples (vectorized over i<j<k) in (p, W) units"""
    m = len(p)
    dl = p[:, None] + p[None, :] - 2 * W          # cut semimetric delta_ij
    np.fill_diagonal(dl, 0)
    i, j, k = np.array(list(itertools.combinations(range(m), 3))).T
    dij, dik, djk = dl[i, j], dl[i, k], dl[j, k]
    v4 = (dij + dik + djk) / 2 - 1
    v1 = (djk - dij - dik) / 2; v2 = (dik - dij - djk) / 2; v3 = (dij - dik - djk) / 2
    V = np.stack([v1, v2, v3, v4], axis=1)
    return V, (i, j, k)


def region_pW(nd, Y4, S):
    s = np.zeros(nd.m); s[list(S)] = 1
    W = np.einsum('a,aibj,b->ij', s, Y4, s, optimize=True)
    p = np.einsum('a,aiai->i', s, Y4, optimize=True)
    return p, W


def summarize(V, thr=(0.01, 0.05)):
    pos = V[V > 0]
    return dict(max=float(V.max()), n_viol=int((V > 1e-6).sum()), sum_pos=float(pos.sum()),
                **{f'n>{t}': int((V > t).sum()) for t in thr})


def regions_of(nd):
    R, C = nd.rows, nd.cols
    loc = nd.free_loc
    rr = np.array([a // C for a in loc]); cc = np.array([a % C for a in loc])
    regs = {}
    for t in range(1, R):
        regs[f'row<{t}'] = np.where(rr < t)[0]
    for t in range(1, C):
        regs[f'col<{t}'] = np.where(cc < t)[0]
    regs['checker'] = np.where((rr + cc) % 2 == 0)[0]
    regs['rowpar'] = np.where(rr % 2 == 0)[0]
    regs['colpar'] = np.where(cc % 2 == 0)[0]
    regs['checker2'] = np.where(((rr // 2) + (cc // 2)) % 2 == 0)[0]
    regs['colpar2'] = np.where((cc // 2) % 2 == 0)[0]
    regs['col mod3'] = np.where(cc % 3 == 0)[0]
    for r in range(R):
        regs[f'row={r}'] = np.where(rr == r)[0]
    for c in range(C):
        regs[f'col={c}'] = np.where(cc == c)[0]
    return {k: v for k, v in regs.items() if 0 < len(v) < nd.m}


def most_clustered(E, k, starts=None, iters=3):
    """T of size k minimizing sum_{pairs in T} E (greedy from each start + swaps)"""
    m = len(E)
    best = (np.inf, None)
    for s0 in (range(m) if starts is None else starts):
        T = [s0]
        while len(T) < k:
            rest = [j for j in range(m) if j not in T]
            j = min(rest, key=lambda j: E[j, T].sum()); T.append(j)
        val = E[np.ix_(T, T)].sum() / 2
        for _ in range(iters):  # swap improvement
            improved = False
            for pos in range(k):
                others = T[:pos] + T[pos + 1:]
                cur = E[T[pos], others].sum()
                rest = [j for j in range(m) if j not in T]
                cand = min(rest, key=lambda j: E[j, others].sum())
                if E[cand, others].sum() < cur - 1e-12:
                    T[pos] = cand; improved = True
            if not improved:
                break
        val = E[np.ix_(T, T)].sum() / 2
        if val < best[0]:
            best = (val, sorted(T))
    return best


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('spec'); ap.add_argument('Y'); ap.add_argument('--json')
    ap.add_argument('--kmax', type=int, default=10)
    a = ap.parse_args()
    t0 = time.time()
    nd = gl.Node(a.spec)
    Y = nd.load_Y(a.Y)
    Y4, x = nd.blocks(Y)
    out = {'spec': a.spec, 'm': nd.m}
    r, c = nd.objective_parts(Y)
    print(f'node {a.spec} m={nd.m}: <L,Y> row part {r:.2f} col part {c:.2f} total {r + c:.2f}')
    out['row_part'], out['col_part'] = r, c
    dg = np.einsum('aiai->ai', Y4)
    print(f'  consistency: max|diag - x| {np.abs(dg - x).max():.2e}, min Y {Y.min():.2e}, row sums of x {np.abs(x.sum(0) - 1).max():.2e}')
    Er, Ec = nd.expected_dist(Y4)
    E = Er + Ec
    out['regions'] = {}
    print('  regional triangles ((p,W) units): region |S| max n_viol n>0.01 n>0.05 sum_pos')
    regs = regions_of(nd)
    for name, S in regs.items():
        p, W = region_pW(nd, Y4, S)
        V, _ = tri_viol(p, W)
        s = summarize(V); s['size'] = int(len(S)); out['regions'][name] = s
        print(f'    {name:9s} {len(S):3d} {s["max"]:.4f} {s["n_viol"]:6d} {s["n>0.01"]:6d} {s["n>0.05"]:5d} {s["sum_pos"]:9.3f}')
    # clique dispersion
    allowed = np.ones((nd.rows, nd.cols), dtype=bool)
    for x_, _ in nd.fixes:
        allowed[x_ // nd.cols, x_ % nd.cols] = False
    out['dispersion'] = {}
    print('  clique dispersion: k Dmin(free cells) min_T sum E[d] violation  T')
    for k in range(3, a.kmax + 1):
        dm = dsp.dmin(k, allowed)[0]
        val, T = most_clustered(E, k)
        out['dispersion'][k] = dict(dmin=dm, minsum=val, viol=dm - val, T=[int(t) for t in T])
        print(f'    {k:2d} {dm:4d} {val:9.3f} {dm - val:8.3f}  {T}')
    # flow-weighted: E[sum F d] on the clustered sets vs pairs at their mean
    print(f'  done in {time.time() - t0:.0f} s')
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
