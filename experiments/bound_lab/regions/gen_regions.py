#!/usr/bin/env python3
"""Region files for `dnncuts -regions FILE -side S` (bound_lab/regions).

  python3 -I gen_regions.py INSTANCE.dat FAMILY[,FAMILY...] [-o OUT] [-seed N] [-stats]

One region per line: original 0-based indices of one side. In tai*b.dat the first matrix
is the symmetric distance (side 0, locations) and the second the asymmetric flow (side 1,
facilities; the solver symmetrises it).

Location families (side 0), from the first matrix D:
  hier        every cluster of the average-linkage dendrogram of D (sizes 2..n-2)
  top:K       the K top clusters of that dendrogram and all unions of them (K <= 8)
  line:K      half-plane prefixes: classical MDS to 2-D, K directions theta = pi j / K
              (j = 0 is the first principal axis), prefixes of the locations sorted by
              projection, one per gap between distinct projections (n-1 per direction)
  linegap:K:G as line:K, but only prefixes whose gap to the next projection is at least
              G times the diameter (cuts between clusters only)
  ball:K      balls {b : d(a,b) <= r}, r = the distance from a to its k-th neighbour,
              k = K values spread over 1..n-2
  random:R    R random subsets (validity tests)
Facility families (side 1), from the symmetrised second matrix F:
  fhier       dendrogram clusters of the flow graph (average linkage on 1/(1+flow))
  fhub:K      prefixes of the facilities sorted by total flow (the K top hubs: K regions)
  fspec:K     prefixes along the first K nontrivial eigenvectors of the normalised
              Laplacian of the flow graph
"""
import sys
import numpy as np


def read_qaplib(path):
    t = open(path).read().split()
    n = int(t[0])
    v = np.array(t[1:1 + 2 * n * n], dtype=float)
    return n, v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)


def mds(D, dim=2):
    n = len(D)
    J = np.eye(n) - 1.0 / n
    G = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(G)
    idx = np.argsort(-w)
    w, V = w[idx], V[:, idx]
    X = V[:, :dim] * np.sqrt(np.maximum(w[:dim], 0))
    D2 = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    stress = np.sqrt(((D - D2) ** 2).sum() / (D ** 2).sum())
    return X, w, stress


def dendro_clusters(D, method='average'):
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    n = len(D)
    Z = linkage(squareform(D, checks=False), method)
    members = {i: [i] for i in range(n)}
    out = []
    for k, (a, b, h, sz) in enumerate(Z):
        members[n + k] = members[int(a)] + members[int(b)]
        out.append((h, sorted(members[n + k])))
    return Z, out


def top_clusters(D, K):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Z = linkage(squareform(D, checks=False), 'average')
    lab = fcluster(Z, K, 'maxclust')
    return [list(np.where(lab == c)[0]) for c in range(1, lab.max() + 1)]


def line_regions(D, K, mingap=0.0):
    X, w, stress = mds(D, 2)
    n = len(D)
    diam = D.max()
    out = []
    for j in range(K):
        th = np.pi * j / K
        pr = X[:, 0] * np.cos(th) + X[:, 1] * np.sin(th)
        order = np.argsort(pr, kind='stable')
        for t in range(1, n):
            gap = pr[order[t]] - pr[order[t - 1]]
            if gap <= 1e-9 * diam or gap < mingap * diam:
                continue
            out.append(sorted(order[:t].tolist()))
    return out


def canon(regions, n):
    """drop empty/full sets, duplicates and complements (the cuts of S and its complement coincide)"""
    seen, out = set(), []
    for r in regions:
        r = tuple(sorted(set(int(x) for x in r)))
        c = tuple(sorted(set(range(n)) - set(r)))
        if not r or not c or r in seen or c in seen:
            continue
        seen.add(r)
        out.append(list(r))
    return out


def family(spec, n, D, F, rng):
    a = spec.split(':')
    name = a[0]
    Fs = F + F.T
    if name == 'hier':
        return [c for h, c in dendro_clusters(D)[1] if 2 <= len(c) <= n - 2]
    if name == 'top':
        K = int(a[1])
        cl = top_clusters(D, K)
        out = []
        for mask in range(1, 1 << len(cl)):
            s = []
            for k in range(len(cl)):
                if mask >> k & 1:
                    s += cl[k]
            out.append(s)
        return out
    if name == 'line':
        return line_regions(D, int(a[1]))
    if name == 'linegap':
        return line_regions(D, int(a[1]), float(a[2]))
    if name == 'ball':
        K = int(a[1])
        ks = sorted(set(np.linspace(1, n - 2, K).round().astype(int).tolist()))
        out = []
        for c in range(n):
            order = np.argsort(D[c], kind='stable')
            for k in ks:
                r = D[c, order[k]]
                out.append(np.where(D[c] <= r)[0].tolist())
        return out
    if name == 'random':
        R = int(a[1])
        out = []
        for _ in range(R):
            k = rng.integers(1, n)
            out.append(rng.choice(n, size=k, replace=False).tolist())
        return out
    if name == 'fhier':
        Dd = 1.0 / (1.0 + Fs)
        np.fill_diagonal(Dd, 0)
        return [c for h, c in dendro_clusters(Dd)[1] if 2 <= len(c) <= n - 2]
    if name == 'fhub':
        K = int(a[1])
        order = np.argsort(-Fs.sum(1), kind='stable')
        return [order[:t].tolist() for t in range(1, min(K, n - 1) + 1)]
    if name == 'fspec':
        K = int(a[1])
        W = Fs.copy()
        np.fill_diagonal(W, 0)
        d = W.sum(1) + 1e-9
        L = np.eye(n) - W / np.sqrt(d[:, None] * d[None])
        w, V = np.linalg.eigh(L)
        out = []
        for k in range(1, K + 1):
            v = V[:, k] / np.sqrt(d)
            order = np.argsort(v, kind='stable')
            out += [order[:t].tolist() for t in range(1, n)]
        return out
    raise SystemExit('unknown family ' + spec)


def stats(n, D, F):
    Fs = F + F.T
    X, w, stress = mds(D, 2)
    print(f'n {n}  MDS eigenvalues {w[:4].round(0)} min {w[-1]:.0f}  2-D stress {stress:.4f}  '
          f'1-D share {w[0] / w[:2].sum():.3f}')
    Z, cl = dendro_clusters(D)
    h = Z[:, 2]
    print('  last merge heights', h[-8:].round(0))
    # the natural number of clusters: the largest ratio of consecutive merge heights
    ratios = h[1:] / np.maximum(h[:-1], 1e-9)
    k = int(np.argmax(ratios[len(h) // 2:]) + len(h) // 2)
    K = n - 1 - k
    cls = top_clusters(D, K)
    lab = np.empty(n, int)
    for c, s in enumerate(cls):
        lab[s] = c
    same = lab[:, None] == lab[None]
    off = ~np.eye(n, dtype=bool)
    print(f'  natural clusters {K} sizes {[len(s) for s in cls]}  (merge height jumps {h[k]:.0f} -> {h[k + 1]:.0f})')
    print(f'  mean distance within {D[same & off].mean():.0f}  between {D[~same].mean():.0f}  '
          f'share of sum(D) between clusters {D[~same].sum() / D.sum():.3f}')
    nz = Fs[np.triu_indices(n, 1)]
    nzs = np.sort(nz[nz > 0])[::-1]
    print(f'  flows: {len(nzs)} of {n * (n - 1) // 2} pairs nonzero, top 10% of them carry '
          f'{nzs[:max(1, len(nzs) // 10)].sum() / nzs.sum():.2f} of the total')
    return cls


if __name__ == '__main__':
    args = sys.argv[1:]
    out, seed, want_stats, pos = None, 1, False, []
    i = 0
    while i < len(args):
        if args[i] == '-o':
            out = args[i + 1]; i += 2
        elif args[i] == '-seed':
            seed = int(args[i + 1]); i += 2
        elif args[i] == '-stats':
            want_stats = True; i += 1
        else:
            pos.append(args[i]); i += 1
    n, D, F = read_qaplib(pos[0])
    if want_stats:
        stats(n, D, F)
    if len(pos) > 1:
        rng = np.random.default_rng(seed)
        regs = []
        for spec in pos[1].split(','):
            regs += family(spec, n, D, F, rng)
        regs = canon(regs, n)
        f = open(out, 'w') if out else sys.stdout
        f.write(f'# {pos[0].split("/")[-1]} {pos[1]}: {len(regs)} regions\n')
        for r in regs:
            f.write(' '.join(map(str, r)) + '\n')
        if out:
            f.close()
            print(f'{out}: {len(regs)} regions', file=sys.stderr)
