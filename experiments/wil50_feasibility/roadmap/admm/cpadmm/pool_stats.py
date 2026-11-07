#!/usr/bin/env python3
"""Structure of a regional triangle pool (dnn_pool_encode's token) for the
box-cone question (README.md): how the pool's triangles would group into the
shared-slack cones of Charbonneau-Hu-Sotirov's Appendix A (Max-Cut triangles
grouped by a base pair: the four triangles of a triple (i, j, k) are
|Z_ik + Z_jk| <= 1 + Z_ij and |Z_ik - Z_jk| <= 1 - Z_ij in the region's
correlation coordinates Z = 1 - 2p_i - 2p_j + 4W_ij, so the triples of one
region sharing the pair (i, j) share the two slacks 1 +- Z_ij).

  python3 -I pool_stats.py POOLFILE
"""
import collections
import sys


def main(path):
    tok = open(path).read().strip()
    f = tok.split("/")
    nreg, ncut = int(f[2]), int(f[4])
    regions = [r for r in f[3].split(";") if r]
    sizes = [len(r.split(",")) for r in regions]
    cuts = [c for c in f[5].split(";") if c]
    assert len(cuts) == ncut
    types = collections.Counter()
    triples = collections.defaultdict(set)
    nzero = 0
    for c in cuts:
        x = c.split(",")
        t, s, mu = int(x[0]), int(x[1]), float(x[3])
        mem = tuple(sorted(int(v) for v in x[4:]))
        types[t] += 1
        nzero += mu == 0
        if t <= 3:
            triples[(s, mem)].add(t)
    per_triple = collections.Counter(len(v) for v in triples.values())
    # pair incidences: for each (region, pair), the distinct pool triples containing it
    pair = collections.Counter()
    for (s, (a, b, c)) in triples:
        for p in ((a, b), (a, c), (b, c)):
            pair[(s, p)] += 1
    counts = sorted(pair.values())
    # a greedy base-pair assignment: each triple to its pair with the most triples
    load = collections.Counter()
    for (s, (a, b, c)) in sorted(triples, key=lambda k: -max(pair[(k[0], p)] for p in ((k[1][0], k[1][1]), (k[1][0], k[1][2]), (k[1][1], k[1][2])))):
        best = max(((a, b), (a, c), (b, c)), key=lambda p: (pair[(s, p)], p))
        load[(s, best)] += 1
    g = sorted(load.values())
    m = 47
    print("pool %s: %d regions (sizes %s), %d cuts (%d with mu = 0), types %s" % (path, nreg, sizes, ncut, nzero, dict(types)))
    print("  distinct triples %d; inequalities per triple: %s" % (len(triples), dict(sorted(per_triple.items()))))
    print("  (region, pair) incidences: %d of %d possible pairs; triples per incidence: mean %.2f, median %d, max %d"
          % (len(pair), nreg * m * (m - 1) // 2, sum(counts) / len(counts), counts[len(counts) // 2], counts[-1]))
    print("  greedy base-pair cones: %d cones; third vertices per cone: mean %.2f, median %d, max %d; cones of size 1: %d"
          % (len(g), sum(g) / len(g), g[len(g) // 2], g[-1], sum(1 for v in g if v == 1)))


if __name__ == "__main__":
    main(sys.argv[1])
