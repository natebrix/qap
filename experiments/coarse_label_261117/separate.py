#!/usr/bin/env python3
"""Separate a small fixed pool of local moment inequalities for one saved QAP node."""

import argparse
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform


def data_at_node(instance, fixes, side, grid):
    data = np.fromfile(instance, sep=" ")
    n = int(data[0])
    matrices = data[1:1 + 2 * n * n].reshape(2, n, n)
    assignments = [] if fixes == "-" else [tuple(map(int, p.split(":"))) for p in fixes.split(",")]
    free = [np.array(sorted(set(range(n)) - {p[s] for p in assignments})) for s in (0, 1)]
    region_matrix = matrices[side]
    if grid:
        rows, cols = map(int, grid.split(":"))
        if rows * cols != n:
            raise ValueError("grid must cover the original index set")
        idx = np.arange(n)
        labels = 2 * (idx // cols >= rows // 2) + (idx % cols >= cols // 2)
    else:
        d = (region_matrix + region_matrix.T) / 2
        labels = fcluster(linkage(squareform(d, checks=False), "average"), 4, "maxclust") - 1
    flow = matrices[1 - side][np.ix_(free[1 - side], free[1 - side])]
    flow = np.abs(flow + flow.T)
    return free[1 - side], labels[free[side]], flow, labels


def local_vertices(t, capacity):
    k = len(capacity)
    pairs = list(itertools.combinations(range(t), 2))
    labels = np.array([s for s in itertools.product(range(k), repeat=t)
                       if np.all(np.bincount(s, minlength=k) <= capacity)])
    v = np.zeros((len(labels), t * k + len(pairs) * k * k), dtype=np.int64)
    for i in range(t):
        v[np.arange(len(labels)), i * k + labels[:, i]] = 1
    for q, (i, j) in enumerate(pairs):
        v[np.arange(len(labels)), t * k + q * k * k + k * labels[:, i] + labels[:, j]] = 1
    return v, pairs


def subsets(flow, t, samples, rng):
    m = len(flow)
    edges = sorted(((int(i), int(j)) for i, j in zip(*np.triu_indices(m, 1))), key=lambda e: -flow[e])
    sets = set()
    for edge in edges[:samples]:
        members = list(edge)
        while len(members) < t:
            scores = flow[members].sum(0)
            scores[members] = -np.inf
            members.append(int(scores.argmax()))
        sets.add(tuple(sorted(members)))
    for _ in range(samples):
        sets.add(tuple(sorted(rng.choice(m, t, replace=False).tolist())))
    return sorted(sets)


def aggregate(y, labels, side, k):
    m = len(labels)
    ri, ci = np.indices((m, m))
    members, places = (ri, ci) if side == 1 else (ci, ri)
    groups = (members * k + labels[places]).ravel()
    t = csr_matrix((np.ones(m * m), (groups, np.arange(m * m))), shape=(m * k, m * m))
    p = (t @ np.diag(y)[1:]).reshape(m, k)
    w = (t @ (t @ y[1:, 1:]).T).T
    return p, (0.5 * (w + w.T)).reshape(m, k, m, k)


def separate(y, labels, side, flow, sizes, samples, seed):
    k = int(labels.max()) + 1
    cap = np.bincount(labels, minlength=k)
    p, w = aggregate(y, labels, side, k)
    rng = np.random.default_rng(seed)
    candidates = []
    per_size = []
    for t in sizes:
        if t > len(labels):
            continue
        v, pairs = local_vertices(t, cap)
        constraints = csr_matrix(np.column_stack((v, -np.ones(len(v)))))
        origin = v.mean(axis=0)
        _, singular, vt = np.linalg.svd(v - origin, full_matrices=False)
        basis = vt[singular > 1e-9].T
        norms = np.r_[np.tile(cap, t),
                      np.tile((np.outer(cap, cap) - np.diag(cap)).ravel() / 2, len(pairs))]
        members_sets = subsets(flow, t, samples, rng)
        start_size = len(candidates)
        for members in members_sets:
            raw = np.r_[p[list(members)].ravel(),
                        np.concatenate([w[members[i], :, members[j], :].ravel() for i, j in pairs])]
            repaired = origin + basis @ (basis.T @ (raw - origin))
            result = linprog(np.r_[-repaired, 1], A_ub=constraints, b_ub=np.zeros(len(v)),
                             bounds=[(-1, 1)] * v.shape[1] + [(None, None)], method="highs-ds")
            if not result.success:
                raise RuntimeError(result.message)
            if result.fun >= -1e-8:
                continue
            # Affine-equivalent directions differ on an inexact ADMM iterate.
            # Choose the one orthogonal to the affine equations before rounding.
            direction = basis @ (basis.T @ result.x[:-1])
            coef = np.rint(direction * (1024 / np.max(np.abs(direction)))).astype(np.int64)
            rhs = int(np.max(v @ coef))
            norm = float(np.sqrt(norms @ (coef.astype(float) ** 2)))
            if norm == 0:
                continue
            repaired_v = float(repaired @ coef - rhs) / norm
            raw_v = float(raw @ coef - rhs) / norm
            if min(repaired_v, raw_v) <= 1e-7:
                continue
            candidates.append({"members": members, "coef": coef.tolist(), "rhs": rhs,
                               "lift_norm": norm / 1024, "raw_violation_per_norm": raw_v,
                               "repaired_violation_per_norm": repaired_v,
                               "rank": min(raw_v, repaired_v)})
        per_size.append({"t": t, "configurations": len(v), "subsets": len(members_sets),
                         "qualified": len(candidates) - start_size})
    candidates.sort(key=lambda c: -c["rank"])
    return candidates, per_size


def write_pool(path, side, labels, cuts):
    k = int(labels.max()) + 1
    with path.open("w") as out:
        print("QAP_COARSE_POOL 1", side, len(labels), k, len(cuts), file=out)
        print(*labels, file=out)
        for cut in cuts:
            print(len(cut["members"]), *cut["members"], *cut["coef"], file=out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instance", type=Path)
    parser.add_argument("y", type=Path)
    parser.add_argument("pool", type=Path)
    parser.add_argument("--fixes", default="-")
    parser.add_argument("--side", type=int, choices=(0, 1), required=True)
    parser.add_argument("--grid", help="quadrant partition R:C, otherwise four average-linkage clusters")
    parser.add_argument("--sizes", type=int, choices=range(2, 7), nargs="+", default=[4, 5])
    parser.add_argument("--samples", type=int, default=35, help="number of flow seeds and number of random subsets per size")
    parser.add_argument("--cut-count", type=int, choices=range(1, 51), default=40)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    start_cpu, start_wall = time.process_time(), time.monotonic()
    members, labels, flow, original_labels = data_at_node(args.instance, args.fixes, args.side, args.grid)
    n1 = len(labels) ** 2 + 1
    y = np.memmap(args.y, np.float64, mode="r", shape=(n1, n1))
    candidates, sizes = separate(y, labels, args.side, flow, args.sizes, args.samples, args.seed)
    cuts = candidates[:args.cut_count]
    for c in cuts:
        c["original_members"] = members[list(c["members"])].tolist()
    write_pool(args.pool, args.side, labels, cuts)
    report = {"instance": str(args.instance), "fixes": args.fixes, "side": args.side, "grid": args.grid,
              "seed": args.seed, "samples": args.samples, "sizes": sizes, "n_cuts": len(cuts),
              "capacity": np.bincount(labels).tolist(), "original_labels": original_labels.tolist(),
              "cpu_s": time.process_time() - start_cpu, "wall_s": time.monotonic() - start_wall}
    print(json.dumps(report, indent=2))
    args.pool.with_suffix(".json").write_text(json.dumps({**report, "cuts": cuts}, indent=2) + "\n")


if __name__ == "__main__":
    main()
