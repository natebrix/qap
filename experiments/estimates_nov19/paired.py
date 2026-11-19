#!/usr/bin/env python3
"""Compare the arms of a sampled A/B node by node (coord/sample.py: every arm solves the same sample of
a paused run's open nodes, each to the end, with its own params).

    python3 experiments/estimates_nov19/paired.py RECORD_DIR RECORD_DIR ...       # collected records
    python3 experiments/estimates_nov19/paired.py --backend B RUN RUN ...         # from a store

A record directory is infra/phase1/runs/RUN as collect.sh leaves it: its sample.json and its node trace
(nodes.txt.gz, or export/nodes.txt.gz). The first arm is the baseline. For every sampled node, its
subtree in each arm is the traced nodes whose first d fixes (in the order they were made) are the
node's (the trace and the sample name nodes alike: latest/estimate.C, fixes_of), and its weight w is
the number of the source's open nodes at its depth it stands for (1 / (pi_dives x pi_sample)). Per arm:
the weighted CPU and nodes of the sampled subtrees, which estimate those of all the
source's open nodes at the sampled depths; against the baseline, the ratio of the weighted sums, with a
bootstrap over the sampled nodes (paired: a node is resampled in every arm at once) for its 5-95% band,
by depth and in all. A node whose subtree is unfinished in an arm (its run not finished) makes the
comparison part-way: it says so.
"""
import argparse
import collections
import gzip
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)


def read_record(path):
    sample = json.load(open(os.path.join(path, "sample.json")))
    for name in ("nodes.txt.gz", os.path.join("export", "nodes.txt.gz")):
        f = os.path.join(path, name)
        if os.path.exists(f):
            return os.path.basename(path.rstrip("/")), sample, gzip.open(f, "rt").read().splitlines(), None
    raise SystemExit(f"{path}: no nodes.txt.gz (collect.sh fetches the export)")


def read_store(b, run_id):
    from coord.common import run_key
    from coord.export import export_run
    import tempfile
    data = b.get_bytes(run_key(run_id, "sample.json"))
    if not data:
        raise SystemExit(f"{run_id}: no sample.json (not a sampled run)")
    d = tempfile.mkdtemp(prefix=f"paired-{run_id}-")
    export_run(b, run_id, d, log=lambda *_: None)
    r = b.get_run(run_id) or {}
    return run_id, json.loads(data), gzip.open(os.path.join(d, "nodes.txt.gz"), "rt").read().splitlines(), r.get("status")


def subtrees(lines, roots):
    """{root key: [cpu, nodes]} over the traced nodes under each sampled root"""
    by_depth = collections.defaultdict(set)
    for k, d in roots.items():
        by_depth[d].add(k)
    acc = {k: [0.0, 0] for k in roots}
    seen = set()
    for l in lines:
        f = l.split()
        if not f or f[0] != "node":
            continue
        d, fx = int(f[1]), f[2]
        if (d, fx) in seen:       # a node traced twice (a retried unit) counts once
            continue
        seen.add((d, fx))
        fixes = fx.split(",") if fx != "-" else []
        for rd, keys in by_depth.items():
            if d >= rd:
                k = ",".join(fixes[:rd])
                if k in keys:
                    a = acc[k]
                    a[0] += float(f[-1])
                    a[1] += 1
    return acc


def cpu(s):
    return f"{s / 3600:,.1f} CPU h" if s >= 3600 else f"{s:,.1f} CPU s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arms", nargs="+")
    ap.add_argument("--backend")
    ap.add_argument("--reps", type=int, default=2000)
    a = ap.parse_args()
    if len(a.arms) < 2:
        raise SystemExit("two or more arms")
    if a.backend:
        from coord.backend import open_backend
        b = open_backend(a.backend)
        arms = [read_store(b, r) for r in a.arms]
    else:
        arms = [read_record(p) for p in a.arms]
    base = arms[0][1]
    keyset = [x["key"] for x in base["nodes"]]
    for name, s, _, _ in arms[1:]:
        if [x["key"] for x in s["nodes"]] != keyset or s["meta"] != base["meta"]:
            raise SystemExit(f"{name}: not the same sample as {arms[0][0]} (its meta or nodes differ)")
    nodes = base["nodes"]
    roots = {x["key"]: x["depth"] for x in nodes}
    w = {x["key"]: x["w"] for x in nodes}
    acc = [subtrees(lines, roots) for _, _, lines, _ in arms]
    m = base["meta"]
    print(f"a sample of {m['src']}'s open nodes ({len(nodes)}: " +
          ", ".join(f"depth {d}: {sum(1 for x in nodes if x['depth'] == int(d))}" for d in sorted(m["plan"], key=int)) +
          f"), sized by {m['by']}, seed {m['seed']}, mix {m['mix']}")
    for (name, _, _, st), ac in zip(arms, acc):
        empty = sum(1 for k in roots if ac[k][1] == 0)
        print(f"  {name}: {sum(v[1] for v in ac.values()):,} traced nodes under the sampled ones"
              + (f"; {empty} sampled nodes not traced yet" if empty else "") + (f" (run {st})" if st and st != "finished" else ""))
    rng = random.Random(1)
    groups = [("all", list(roots))] + [(f"depth {d}", [k for k in roots if roots[k] == d]) for d in sorted(set(roots.values()))]
    print()
    print("| | " + " | ".join(n for n, *_ in arms) + " |")
    print("|---|" + "---|" * len(arms))
    for gname, ks in groups:
        tot = [sum(w[k] * ac[k][0] for k in ks) for ac in acc]
        ntot = [sum(w[k] * ac[k][1] for k in ks) for ac in acc]
        cells = []
        for i, ac in enumerate(acc):
            if i == 0:
                cells.append(f"{cpu(tot[0])}, {ntot[0]:,.0f} nodes")
                continue
            boots = []
            for _ in range(a.reps):
                draw = [ks[rng.randrange(len(ks))] for _ in ks]
                b0 = sum(w[k] * acc[0][k][0] for k in draw)
                bi = sum(w[k] * ac[k][0] for k in draw)
                if b0 > 0:
                    boots.append(bi / b0)
            boots.sort()
            q = lambda f: boots[min(int(f * len(boots)), len(boots) - 1)] if boots else float("nan")
            r = tot[i] / tot[0] if tot[0] else float("nan")
            rn = ntot[i] / ntot[0] if ntot[0] else float("nan")
            cells.append(f"{cpu(tot[i])} ({r:.2f}x, 5-95% {q(.05):.2f}-{q(.95):.2f}), nodes {rn:.2f}x")
        print(f"| {gname} ({len(ks)}) | " + " | ".join(cells) + " |")
    print()
    print("(weighted: each sampled node stands for its share of the source's open nodes at its depth; the "
          "CPU of two arms is comparable only when they ran side by side on the same instance types)")


if __name__ == "__main__":
    main()
