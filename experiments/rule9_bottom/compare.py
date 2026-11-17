#!/usr/bin/env python3
"""Two or more whole-run node traces (QAP_NODE_TRACE, as qapctl export or collect.sh leave them,
nodes.txt or nodes.txt.gz) side by side, by depth: nodes, CPU, how many branched, how often the
children of a depth's branchings branch in turn, and how many distinct lines that depth's branchings
used (rule 2 picks one location per depth on the tai-b instances: README.md).

  python3 experiments/rule9_bottom/compare.py base=infra/phase1/runs/tai50b-261116/nodes.txt.gz \
      r9=infra/phase1/runs/tai50b-r9/nodes.txt.gz

Run from the repository root (it imports coord.tree, which rebuilds the tree as the website does).
"""
import collections
import gzip
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from coord.tree import LiveTree  # noqa: E402


def load(path):
    t = LiveTree()
    with (gzip.open if path.endswith(".gz") else open)(path, "rb") as f:
        t.ingest(f.read())
    return t


def by_depth(t):
    """depth -> nodes, cpu, branched, children of branchings, of them branched, lines used (Counter)"""
    D = collections.defaultdict(lambda: {"n": 0, "cpu": 0.0, "b": 0, "kids": 0, "kidb": 0,
                                         "lines": collections.Counter()})
    for i in range(len(t)):
        d = D[t.dep[i]]
        d["n"] += 1
        d["cpu"] += t.cpu[i]
        if t.out[i] == ord("B"):
            d["b"] += 1
            ln = t.line[i]
            d["lines"][("r%d" % (ln - 1)) if ln > 0 else ("c%d" % (-ln - 1)) if ln < 0 else "-"] += 1
            for c in t.children(i):
                d["kids"] += 1
                d["kidb"] += t.out[c] == ord("B")
    return D


def main(args):
    if len(args) < 1 or any("=" not in a for a in args):
        sys.exit(__doc__)
    runs = [(a.split("=", 1)[0], by_depth(load(a.split("=", 1)[1]))) for a in args]
    for name, D in runs:
        n = sum(d["n"] for d in D.values()); cpu = sum(d["cpu"] for d in D.values())
        deep_n = sum(d["n"] for k, d in D.items() if k >= 7); deep_c = sum(d["cpu"] for k, d in D.items() if k >= 7)
        print(f"{name}: {n:,} nodes, {cpu / 3600:,.1f} core-hours; depth 7 and below {deep_n:,} nodes "
              f"({100 * deep_n / n:.0f}%), {deep_c / 3600:,.1f} core-hours ({100 * deep_c / cpu:.0f}%)")
    print()
    head = "depth"
    for name, _ in runs:
        head += f" | {name}: nodes  CPU h  branched  kids-branch%  lines (top share)"
    print(head)
    for k in range(max(max(D) for _, D in runs) + 1):
        row = f"{k:5d}"
        for _, D in runs:
            d = D.get(k)
            if not d:
                row += " |" + " " * 40
                continue
            top = d["lines"].most_common(1)
            share = f"{len(d['lines'])} ({top[0][0]} {100 * top[0][1] / d['b']:.0f}%)" if top else "-"
            kb = f"{100 * d['kidb'] / d['kids']:5.1f}" if d["kids"] else "    -"
            row += f" | {d['n']:6,d} {d['cpu'] / 3600:6.1f} {d['b']:5d} {kb}  {share}"
        print(row)


if __name__ == "__main__":
    main(sys.argv[1:])
