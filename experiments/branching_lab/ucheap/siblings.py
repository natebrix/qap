#!/usr/bin/env python3
"""ucheap/siblings.py: within a node's children, does the child's own QPB bound say which subtrees are large?

  python3 siblings.py TRACE [TRACE ...]      node traces (nodes.txt.gz) of whole trees

For every branched node, its children as the search bounded them: each child's QPB bound and first DNN bound (the
trace's columns) against the size of the subtree under it (trace lines). Spearman within each sibling group with at
least 5 children carrying both bounds and not all the same size, averaged by the parent's depth. This is the
information a QPB strong branching (one QPB solve a child) would have, against the DNN's.
"""
import gzip, sys, collections
sys.path.insert(0, "/home/user/wt-port/experiments/line_signals")
from signals import spearman

def load(path):
    nodes = {}
    for ln in (gzip.open(path, "rt") if path.endswith(".gz") else open(path)):
        f = ln.split()
        if not f or f[0] != "node": continue
        key = frozenset(f[2].split(",")) if f[2] != "-" else frozenset()
        nodes[key] = (int(f[1]), f[3], f[4], f[5], float(f[6]), f[7])
    return nodes

def main():
    by_depth = collections.defaultdict(lambda: [[], [], [], 0, 0])   # rho qpb, rho dnn, worst-child hits (q, d), groups
    for path in sys.argv[1:]:
        nodes = load(path)
        size = collections.Counter()
        for key in nodes:                       # every ancestor's subtree gets this node
            pass
        # subtree sizes: count each node into all of its ancestors that are in the trace, by walking prefixes is
        # not possible on sets, so bucket children by parent first
        kids = collections.defaultdict(list)
        keys_by_len = collections.defaultdict(list)
        for key in nodes: keys_by_len[len(key)].append(key)
        sub = {key: 1 for key in nodes}
        parent = {}
        for L in sorted(keys_by_len, reverse=True):
            if L == 0: continue
            up = set(keys_by_len.get(L - 1, []))
            for key in keys_by_len[L]:
                p = next((key - {x} for x in key if key - {x} in up), None)
                if p is None: continue
                parent[key] = p
        for L in sorted(keys_by_len, reverse=True):
            for key in keys_by_len[L]:
                if key in parent:
                    sub[parent[key]] += sub[key]; kids[parent[key]].append(key)
        for p, ch in kids.items():
            rows = []
            for k in ch:
                d, q, dn, cu, cut, out = nodes[k]
                if q == "-" or dn == "-": continue
                rows.append((cut - float(q), cut - float(dn), sub[k]))
            if len(rows) < 5 or len({r[2] for r in rows}) < 2: continue
            b = by_depth[len(p)]
            b[0].append(spearman([r[0] for r in rows], [r[2] for r in rows]))
            b[1].append(spearman([r[1] for r in rows], [r[2] for r in rows]))
            big = max(rows, key=lambda r: r[2])
            b[3] += max(rows, key=lambda r: r[0]) is big
            b[4] += max(rows, key=lambda r: r[1]) is big
            b[2].append(len(rows))
    print("| parent's depth | sibling groups | children a group | Spearman with the subtree: QPB shortfall | DNN shortfall | largest subtree is the largest QPB shortfall | the largest DNN shortfall |")
    print("|---|---|---|---|---|---|---|")
    for d in sorted(by_depth):
        q, dn, n, hq, hd = by_depth[d]
        print(f"| {d} | {len(q)} | {sum(n) / len(n):.1f} | {sum(q) / len(q):+.2f} | {sum(dn) / len(dn):+.2f} | {hq}/{len(q)} | {hd}/{len(q)} |")

if __name__ == "__main__":
    main()
