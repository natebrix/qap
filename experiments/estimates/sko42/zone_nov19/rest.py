#!/usr/bin/env python3
"""What is left of a stopped run's tree, from its node trace. An open node is a child a branched node made
that is not traced; each is imputed the mean CPU and nodes of the closed subtrees rooted at its depth, or of
its closed siblings' (three or more; else its depth's). Closed subtrees are the small ones, so both run low:
backtest.py measures by how much. In all and by depth-1 subtree.

    python3 experiments/estimates/sko42/zone_nov19/rest.py [TRACE ...]
"""
import collections
import statistics as st
import sys

from common import ARMS, load


def impute(t):
    nodes, ch, cpu, size, closed = t["nodes"], t["ch"], t["cpu"], t["size"], t["closed"]
    bydep = collections.defaultdict(list)
    for k, n in nodes.items():
        if closed[k]:
            bydep[n["depth"]].append((cpu[k], size[k]))
    mean = {d: (st.fmean(c for c, _ in v), st.fmean(s for _, s in v)) for d, v in bydep.items()}
    d1 = collections.defaultdict(lambda: dict(done=0.0, open=0, dep=0.0, sib=0.0, sib_nodes=0.0, closed=False, maxd=0))
    opened = collections.Counter()
    for k, n in nodes.items():
        g = d1[k.split(",")[0] if k else ""]
        g["done"] += n["cpu"]
        g["maxd"] = max(g["maxd"], n["depth"])
        if n["out"] != "B":
            continue
        miss = n["kids"] - len(ch.get(k, []))
        if miss <= 0:
            continue
        d = n["depth"] + 1
        md, ms = mean.get(d, (0.0, 1.0))
        cs = [(cpu[x], size[x]) for x in ch.get(k, []) if closed[x]]
        sc, ss = (st.fmean(c for c, _ in cs), st.fmean(s for _, s in cs)) if len(cs) >= 3 else (md, ms)
        g["open"] += miss
        g["dep"] += miss * md
        g["sib"] += miss * sc
        g["sib_nodes"] += miss * ss
        opened[d] += miss
    for r in ch[""]:
        d1[r]["closed"] = closed[r]
    return d1, opened, mean


def main():
    for path in sys.argv[1:] or [ARMS["z10"], ARMS["z25"]]:
        t = load(path)
        d1, opened, mean = impute(t)
        done = t["cpu"][""] / 3600
        dep = sum(g["dep"] for g in d1.values()) / 3600
        sib = sum(g["sib"] for g in d1.values()) / 3600
        sn = sum(g["sib_nodes"] for g in d1.values())
        print(f"{path.split('/')[-2]}: done {done:,.0f} CPU h, {len(t['nodes']):,} nodes; open nodes {sum(opened.values()):,}")
        print("  open by depth: " + ", ".join(f"{d}: {opened[d]:,}" for d in sorted(opened)))
        print(f"  rest at the depth means {dep:,.0f} CPU h; at the siblings' {sib:,.0f} CPU h, {sn:,.0f} nodes")
        print(f"  tree {done + dep:,.0f} to {done + sib:,.0f} CPU h, {len(t['nodes']) + sn:,.0f} nodes (siblings')")
        print("  | depth-1 node | done CPU h | open nodes | rest CPU h (depth / siblings) | deepest |")
        print("  |---|---|---|---|---|")
        for k, g in sorted(((k, g) for k, g in d1.items() if k), key=lambda x: -x[1]["done"]):
            print(f"  | {k} | {g['done'] / 3600:.0f} | {'closed' if g['closed'] else format(g['open'], ',')} | "
                  f"{g['dep'] / 3600:.0f} / {g['sib'] / 3600:.0f} | {g['maxd']} |")


if __name__ == "__main__":
    main()
