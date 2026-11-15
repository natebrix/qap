"""Whole-tree estimate of tho40-b from the laptop's Knuth dives (Nov 15): one bnbdrv -e
process per depth-1 node 26:c (s<c>/), tho40-b's recipe and incumbent, start level 3
(the node and its depth-2 children bounded exactly, dives from depth 3), q = 10.

Each dive's line in s<c>/dives.txt: its last depth, its estimate of the stratum's nodes
(the exact top included), its estimate of the stratum's CPU in laptop seconds. The tree
is the root plus the 12 strata. Laptop CPU is put in c7a seconds by the same nodes in
tho40-b (the 12 depth-1 nodes and the depth-2 nodes it has merged, from the live tree's
root shard). The bands are a bootstrap over each stratum's dives, which a heavy tail
makes too narrow: a stratum whose sum one dive carries reads low more often than high.

usage: python3 strata.py [ROOT_SHARD.json]   (default: fetched from the site)
"""
import json
import os
import random
import re
import statistics
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
STRATA = (0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else None
    sh = json.load(open(src)) if src else json.loads(urllib.request.urlopen(
        "https://nathanbrixius.com/projects/qap/tree/tho40-b/data/0.json", timeout=60).read())
    C = {c: i for i, c in enumerate(sh["columns"])}
    rows = sh["rows"]
    root = [r for r in rows if r[C["depth"]] == 0][0]
    d1 = {r[C["fix"]]: r for r in rows if r[C["depth"]] == 1}
    d2 = [r for r in rows if r[C["depth"]] == 2]
    rng = random.Random(11)
    lap = run = 0.0
    res = []
    for c in STRATA:
        out = open(f"{HERE}/s{c}/est.out").read()
        m = re.search(r"level    nodes\t%fathom\t%c.elim\t    time\n 1\s+1\s+\S+\s+\S+\s+(\S+)\n 2\s+(\d+)\s+\S+\s+\S+\s+(\S+)", out)
        t1, n2, t2 = float(m.group(1)), int(m.group(2)), float(m.group(3))
        p = d1[f"26:{c}"]
        kids = [r for r in d2 if r[C["parent"]] == p[C["id"]]]
        lap += t1 + (t2 / n2 * len(kids) if kids else 0)
        run += p[C["cpu"]] + sum(r[C["cpu"]] for r in kids)
        dv = [tuple(map(float, l.split())) for l in open(f"{HERE}/s{c}/dives.txt") if l.strip()]
        N, T = [x[1] for x in dv], [x[2] for x in dv]
        res.append(dict(c=c, N=N, T=T, merged=p[C["size"]]))
    f = run / lap
    boot = lambda v: statistics.fmean(rng.choices(v, k=len(v)))
    print(f"c7a CPU per laptop CPU, same nodes: {f:.2f}")
    print(f"{'node':>6} {'dives':>5} {'nodes: mean (80% band)':>34} {'one dive':>8} {'merged':>9} {'CPU, c7a days':>13}")
    for r in res:
        bs = sorted(boot(r["N"]) for _ in range(2000))
        mean = statistics.fmean(r["N"])
        print(f"{'26:' + str(r['c']):>6} {len(r['N']):>5} {mean:>13,.0f} ({bs[200]:,.0f}-{bs[1800]:,.0f}) "
              f"{max(r['N']) / sum(r['N']):>8.2f} {r['merged']:>9,} {f * statistics.fmean(r['T']) / 86400:>13.1f}"
              + ("  below merged" if mean < r["merged"] else ""))
    B = sorted((1 + sum(boot(r["N"]) for r in res), root[C["cpu"]] + f * sum(boot(r["T"]) for r in res))
               for _ in range(2000))
    bn, bt = sorted(b[0] for b in B), sorted(b[1] for b in B)
    n_tot = 1 + sum(statistics.fmean(r["N"]) for r in res)
    t_tot = root[C["cpu"]] + f * sum(statistics.fmean(r["T"]) for r in res)
    top = max(max(r["N"]) / len(r["N"]) for r in res) / n_tot
    print(f"tree: {n_tot:,.0f} nodes (80% {bn[200]:,.0f}-{bn[1800]:,.0f}), {t_tot / 86400:,.0f} c7a CPU days "
          f"(80% {bt[200] / 86400:,.0f}-{bt[1800] / 86400:,.0f}); the largest single dive carries {top:.0%} of the nodes")
    print(f"tho40-b so far: {root[C['size']]:,} nodes, {root[C['subCpu']] / 86400:,.1f} CPU days")


if __name__ == "__main__":
    main()
