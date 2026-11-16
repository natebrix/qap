"""Each phase of a cut arm against the base arm's cost under the same node.
  python3 matched.py ARM [--list]
A node is named by its set of fixes. For every phase of ARM whose node the base arm also visited:
  base  = base iterations in the node's subtree, the node's own first solve excluded
  arm   = the phase's iterations + ARM's iterations below the node
(iterations weighted: a phase iteration counts W first-solve iterations, W = CPU ratio measured in ARM at that depth).
"""
import sys, collections
from summarize import load

def key(e): return frozenset(e["fixes"].split())

def subtree(dnn):
    """iterations (first, phase) strictly below each node, by fix set"""
    idx = {key(e): e for e in dnn}
    below = collections.defaultdict(lambda: [0, 0, 0])
    for e in dnn:
        k = key(e)
        # walk up: every ancestor fix-subset that is a node. Ancestors are found by size, so test the known nodes of smaller depth
        for a in anc.get(e["depth"], ()):  # filled below
            pass
    return idx, below

def main():
    arm = sys.argv[1]
    _, b = load("base"); _, c = load(arm)
    def below(dnn):
        bydepth = collections.defaultdict(dict)
        for e in dnn: bydepth[e["depth"]][key(e)] = e
        out = collections.defaultdict(lambda: [0, 0, 0])   # first it, phase it, nodes
        for e in dnn:
            k = key(e)
            for d in range(1, e["depth"]):
                for a in bydepth[d]:
                    if a < k:
                        o = out[a]; o[0] += e["iters"]; o[1] += e["cit"]; o[2] += 1
        return bydepth, out
    bd, bb = below(b); cd, cb = below(c)
    # cost of a phase iteration against a first-solve iteration, per depth, in ARM
    W = {}
    for d in cd:
        f = [(e["cpu"], e["iters"]) for e in cd[d].values()]; p = [(e["ccpu"], e["cit"]) for e in cd[d].values() if e["cut"]]
        if p and sum(x[1] for x in p):
            W[d] = (sum(x[0] for x in p) / sum(x[1] for x in p)) / (sum(x[0] for x in f) / sum(x[1] for x in f))
    print("phase iteration cost / first-solve iteration cost by depth:", {d: round(w, 2) for d, w in sorted(W.items())})
    rows = []
    for e in c:
        if not e["cut"]: continue
        k = key(e)
        if k not in bd[e["depth"]]: continue
        w = W[e["depth"]]
        base = bb[k][0]; armc = e["cit"] * w + cb[k][0] + cb[k][1] * w
        rows.append((e["cutoff"] - e["lb"], e["depth"], e["kept"] > 0, e["cstatus"], e["cit"], e["cutoff"] - e["cbest"], bb[k][2], base, cb[k][2], armc, e["iters"], bd[e["depth"]][k]["iters"]))
    rows.sort()
    if "--list" in sys.argv:
        print("short depth pooled status phase_it left | base: nodes below, it below | arm: nodes below, weighted it (phase + below) | first-solve it arm/base")
        for r in rows:
            print("%6.1f d%d %s %-9s %5d %6.1f | %4d %6d | %4d %7.0f | %d/%d %s" % (r[0], r[1], "pool" if r[2] else "    ", r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], "WIN" if r[9] < r[7] else "loss"))
    print("\n| first-solve shortfall | phases | fathomed | base iterations below these nodes | arm: phase + below (weighted) | arm / base | phases that paid |\n|---|---|---|---|---|---|---|")
    bins = [0, 5, 10, 15, 20, 30, 40, 61.3, 100, 184, 400]
    for lo, hi in zip(bins, bins[1:]):
        r = [x for x in rows if lo <= x[0] < hi]
        if not r: continue
        sb = sum(x[7] for x in r); sa = sum(x[9] for x in r)
        print("| %g-%g | %d | %d | %s | %s | %.2f | %d |" % (lo, hi, len(r), sum(x[3] == "fathomed" for x in r), f"{sb:,}", f"{sa:,.0f}", sa / max(1, sb), sum(x[9] < x[7] for x in r)))
    sb = sum(x[7] for x in rows); sa = sum(x[9] for x in rows)
    print("| all matched phases | %d | %d | %s | %s | %.2f | %d |" % (len(rows), sum(x[3] == "fathomed" for x in rows), f"{sb:,}", f"{sa:,.0f}", sa / max(1, sb), sum(x[9] < x[7] for x in rows)))
    print("oracle (phase only where it paid): saves %s weighted iterations of base's %s" % (f"{sum(x[7]-x[9] for x in rows if x[9] < x[7]):,.0f}", f"{sum(e['iters'] for e in b):,}"))

main()
