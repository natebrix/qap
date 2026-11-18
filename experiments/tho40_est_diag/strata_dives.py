"""Where tho40-est's estimate comes from: strata by depth, dives by length, the heavy tail.

    python3 -I strata_dives.py
"""
import collections, math, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from load import load

units, strata, dives, samples, est = load()
H = 3600.0

# 1. the sample: open nodes per unit, kept, by the units' depth
open_n = sum(s["open"] for s in samples)
kept = sum(s["kept"] for s in samples)
print(f"open nodes in merged units {open_n:,}, kept (strata) {kept:,} ({kept / open_n:.1%})")

# 2. strata by depth
print("\n| stratum depth | strata | sum 1/pi (open nodes they stand for) | pi median | branched (kids) | dives | dives/branched stratum |")
print("|---|---|---|---|---|---|---|")
byd = collections.defaultdict(list)
for k, s in strata.items():
    byd[s["depth"]].append(s)
for d in sorted(byd):
    ss = byd[d]
    br = [s for s in ss if s["survivors"] > 0]
    nd = sum(len(dives.get(s["key"], [])) for s in ss)
    print(f"| {d} | {len(ss):,} | {sum(1 / s['pi'] for s in ss):,.0f} | {statistics.median(s['pi'] for s in ss):.3f} | "
          f"{len(br):,} | {nd:,} | {nd / max(len(br), 1):.1f} |")

# 3. dives: length (deepest - stratum depth), by stratum depth
print("\n| stratum depth | dives | ended at depth+1 | +2 | +3 | +4 | deepest |")
print("|---|---|---|---|---|---|---|")
for d in sorted(byd):
    ls = collections.Counter()
    deep = 0
    for s in byd[d]:
        for x in dives.get(s["key"], []):
            ls[x["deepest"] - d] += 1
            deep = max(deep, x["deepest"])
    n = sum(ls.values())
    if n:
        print(f"| {d} | {n:,} | {ls[1] / n:.1%} | {ls[2] / n:.2%} | {ls[3] / n:.3%} | {ls[4] / n:.3%} | {deep} |")

# 4. each dive's contribution to the estimate (by_depth / dives of its stratum / pi), the tail
contrib = []    # (nodes, cpu, key, deepest, by_depth contribution)
for k, ds in dives.items():
    w = strata[k]["pi"]
    for x in ds:
        f = 1.0 / len(ds) / w
        contrib.append((x["nodes"] * f, x["cpu"] * f, k, x["deepest"], {d: v * f for d, v in x["by_depth"].items()},
                        {d: v * f for d, v in x["cpu_by_depth"].items()}))
exact_n = sum(s["exact_nodes"] / s["pi"] for s in strata.values())
exact_c = sum(s["exact_cpu"] / s["pi"] for s in strata.values())
dive_n = sum(c[0] for c in contrib)
dive_c = sum(c[1] for c in contrib)
print(f"\nrest = strata's own nodes {exact_n:,.0f} ({exact_c / H:,.1f} CPU h) + dives {dive_n:,.0f} ({dive_c / H:,.1f} CPU h)"
      f" = {exact_n + dive_n:,.0f} nodes, {(exact_c + dive_c) / H:,.1f} CPU h  (estimate.json rest: "
      f"{est['rest']['nodes']:,.0f}, {est['rest']['cpu'] / H:,.1f})")
for what, i, tot in (("nodes", 0, exact_n + dive_n), ("CPU", 1, exact_c + dive_c)):
    cs = sorted((c[i] for c in contrib), reverse=True)
    acc, out = 0.0, []
    for top in (1, 10, 100, 1000):
        out.append(f"top {top}: {sum(cs[:top]) / tot:.1%}")
    print(f"{what}: largest dive contributions as a share of the rest: " + ", ".join(out))

# the depth-6 mass: which dives carry it
print("\n| depth | rest nodes (dives) | dives contributing | largest one's share | top-10 share | rest CPU h (dives) |")
print("|---|---|---|---|---|---|")
alld = sorted({d for c in contrib for d in c[4]})
for d in alld:
    v = sorted((c[4].get(d, 0.0) for c in contrib), reverse=True)
    tot = sum(v)
    nz = sum(1 for x in v if x > 0)
    cpu = sum(c[5].get(d, 0.0) for c in contrib)
    print(f"| {d} | {tot:,.0f} | {nz:,} | {v[0] / tot:.1%} | {sum(v[:10]) / tot:.1%} | {cpu / H:,.1f} |")

# 5. the dives' per-step weights: 1/p of the chosen child against the number of children (what a
#    uniform choice would give): QPB weights far from uniform show as a spread of 1/p around n

# 6. per open node, by the stratum's depth: what an open node at that depth stands for (its own
#    bound plus its subtree), and so what the 28 never-leased units' 585 open nodes would add
print("\n| open node depth | open nodes (sum 1/pi) | subtree nodes per open node | CPU s per open node | branched (weighted) |")
print("|---|---|---|---|---|")
per = {}
for d in sorted(byd):
    ss = byd[d]
    npi = sum(1 / s["pi"] for s in ss)
    n = c = 0.0
    for s in ss:
        ds = dives.get(s["key"], [])
        n += (s["exact_nodes"] + (statistics.fmean(x["nodes"] for x in ds) if ds else 0.0)) / s["pi"]
        c += (s["exact_cpu"] + (statistics.fmean(x["cpu"] for x in ds) if ds else 0.0)) / s["pi"]
    br = sum(1 / s["pi"] for s in ss if s["survivors"] > 0) / npi
    per[d] = (n / npi, c / npi)
    print(f"| {d} | {npi:,.0f} | {n / npi:,.1f} | {c / npi:,.0f} | {br:.1%} |")
pend = [u for u in units if u["status"] != "merged"]
print("\nnever-leased units: " + ", ".join(f"{u['uid']} depth {u['depth_min']}-{u['depth_max']} {u['nodes']} nodes" for u in pend[:5]) + ", ...")
# a unit's open nodes split between its depths: depth_min where min == max, else half and half
add_n = add_c = 0.0
for u in pend:
    ds = [u["depth_min"]] if u["depth_min"] == u["depth_max"] else [u["depth_min"], u["depth_max"]]
    for d in ds:
        add_n += u["nodes"] / len(ds) * per[d][0]
        add_c += u["nodes"] / len(ds) * per[d][1]
print(f"their {sum(u['nodes'] for u in pend)} open nodes would add about {add_n:,.0f} nodes and {add_c / H:,.1f} CPU h "
      f"({add_n / est['total']['nodes']:.1%} and {add_c / est['total']['cpu']:.1%} of the estimate)")
