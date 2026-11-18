"""The heavy tail of tho40-est's dives: the largest contributions, the estimate without them, and a
delete-a-group jackknife over units (units are the sampling blocks: each unit's strata are kept or not
by its own draw), against combine()'s bootstrap.

    python3 -I tail.py
"""
import collections, math, os, random, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from load import load

units, strata, dives, samples, est = load()
H = 3600.0
done_n, done_c = est["done"]["nodes"], est["done"]["cpu"]


def part(k, ds=None):
    s = strata[k]
    ds = dives.get(k, []) if ds is None else ds
    n = s["exact_nodes"] + (statistics.fmean(x["nodes"] for x in ds) if ds else 0.0)
    c = s["exact_cpu"] + (statistics.fmean(x["cpu"] for x in ds) if ds else 0.0)
    return n / s["pi"], c / s["pi"]


def depth_parts(k):
    s = strata[k]
    ds = dives.get(k, [])
    n, c = collections.Counter(), collections.Counter()
    for d, v in s["exact_by_depth"].items():
        n[int(d)] += v / s["pi"]
    for d, v in s["exact_cpu_by_depth"].items():
        c[int(d)] += v / s["pi"]
    for x in ds:
        for d, v in x["by_depth"].items():
            n[d] += v / len(ds) / s["pi"]
        for d, v in x["cpu_by_depth"].items():
            c[d] += v / len(ds) / s["pi"]
    return n, c


rest_n = sum(part(k)[0] for k in strata)
rest_c = sum(part(k)[1] for k in strata)
print(f"rest {rest_n:,.0f} nodes, {rest_c / H:,.1f} CPU h; total {done_n + rest_n:,.0f} nodes, {(done_c + rest_c) / H:,.1f} CPU h")

# the largest dives
rows = []
for k, ds in dives.items():
    s = strata[k]
    for x in ds:
        f = 1.0 / len(ds) / s["pi"]
        rows.append((x["nodes"] * f, x["cpu"] * f, s["depth"], x["deepest"], len(ds), s["pi"], x["nodes"], k))
rows.sort(reverse=True)
print("\nthe 12 largest dives by nodes: | contribution, nodes | CPU h | stratum depth | deepest | dives in its stratum | pi | the dive's own estimate |")
for r in rows[:12]:
    print(f"| {r[0]:,.0f} | {r[1] / H:,.1f} | {r[2]} | {r[3]} | {r[4]} | {r[5]:.3f} | {r[6]:,.0f} |")

# without the largest dives: each dropped dive leaves its stratum's mean over the others
for drop in (1, 3, 10):
    gone = collections.defaultdict(set)
    for r in sorted(rows, key=lambda r: -r[0])[:drop]:
        ds = dives[r[7]]
        i = max(range(len(ds)), key=lambda i: (ds[i]["nodes"] if i not in gone[r[7]] else -1))
        gone[r[7]].add(i)
    n = c = 0.0
    for k in strata:
        ds = [x for i, x in enumerate(dives.get(k, [])) if i not in gone.get(k, ())]
        if dives.get(k) and not ds:
            ds = []
        pn, pc = part(k, ds)
        n += pn
        c += pc
    print(f"without the {drop} largest dive(s) by nodes: total {done_n + n:,.0f} nodes, {(done_c + c) / H:,.1f} CPU h")

# delete-a-group jackknife over units: G groups of units at random
G = 30
rng = random.Random(7)
merged = [u for u in units if u["status"] == "merged"]
rng.shuffle(merged)
grp = {u["uid"]: i % G for i, u in enumerate(merged)}
gn, gc = [0.0] * G, [0.0] * G
gdn = [collections.Counter() for _ in range(G)]
gdc = [collections.Counter() for _ in range(G)]
for k, s in strata.items():
    g = grp[s["unit"]]
    pn, pc = part(k)
    gn[g] += pn
    gc[g] += pc
    dn, dc = depth_parts(k)
    gdn[g].update(dn)
    gdc[g].update(dc)
tn, tc = sum(gn), sum(gc)
jn = [(tn - gn[g]) * G / (G - 1) for g in range(G)]
jc = [(tc - gc[g]) * G / (G - 1) for g in range(G)]
se = lambda xs, t: math.sqrt((G - 1) / G * sum((x - statistics.fmean(xs)) ** 2 for x in xs))
sn, sc = se(jn, tn), se(jc, tc)
print(f"\njackknife over units ({G} groups): rest nodes se {sn:,.0f} ({sn / tn:.0%}), CPU se {sc / H:,.0f} h ({sc / tc:.0%})")
print(f"  total +-2 se: nodes {done_n + tn - 2 * sn:,.0f} to {done_n + tn + 2 * sn:,.0f}; "
      f"CPU {(done_c + tc - 2 * sc) / H:,.0f} to {(done_c + tc + 2 * sc) / H:,.0f} h")
print(f"  combine()'s bootstrap 5-95%: nodes {est['total']['nodes_5_95'][0]:,.0f} to {est['total']['nodes_5_95'][1]:,.0f}; "
      f"CPU {est['total']['cpu_5_95'][0] / H:,.0f} to {est['total']['cpu_5_95'][1] / H:,.0f} h")
print("\n| depth | rest nodes | jackknife se | rest CPU h | se |")
print("|---|---|---|---|---|")
for d in sorted(set().union(*[set(x) for x in gdn])):
    vn = [gdn[g][d] for g in range(G)]
    vc = [gdc[g][d] for g in range(G)]
    t1, t2 = sum(vn), sum(vc)
    j1 = [(t1 - v) * G / (G - 1) for v in vn]
    j2 = [(t2 - v) * G / (G - 1) for v in vc]
    s1, s2 = se(j1, t1), se(j2, t2)
    print(f"| {d} | {t1:,.0f} | {s1:,.0f} ({s1 / t1 if t1 else 0:.0%}) | {t2 / H:,.1f} | {s2 / H:,.1f} ({s2 / t2 if t2 else 0:.0%}) |")
