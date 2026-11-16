#!/usr/bin/env python3
"""The headline table: the cut arm's tree against the base, in c7a CPU (python3 scripts/final.py [--all]).

Base reference: the laptop study ../../tho40b_analysis/knuth_tho40 (3,322 dives, incumbent 240,518, dives from depth 3),
its laptop seconds x 2.23 = c7a seconds (its README). This box under the run's load: the reproduction run's top (40 nodes,
the same iterations) took 1.73 x the laptop's CPU, so c7a seconds = 2.23 / 1.73 = 1.29 x this box's."""
import math, os, random, re, statistics as st, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import load, estimate, STRATA, HERE, fmt
LAP = HERE + "/../../tho40b_analysis/knuth_tho40"
C7A_PER_LAPTOP = 2.23
def box_factor():
    g = lambda p: float(re.search(r"Total Time: (\S+) seconds", open(p, errors="replace").read())[1])
    return g(HERE + "/runs/repro_s0/est.out") / g(LAP + "/s0/est.out")
F = C7A_PER_LAPTOP / box_factor()      # c7a s per s here
D = 86400.0
use_all = "--all" in sys.argv
lap = {}
for c in STRATA:
    d = [tuple(map(float, l.split())) for l in open(f"{LAP}/s{c}/dives.txt") if l.strip()]
    N, T = [x[1] for x in d], [x[2] * C7A_PER_LAPTOP for x in d]
    lap[c] = dict(n=len(d), N=st.fmean(N), T=st.fmean(T), seN=st.stdev(N) / math.sqrt(len(N)), seT=st.stdev(T) / math.sqrt(len(T)))
cut, base = load("cut06", use_all), load("base", use_all)
res = {}
for name, S in (("cut", cut), ("base", base)):
    res[name] = {k: estimate(S, k) for k in ("nodes", "cpu", "it1", "itc")}
print(f"c7a CPU s per CPU s here (under load): {F:.2f}\n")
print("| | dives | nodes | standard error | 80% band | CPU s here | standard error | 80% band | c7a CPU-days | standard error | 80% band |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for name, S in (("cut arm (zone 0.06)", "cut"), ("base, this box", "base")):
    n, c = res[S]["nodes"], res[S]["cpu"]
    nd = sum(p[2] for p in n[3].values())
    print(f"| {name} | {nd} | {fmt(n[0] + 1)} | {fmt(n[1])} | {fmt(n[2][0])} to {fmt(n[2][1])} | {fmt(c[0])} | {fmt(c[1])} | {fmt(c[2][0])} to {fmt(c[2][1])} | "
          f"{c[0] * F / D:,.0f} | {c[1] * F / D:,.0f} | {c[2][0] * F / D:,.0f} to {c[2][1] * F / D:,.0f} |")
LN, LT = sum(l["N"] for l in lap.values()), sum(l["T"] for l in lap.values())
LseN, LseT = math.sqrt(sum(l["seN"] ** 2 for l in lap.values())), math.sqrt(sum(l["seT"] ** 2 for l in lap.values()))
print(f"| base, the laptop study | {sum(l['n'] for l in lap.values())} | {fmt(LN)} | {fmt(LseN)} | | | | | {LT / D:,.0f} | {LseT / D:,.0f} | |")
def ratio(a, sa, b, sb):
    r = a / b
    return r, r * math.sqrt((sa / a) ** 2 + (sb / b) ** 2)
c, n = res["cut"]["cpu"], res["cut"]["nodes"]
print("\n| ratio cut / base | value | standard error (delta method) | base / cut |\n|---|---|---|---|")
r, s = ratio(c[0] * F, c[1] * F, LT, LseT); print(f"| CPU, against the laptop study | {r:.3f} | {s:.3f} | {1 / r:.1f}x |")
r, s = ratio(n[0], n[1], LN, LseN); print(f"| nodes, against the laptop study | {r:.4f} | {s:.4f} | {1 / r:.0f}x |")
if all(k in base and base[k]["dives"] for k in STRATA):
    b, bn = res["base"]["cpu"], res["base"]["nodes"]
    r, s = ratio(c[0], c[1], b[0], b[1]); print(f"| CPU, against the base dives here | {r:.3f} | {s:.3f} | {1 / r:.1f}x |")
    r, s = ratio(n[0], n[1], bn[0], bn[1]); print(f"| nodes, against the base dives here | {r:.4f} | {s:.4f} | {1 / r:.0f}x |")
print("\n| stratum | cut: dives | nodes (se) | c7a CPU-days (se) | base here: dives | nodes (se) | c7a CPU-days (se) | base, laptop study: dives | nodes (se) | c7a CPU-days (se) | CPU cut / base here | CPU cut / laptop study | nodes cut / laptop study |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for k in STRATA:
    pc, pn = res["cut"]["cpu"][3].get(k), res["cut"]["nodes"][3].get(k)
    bc, bn = res["base"]["cpu"][3].get(k), res["base"]["nodes"][3].get(k)
    l = lap[k]
    f = lambda p, s=1.0: f"{fmt(p[0] * s)} ({fmt(p[1] * s)})" if p and p[2] else "-"
    f2 = lambda p, s=1.0: f"{p[0] * s:.1f} ({p[1] * s:.1f})" if p and p[2] else "-"
    print(f"| 26:{k} | {pn[2] if pn else 0} | {f(pn)} | {f2(pc, F / D)} | {bn[2] if bn else 0} | {f(bn)} | {f2(bc, F / D)} | {l['n']} | {fmt(l['N'])} ({fmt(l['seN'])}) | {l['T'] / D:.1f} ({l['seT'] / D:.1f}) | "
          + (f"{pc[0] / bc[0]:.2f}" if pc and bc and pc[2] and bc[2] else "-") + " | " + (f"{pc[0] * F / l['T']:.2f}" if pc and pc[2] else "-") + " | " + (f"{pn[0] / l['N']:.4f}" if pn and pn[2] else "-") + " |")
c1, c2 = res["cut"]["it1"][0], res["cut"]["itc"][0]
print(f"\ncut arm iterations: first solves {fmt(c1)}, phases {fmt(c2)}: {c2 / (c1 + c2):.0%} are phase iterations")
