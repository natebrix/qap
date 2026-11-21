#!/usr/bin/env python3
"""The solved-sample estimate (coord/sample.py's design) simulated on sko42's known tree: the strata of
sko42-est-261121b's dives sized by their dive estimates, drawn with pi ~ size^0.5 (mix 0.2), each sampled
node priced at its true subtree CPU in sko42-261120b; the Horvitz-Thompson total against the true rest.

    python3 experiments/estimates_nov19/sample_sim.py
"""
import collections, gzip, json, os, random, statistics as st, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "experiments", "estimates", "sko42", "zone_nov19"))
from coord.sample import inclusion
from common import RUNS, load
T = load(os.path.join(RUNS, "sko42-261120b", "nodes.txt.gz")); cpu = T["cpu"]
strata, dives = {}, collections.defaultdict(list)
for l in gzip.open(os.path.join(RUNS, "sko42-est-261121b", "dives", "dives.jsonl.gz"), "rt"):
    x = json.loads(l)
    if x["type"] == "stratum": strata[x["key"]] = x
    elif x["type"] == "dive": dives[x["key"]].append(x)
rows = []
for k, s in strata.items():
    if k not in cpu: continue
    ds = dives.get(k, [])
    size = s["exact_cpu"] + (st.fmean(d["cpu"] for d in ds) if ds else 0.0)
    rows.append((s["depth"], k, max(size, 0.0), float(s.get("pi", 1)) or 1.0, cpu[k]))
true_rest = 539.3 - 46.7
by = collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
# the part no sample needs: strata at depths not sampled count at their dive estimate (as the estimator would), or exactly
def design(plan, draws=2000, seed=1):
    rng = random.Random(seed)
    pis = {d: inclusion([r[2] ** 0.5 for r in by[d]], min(n, len(by[d])), 0.2) for d, n in plan.items()}
    rest_other = sum(r[4] / r[3] for d in by if d not in plan for r in by[d])   # unsampled depths at their TRUE cost
    tots, costs = [], []
    for _ in range(draws):
        t, c = rest_other, 0.0
        for d, pl in pis.items():
            for r, p in zip(by[d], pl):
                if rng.random() < p:
                    t += r[4] / (r[3] * p); c += r[4]
        tots.append(t / 3600); costs.append(c / 3600)
    tots.sort()
    return st.fmean(tots), tots[int(.05 * draws)], tots[int(.95 * draws)], st.fmean(costs)
print(f"strata {len(rows):,}; the strata's true subtrees weighted: {sum(r[4]/r[3] for r in rows)/3600:,.0f} CPU h; the true rest {true_rest:,.0f}")
print("depths 5-6 at their true cost in every design (2 CPU h)" )
for plan in ({2: 20, 3: 50, 4: 50}, {2: 30, 3: 60, 4: 60}, {2: 40, 3: 80, 4: 80}, {2: 60, 3: 100, 4: 100}, {2: 169, 3: 100, 4: 100}):
    m, lo, hi, c = design(plan)
    print(f"  plan {plan}: total mean {m:,.0f} CPU h (5-95% {lo:,.0f}-{hi:,.0f}) against {true_rest:,.0f}; the sample's own CPU {c:,.0f} h")
