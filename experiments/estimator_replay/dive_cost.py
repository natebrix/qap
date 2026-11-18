"""What a dives run costs after a pause, on a finished tree (docs/estimates_plan.md, step 3 sizing).

The solve is approximated as breadth-first by depth (shallowest first, then by shortfall) until it
has spent FRACTION of the tree's CPU; its open nodes are the strata. Per stratum the estimator
always bounds the node and, if it branches, all its children (setup), then at least PILOT dives
before its budget is looked at. A dive re-processes its survivor and bounds all of its children,
level after level. Costs are the recorded CPU of those nodes; dives are drawn as the estimator
draws them (exp(b x) over the survivors, 0.1 uniform mixed in), REPS times.
    python3 -I dive_cost.py TRACE --fraction 0.16 --pilot 2 --b 0.4 [--slots 192]"""
import argparse, math, os, random, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay

ap = argparse.ArgumentParser()
ap.add_argument("trace")
ap.add_argument("--fraction", type=float, default=0.16)
ap.add_argument("--pilot", type=int, default=2)
ap.add_argument("--b", type=float, default=0.4)
ap.add_argument("--slots", type=int, default=192)
ap.add_argument("--reps", type=int, default=20)
a = ap.parse_args()
T = replay.Tree(replay.load_shards(a.trace) if os.path.isdir(a.trace) else replay.load_trace(a.trace))
C = T.C[T.root]
order = sorted(range(len(T.depth)), key=lambda i: (T.depth[i], -(T.cutoff[i] - T.best[i])))
# the solve: nodes processed in that order, a node only after its parent
done, spent = set(), 0.0
frontier = [T.root]
import heapq
h = [(T.depth[T.root], 0.0, T.root)]
while h and spent < a.fraction * C:
    _, _, i = heapq.heappop(h)
    done.add(i); spent += T.cpu[i]
    for j in T.kids[i]:
        heapq.heappush(h, (T.depth[j], -(T.cutoff[j] - T.best[j]), j))
open_nodes = [i for _, _, i in h]
setup = sum(T.cpu[s] + (sum(T.cpu[c] for c in T.kids[s]) if T.branched[s] else 0.0) for s in open_nodes)
strata_surv = [s for s in open_nodes if T.branched[s] and any(T.branched[c] for c in T.kids[s])]
q = ("x", a.b)
rng = random.Random(1)
def dive_cost(s):
    surv, p, cp, _, _ = T.surv_probs(s, q)
    k = replay.draw(rng, cp); j = surv[k]; cost = 0.0
    while True:
        cost += T.cpu[j] + sum(T.cpu[c] for c in T.kids[j])   # j again, and its children
        surv, p, cp, _, _ = T.surv_probs(j, q)
        if not surv: return cost
        j = surv[replay.draw(rng, cp)]
pilots = []
for _ in range(a.reps):
    pilots.append(sum(dive_cost(s) for s in strata_surv for _ in range(a.pilot)))
pm = statistics.fmean(pilots)
hh = lambda s: s / 3600
print(f"{os.path.basename(a.trace)}: {C / 3600:,.0f} CPU h. Solve to {a.fraction:.0%}: {hh(spent):,.0f} CPU h, "
      f"{len(open_nodes):,} open nodes, {len(strata_surv):,} with children that branch")
print(f"setup (each open node and, where it branches, its children): {hh(setup):,.0f} CPU h ({setup / C:.1%} of the tree)")
print(f"{a.pilot} pilot dives a stratum with survivors: {hh(pm):,.0f} CPU h on average ({pm / C:.1%}); "
      f"10-90%: {hh(sorted(pilots)[int(.1 * len(pilots))]):,.0f}-{hh(sorted(pilots)[int(.9 * len(pilots))]):,.0f}")
tot = spent + setup + pm
print(f"solve + setup + pilot: {hh(tot):,.0f} CPU h ({tot / C:.0%} of the tree), {hh(tot) / a.slots:.1f} h on {a.slots} slots")
