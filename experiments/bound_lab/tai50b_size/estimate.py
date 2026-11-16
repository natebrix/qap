"""The estimate of tai50b's tree from the stage runs:  python3 estimate.py   (mkreport.sh runs it -> logs/estimate.txt)

Measured: the root (runs/root_long), the depth-1 child on the path of the best known solution, 31:7 (stratum A),
and a random 17 of the other 49 depth-1 children (stratum B; runs/L1_*), and the path of the best known solution
(the "spine", runs/spine_D) down to depth DS = 8.
Estimated by Knuth dives (bnbdrv -e -r STATE, one finished dive per line of dives.txt = an unbiased estimate of
the nodes / CPU-s of all subtrees rooted at STATE's nodes):
  LA_offD   the off-path children at depth D of the spine node at depth D-1 (D = 2..8)
  LA_tail8  all children (depth 9) of the spine node at depth 8
  LB, LB2   the depth-2 children of the sampled B nodes that branched (two groups of parents)
  T_A = sum_D mean(LA_offD) + 7 spine nodes + mean(LA_tail8)
  T_i (a B parent) = mean over its set's dives of (dive estimate x [the dive started under i])

E1 (the sampling design alone):   51 + T_A + (49/17) sum_{i in B sample} T_i
E2 (post-stratified by landscape.py: the depth-1 children holding a solution within 1M of the best known value,
    31:0 and 31:32, are "twins" of 31:7 and are given T_A each [e]; 31:43 and 31:42 (within 2.6M) are given
    31:43's T_i; the other 45, of which 15 are in the sample, are expanded 45/15):
                                  51 + 3 T_A + 2 T_31:43 + (45/15) sum_{i in B sample, not 31:32, 31:43} T_i
Intervals: bootstrap, 5th-95th percentile: dives resampled within every set, B parents resampled with
replacement (no finite-population credit)."""
import sys, glob, os, random, csv, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, dives as DV
H = os.path.dirname(os.path.abspath(__file__))
cpu = lambda d: d["fcpu"] + d["ccpu"]
TWIN, NEAR = ["31:32"], ["31:43"]            # in the sample; landscape: twins 31:0, 31:32; near 31:43, 31:42

def dive_rows(setname):
    """(depth-1 fix, node estimate, CPU estimate, depth reached) per finished dive of a set"""
    rows = []
    for d in sorted(glob.glob(f"{H}/runs/{setname}/s*")):
        if not os.path.exists(d + "/dives.txt"): continue
        est = [ln.split() for ln in open(d + "/dives.txt")]
        ch = [c for f, c in DV.dives(f"{H}/runs/{setname}") if os.path.dirname(f) == d][:len(est)]
        assert len(ch) == len(est), (d, len(ch), len(est))
        for e, c in zip(est, ch):
            rows.append((next(x for x in c[0]["fixes"] if x.startswith("31:")), float(e[1]), float(e[2]), int(e[0])))
    return rows

def main():
    root = N.load([f"{H}/runs/root_long/dnn.txt"])[0]
    d1 = [d for d in N.load(sorted(glob.glob(f"{H}/runs/L1_*/dnn.txt"))) if d["depth"] == 1]
    A = [d for d in d1 if d["onpath"]][0]; B = [d for d in d1 if not d["onpath"]]; nB = len(B)
    Bfix = [d["fixes"][0] for d in B]; cpuB = {d["fixes"][0]: cpu(d) for d in B}
    B2 = [b for b in Bfix if b not in TWIN + NEAR]
    sp = {d["depth"]: d for d in N.load(sorted(glob.glob(f"{H}/runs/spine_*/dnn.txt"))) if d["onpath"]}
    la_sets = [f"LA_off{d}" for d in range(2, 9)] + ["LA_tail8"]
    la = {x: dive_rows(x) for x in la_sets}
    lb = {os.path.basename(p): dive_rows(os.path.basename(p)) for p in sorted(glob.glob(f"{H}/runs/LB*"))}
    missing = [x for x in la_sets if not la[x]]
    spine = [sp[d] for d in range(2, 9) if d in sp]
    print(f"root: {root['fit'] + root['cit']} iterations, {cpu(root):,.0f} CPU-s, certified {root['final']:,.0f}")
    print(f"depth 1 [m]: 31:7 branched ({cpu(A):.0f} CPU-s); B sample {nB} of 49: {sum(d['fath'] for d in B)} fathomed, "
          f"{sum(1 - d['fath'] for d in B)} branched; mean CPU-s {st.mean(cpuB.values()):.0f}")
    print("spine [m]: depth: shortfall", {d: round(x["cutoff"] - x["final"]) for d, x in sorted(sp.items())})
    print("dives per set:", {x: len(r) for x, r in {**la, **lb}.items()}, "| sets without dives:", missing)
    print("deepest dives:", {x: max((r[3] for r in rows), default=None) for x, rows in {**la, **lb}.items()})

    def parts(la, lb, k):            # k = 1 nodes, 2 CPU-s
        TA = sum(st.mean(r[k] for r in x) for x in la.values() if x) + (len(spine) if k == 1 else sum(cpu(d) for d in spine))
        share = {}
        for x in lb.values():
            for r in x: share[r[0]] = share.get(r[0], 0.0) + r[k] / len(x)
        return TA, share
    def est(la, lb, Bs, B2s):
        out = []
        for k in (1, 2):
            TA, sh = parts(la, lb, k)
            top1 = 51.0 if k == 1 else cpu(root) + cpu(A) + 49.0 / nB * sum(cpuB[b] for b in Bs)
            e1 = top1 + TA + 49.0 / nB * sum(sh.get(b, 0.0) for b in Bs)
            top2 = 51.0 if k == 1 else cpu(root) + cpu(A) + sum(cpuB[b] for b in TWIN + NEAR) * 2 + 45.0 / len(B2) * sum(cpuB[b] for b in B2s)
            e2 = top2 + 3 * TA + 2 * sum(sh.get(b, 0.0) for b in NEAR) + 45.0 / len(B2) * sum(sh.get(b, 0.0) for b in B2s)
            out += [e1, e2, TA]
        return out
    e = est(la, lb, Bfix, B2)
    TAn, shn = parts(la, lb, 1); TAc, shc = parts(la, lb, 2)
    print(f"\nT_A (below 31:7) [e]: {TAn:,.0f} nodes, {TAc:,.0f} CPU-s; by set (nodes, CPU-s):",
          {x: (round(st.mean(r[1] for r in rows)), round(st.mean(r[2] for r in rows))) for x, rows in la.items() if rows})
    print("T_i of the branched B parents [e] (nodes, CPU-s):", {b: (round(shn.get(b, 0)), round(shc.get(b, 0))) for b in Bfix if not next(d for d in B if d["fixes"][0] == b)["fath"]})
    rng = random.Random(1); reps = []
    for _ in range(4000):
        la_ = {x: [rng.choice(r) for _ in r] for x, r in la.items()}
        lb_ = {x: [rng.choice(r) for _ in r] for x, r in lb.items()}
        reps.append(est(la_, lb_, [rng.choice(Bfix) for _ in Bfix], [rng.choice(B2) for _ in B2]))
    def line(name, i, unit):
        v = sorted(r[i] for r in reps); q = lambda p: v[int(p * (len(v) - 1))]
        s = f"{name}: {e[i]:,.0f} {unit}  (bootstrap 5%-95%: {q(.05):,.0f} to {q(.95):,.0f})"
        if unit == "CPU-s": s += f"  = {e[i] / 3600:,.0f} core-hours ({q(.05) / 3600:,.0f} to {q(.95) / 3600:,.0f})"
        print(s)
    print()
    line("E1 nodes", 0, "nodes"); line("E1 CPU", 3, "CPU-s"); line("E2 nodes", 1, "nodes"); line("E2 CPU", 4, "CPU-s")
    line("T_A nodes", 2, "nodes"); line("T_A CPU", 5, "CPU-s")

    # per depth, for E2: the estimator's per-level tables as pool.py splits them by dive (runs/SET/levels.csv)
    lev = {0: [1.0, cpu(root)], 1: [50.0, cpu(A) + sum(cpuB[b] for b in TWIN + NEAR) * 2 + 45.0 / len(B2) * sum(cpuB[b] for b in B2)]}
    for d in spine:
        L = lev.setdefault(d["depth"], [0.0, 0.0]); L[0] += 3; L[1] += 3 * cpu(d)
    def add(name, mult, parent=None):
        p = f"{H}/runs/{name}/levels.csv"
        if not os.path.exists(p): return
        for r in csv.DictReader(open(p)):
            L = lev.setdefault(int(r["level"]), [0.0, 0.0]); L[0] += mult * float(r["nodes_mean"]); L[1] += mult * float(r["cpu_mean"])
    for x in la_sets: add(x, 3.0)
    # the B sets by level are not split by parent: LB (31:4, 31:15, 31:32, 31:43) enters with the weight that
    # reproduces E2's total for its parents, LB2 (31:5, 31:11) with 45/15
    tb = sum(shc.get(b, 0.0) for b in Bfix if b in ("31:4", "31:15", "31:32", "31:43"))
    wb = (2 * sum(shc.get(b, 0.0) for b in NEAR) + 45.0 / len(B2) * sum(shc.get(b, 0.0) for b in ("31:4", "31:15"))) / tb if tb else 0.0
    add("LB", wb); add("LB2", 45.0 / len(B2))
    tn = sum(v[0] for v in lev.values()); tc = sum(v[1] for v in lev.values())
    print("\nE2 by depth [e] (depths 0, 1 and the spine measured):\n\n| depth | nodes | CPU-s | share of CPU |\n|---|---|---|---|")
    for dpt in sorted(lev): print(f"| {dpt} | {lev[dpt][0]:,.0f} | {lev[dpt][1]:,.0f} | {100 * lev[dpt][1] / tc:.1f}% |")
    print(f"| all | {tn:,.0f} | {tc:,.0f} | |")

if __name__ == "__main__": main()
