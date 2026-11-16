"""The stratified estimate:  python3 estimate.py  (reads runs/root_long, runs/L1_*, runs/LA, runs/LB)
Depths 0 and 1 are measured: the root, the depth-1 child on the best known solution's path (stratum A, certain)
and a random NB of the other 49 (stratum B). Below depth 1 Knuth dives (bnbdrv -e -r, level 2) from the open
depth-2 nodes: set LA under the on-path child, set LB under the sampled B nodes that branched.
  nodes = 1 + 50 + N_A + (49/NB) * S_B      CPU likewise, with the measured CPU of the root and the depth-1 nodes
Interval: bootstrap (5th-95th percentile): B parents resampled with replacement (no finite-population credit),
dives resampled within each set; a parent's share of S_B is the mean over LB's dives of (dive estimate x [dive under it])."""
import sys, glob, os, random, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, dives as DV
H = os.path.dirname(os.path.abspath(__file__))
def dive_rows(setname):
    """(parent fix, node estimate, cpu estimate, depth reached, chain) per finished dive"""
    rows = []
    for d in sorted(glob.glob(f"{H}/runs/{setname}/s*")):
        if not os.path.exists(d + "/dives.txt"): continue
        est = [ln.split() for ln in open(d + "/dives.txt")]
        ch = [c for f, c in DV.dives(f"{H}/runs/{setname}") if os.path.dirname(f) == d]
        ch = [c for c in ch if c[-1]["fath"]][:len(est)] if len([c for c in ch if c[-1]["fath"]]) >= len(est) else ch[:len(est)]
        assert len(ch) == len(est), (d, len(ch), len(est))
        for e, c in zip(est, ch): rows.append((next(x for x in c[0]["fixes"] if x.startswith("31:")), float(e[1]), float(e[2]), int(e[0]), c))
    return rows
def main():
    root = N.load([f"{H}/runs/root_long/dnn.txt"])[0]
    d1 = [d for d in N.load(sorted(glob.glob(f"{H}/runs/L1_*/dnn.txt"))) if d["depth"] == 1]
    A = [d for d in d1 if d["onpath"]]; B = [d for d in d1 if not d["onpath"]]; nB = len(B)
    cpu = lambda d: d["fcpu"] + d["ccpu"]
    print(f"root: {root['fit'] + root['cit']} iterations, {cpu(root):.0f} CPU-s, bound {root['final']:,.0f}")
    print(f"depth 1: on-path child {'branched' if not A[0]['fath'] else 'FATHOMED'} ({cpu(A[0]):.0f} CPU-s); B sample {nB} of 49: "
          f"{sum(d['fath'] for d in B)} fathomed, {sum(1 - d['fath'] for d in B)} branched; CPU-s mean {st.mean(map(cpu, B)):.0f}")
    # the B stratum's dives may be in several sets (LB, LB2, ...: disjoint groups of parents, one state file each)
    sp = [d for d in N.load(sorted(glob.glob(f"{H}/runs/spine_*/dnn.txt"))) if d["onpath"]]
    la_sets = sorted(os.path.basename(p) for p in glob.glob(f"{H}/runs/LA_*"))
    la = [dive_rows(x) for x in la_sets]
    used = [d for d in sp if f"LA_off{d['depth']}" in la_sets or (f"LA_off{d['depth']+1}" in la_sets) or "LA_tail%d" % d["depth"] in la_sets]
    print("spine nodes measured:", [(d["depth"], round(d["cutoff"] - d["final"])) for d in sp], "; counted:", len(used))
    print("LA sets:", {x: len(r) for x, r in zip(la_sets, la)})
    spine_n = len(used); spine_cpu = sum(cpu(d) for d in used)
    lb = [dive_rows(os.path.basename(p)) for p in sorted(glob.glob(f"{H}/runs/LB*"))]
    print(f"dives: LB sets {[len(x) for x in lb]}; deepest (LA) {sorted(r[3] for x in la for r in x)[-5:]}, (LB) {sorted(r[3] for x in lb for r in x)[-5:]}")
    top_nodes = 51.0
    def total(la, lb, Bs, k):   # k = 1 nodes, 2 CPU: (N_A, S_B over the parents Bs)
        NA = sum(st.mean(r[k] for r in x) for x in la if x) + (spine_n if k == 1 else spine_cpu)
        share = {}
        for x in lb:
            for r in x: share[r[0]] = share.get(r[0], 0.0) + r[k] / len(x)
        return NA, sum(share.get(b, 0.0) for b in Bs)
    Bfix = [d["fixes"][0] for d in B]; cpuB = {d["fixes"][0]: cpu(d) for d in B}
    def est(la, lb, Bs):
        out = []
        for k in (1, 2):
            NA, SB = total(la, lb, Bs, k)
            top = top_nodes if k == 1 else cpu(root) + cpu(A[0]) + 49.0 / nB * sum(cpuB[b] for b in Bs)
            out.append((top, NA, 49.0 / nB * SB))
        return out
    e = est(la, lb, Bfix)
    for name, (top, NA, SBx) in zip(("nodes", "CPU-s"), e):
        print(f"{name}: depths 0-1 {top:,.0f} + under the on-path child {NA:,.0f} + under the other 49 {SBx:,.0f} = {top + NA + SBx:,.0f}")
    rng = random.Random(1); reps = []
    for _ in range(4000):
        la_ = [[rng.choice(x) for _ in x] for x in la]
        lb_ = [[rng.choice(x) for _ in x] for x in lb]
        Bs = [rng.choice(Bfix) for _ in Bfix]
        reps.append([sum(x) for x in est(la_, lb_, Bs)])
    for k, name in ((0, "nodes"), (1, "CPU-s")):
        v = sorted(r[k] for r in reps); q = lambda p: v[int(p * (len(v) - 1))]
        print(f"{name}: bootstrap 5% {q(.05):,.0f}  median {q(.5):,.0f}  95% {q(.95):,.0f}" + (f"   = {q(.05)/3600:,.0f} to {q(.95)/3600:,.0f} core-hours (point {sum(e[1])/3600:,.0f})" if k else ""))
    # per level, from the estimator's own tables (pool.py), is printed by pool.sh
if __name__ == "__main__": main()
