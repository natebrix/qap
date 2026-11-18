"""The estimate run's flow replayed on a finished tree, with the three fixes of README round 3.

1. the solve, to SOLVE of the tree's CPU: breadth-first (bfs: the shallowest open node first, as the
   coordinator's queue roughly goes) or best-first (best: the open node whose parent has the largest
   shortfall first, the heavy part of the tree solved exactly);
2. its open nodes SAMPLED (Poisson sampling, node i kept with probability pi_i = min(1, m w_i / sum w),
   m = SAMPLE x the open nodes): w uniform (unif), or g^2 of the node's own QPB bound (qpb: QPB is
   cheap next to the DNN bound, so the estimator can compute it for every open node), or exp(b x) of
   its parent's shortfall over the depth's median (parent); each sampled node weighted 1/pi_i
   (Horvitz-Thompson: unbiased);
3. each sampled node bounded, and dives below it (knuth: one child a level, g^1.5 of the child's own
   QPB bound; surv: all children a level, exp(0.4 x) over the survivors), a pilot of 1 dive per
   sampled node with children, then Neyman, until the TOTAL budget (solve + everything after it)
   is spent; the dives a node's last one is cut at are kept (the replay cannot cut a dive).

    python3 -I flow.py TRACE [--solve 0.1] [--sample 0.05] [--total 0.2] [--reps 20] [--grid default]

A node bounded more than once (by several dives) is paid once, as with the estimator's memo of
each stratum's children; the dives' own nodes are all tree nodes, so their cost is tree work.
"""
import argparse, heapq, math, os, random, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay


def solve_phase(T, frac, order):
    C = T.C[T.root]
    key = (lambda i, par: (T.depth[i], -(T.cutoff[par] - T.best[par]) if par is not None else 0.0)) if order == "bfs" \
        else (lambda i, par: (-(T.cutoff[par] - T.best[par]) / T.med_short(T.depth[par]) if par is not None else 0.0,
                              T.depth[i]))
    h = [(key(T.root, None), T.root, None)]
    spent = 0.0
    done_n = done_c = 0.0
    while h and spent < frac * C:
        _, i, _ = heapq.heappop(h)
        spent += T.cpu[i]
        done_n += T.cnt[i]
        done_c += T.cpu[i]
        for j in T.kids[i]:
            heapq.heappush(h, (key(j, i), j, i))
    return [(i, par) for _, i, par in h], done_n, done_c, spent


def run_flow(T, rng, solve, sample, total, order, pi_mode, dive_mode):
    C = T.C[T.root]
    openn, done_n, done_c, spent = solve_phase(T, solve, order)
    if not openn:
        return done_n, done_c, spent, 0
    # 2. sample the open nodes
    if pi_mode == "unif":
        w = [1.0] * len(openn)
    elif pi_mode == "qpb":
        w = [max(T.g(T.qpb[i]), 0.0) ** 2 + 1e-12 for i, _ in openn]
    else:   # parent
        w = [math.exp(min(0.4 * max(T.cutoff[p] - T.best[p], 0.0) / T.med_short(T.depth[p]), 50.0)) for _, p in openn]
    m = max(1.0, sample * len(openn))
    # inclusion probabilities min(1, m w / sum w), the cap's excess spread again (a few passes)
    pi = [0.0] * len(w)
    left, cap = list(range(len(w))), set()
    for _ in range(10):
        sw = sum(w[k] for k in left)
        mm = m - len(cap)
        if sw <= 0 or mm <= 0:
            break
        new_cap = [k for k in left if mm * w[k] / sw >= 1.0]
        if not new_cap:
            for k in left:
                pi[k] = mm * w[k] / sw
            break
        for k in new_cap:
            pi[k] = 1.0
            cap.add(k)
        left = [k for k in left if k not in cap]
    for k in cap:
        pi[k] = 1.0
    chosen = [(openn[k][0], pi[k]) for k in range(len(openn)) if pi[k] > 0 and rng.random() < pi[k]]
    paid = set()
    est_n, est_c = done_n, done_c
    strata = []
    q = ("x", 0.4) if dive_mode == "surv" else 1.5
    mode = "surv" if dive_mode == "surv" else "knuth"
    for i, p_i in chosen:
        paid.add(i)
        spent += T.cpu[i]
        est_n += T.cnt[i] / p_i
        est_c += T.cpu[i] / p_i
        if not T.branched[i]:
            continue
        if mode == "surv":
            for c in T.kids[i]:
                if c not in paid:
                    paid.add(c)
                    spent += T.cpu[c]
            surv, pr, cp, fn, fc = T.surv_probs(i, q)
            est_n += fn / p_i
            est_c += fc / p_i
            if not surv:
                continue
        strata.append([i, p_i, 0, 0.0, 0.0, 0.0, 0.0])   # node, pi, dives, sum n, sum c, sum c^2, cost EMA
    def one(s):
        nonlocal spent
        i = s[0]
        if mode == "surv":
            surv, pr, cp, _, _ = T.surv_probs(i, q)
            k = replay.draw(rng, cp)
            en, ec, cost, _ = replay.dive(T, surv[k], 1.0 / pr[k], "surv", q, rng, paid)
        else:
            ks, pr, cp = T.knuth_probs(i, q)
            k = replay.draw(rng, cp)
            en, ec, cost, _ = replay.dive(T, ks[k], 1.0 / pr[k], "knuth", q, rng, paid)
        spent += cost
        s[2] += 1; s[3] += en; s[4] += ec; s[5] += ec * ec
        s[6] = cost if s[2] == 1 else 0.7 * s[6] + 0.3 * cost
    for s in strata:                     # the pilot: one dive each, whatever the budget
        one(s)
    n = len(strata)
    # Knuth's dives go round the strata in turn: Neyman's allocation reads low with their heavy tails
    # (the estimator's QAP_EST_ALLOC=auto; nug15: -6% with Neyman, -0.7% equal)
    while strata and mode == "knuth" and spent < total * C and n < 200000:
        for s in strata:
            one(s); n += 1
            if spent >= total * C:
                break
    while strata and spent < total * C and n < 200000:
        wts = []
        for s in strata:
            mu = s[4] / s[2]; sd = math.sqrt(max(s[5] / s[2] - mu * mu, 0.0)); sd = max(sd, mu)
            wts.append(sd / s[1] / math.sqrt(max(s[6], 1.0)) / s[2])
        cpw = replay.cum(wts)
        for _ in range(max(1, len(strata) // 4)):
            one(strata[replay.draw(rng, cpw)]); n += 1
            if spent >= total * C:
                break
    for s in strata:
        est_n += s[3] / s[2] / s[1]
        est_c += s[4] / s[2] / s[1]
    return est_n, est_c, spent, len(chosen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trace")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--grid", default="default")
    ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args()
    T = replay.Tree(replay.load_shards(a.trace) if os.path.isdir(a.trace) else replay.load_trace(a.trace))
    N, C = T.N[T.root], T.C[T.root]
    print(f"# {a.trace}: {N:,} nodes, {C / 3600:,.1f} CPU h")
    print("| solve | order | sample | pi | dives | total budget | open nodes sampled | median est/true CPU | 10-90% | "
          "within 20% | within 2x | CPU spent / true |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    rng = random.Random(a.seed)
    grids = {
        "default": [(s, o, m, p, d, t) for s in (0.05, 0.1) for o in ("bfs", "best") for m in (0.02, 0.1)
                    for p in ("unif", "qpb", "parent") for d in ("knuth", "surv") for t in (0.2,)],
        "bushy": [(s, o, m, p, "knuth", t) for s in (0.02, 0.05, 0.1) for o in ("bfs",) for m in (0.05, 0.2, 1.0)
                  for p in ("unif", "qpb") for t in (0.1, 0.2, 0.3)],
        "budget": [(s, o, m, p, d, t) for s in (0.1,) for o in ("best",) for m in (0.02, 0.05, 0.1) for p in ("qpb", "parent")
                   for d in ("knuth", "surv") for t in (0.1, 0.2, 0.3)],
    }
    for solve, order, sample, pim, dm, total in grids[a.grid]:
        rc, sp, ch = [], [], []
        for _ in range(a.reps):
            en, ec, spent, nch = run_flow(T, rng, solve, sample, total, order, pim, dm)
            rc.append(ec / C); sp.append(spent / C); ch.append(nch)
        rc.sort()
        pct = lambda xs, f: xs[min(int(f * len(xs)), len(xs) - 1)]
        w20 = sum(0.8 <= x <= 1.2 for x in rc) / len(rc); w2 = sum(0.5 <= x <= 2 for x in rc) / len(rc)
        print(f"| {solve:g} | {order} | {sample:g} | {pim} | {dm} | {total:g} | {statistics.median(ch):,.0f} | "
              f"{pct(rc, .5):.2f} | {pct(rc, .1):.2f}-{pct(rc, .9):.2f} | {w20:.0%} | {w2:.0%} | {statistics.fmean(sp):.3f} |",
              flush=True)


if __name__ == "__main__":
    main()
