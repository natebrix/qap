#!/usr/bin/env python3
"""Whole-tree estimates from the dives of runs/ARM_s<c>_<tag>/ (scripts/run.sh).

  python3 scripts/analyze.py [--all] [--arms cut06,base]

Each process starts from a depth-1 node 26:c (a stratum), bounds it exactly, then dives from its
depth-2 children. trace.txt holds, in order: a "dive i" line per dive, and per node the solver's
DNN trace line ("pid ...") followed by the patch's "dv" line with the node's weight (factor = the
product of 1/probability of every draw above it). A dive's estimate of a quantity over the stratum's
subtree below the depth-1 node is sum over its nodes of factor x quantity; with DIVE_FULL_DEPTH the
nodes below that depth all have their parent's factor (every child taken).

  stratum estimate  = exact top (the depth-1 node) + mean over the stratum's dives
  tree              = root (1 node) + sum of the 12 strata
  standard error    = sqrt(sum over strata of s_c^2 / n_c); the 80% band is a bootstrap within strata

Only processes that ended (cmd.txt has its exit line) count, unless --all (then every dive followed
by another dive line).
"""
import glob, math, os, random, re, statistics, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRATA = (0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19)
QUANT = ("nodes", "cpu", "it1", "itc", "dnn_nodes", "phases", "ph_fath", "ph_stall", "ph_limit", "ph_giveup")
LIMIT = 3400


def parse(path):
    """-> top (list of nodes), dives (list of lists of nodes), complete flags"""
    top, dives, cur, dnn = [], [], None, None
    for ln in open(path):
        if ln.startswith("dive "):
            cur = []
            dives.append(cur)
            dnn = None
        elif ln.startswith("pid "):
            dnn = ln
        elif ln.startswith("dv "):
            h, t = ln.split(" | ")
            hs = h.split()
            f = t.split()
            kv = dict(zip(f[0::2], f[1::2]))
            nd = dict(dive=int(hs[1]), depth=int(hs[3]), fixes=" ".join(hs[5:]), factor=float(kv["factor"]),
                      fate=kv["fate"], kids=int(kv["kids"]), imp=float(kv["imp"]), mode=int(kv["mode"]),
                      cpu=float(kv["cpu"]), qpb=float(kv["qpb"]), bound=float(kv["bound"]),
                      it1=0, itc=0, cut=0, cstat="", lb1=None, cutoff=None, pool=0, c1=0.0, cc=0.0, cbest=None)
            if dnn is not None:
                parts = dnn.split(" | ")
                if parts[0].split("fixes", 1)[1].split() == hs[5:]:
                    a = parts[1].split()
                    k1 = dict(zip(a[1::2], a[2::2]))
                    nd.update(it1=int(k1["iters"]), lb1=float(k1["lb"]), c1=float(k1["cpu"]), pool=int(k1["pool"]),
                              st1=k1["status"])
                    c = parts[2].split()
                    if c[1] == "1":
                        kc = dict(zip(c[2::2], c[3::2]))
                        nd.update(cut=1, itc=int(kc["iters"]), cstat=kc["status"], cc=float(kc["cpu"]),
                                  cbest=float(kc["best"]), stalled=int(kc.get("stalled", 0)),
                                  gaveup=int(kc.get("gaveup", 0)))
                    nd["cutoff"] = float(parts[3].split()[1])
                dnn = None
            (cur if nd["dive"] else top).append(nd)
    return top, dives


def q(nd):
    ph = nd["cut"]
    fath = ph and nd["fate"] == "C"
    lim = ph and not fath and nd["itc"] >= LIMIT
    gu = ph and not fath and not lim and nd.get("gaveup", 0)
    st = ph and not fath and not lim and not gu
    return dict(nodes=1, cpu=nd["cpu"], it1=nd["it1"], itc=nd["itc"], dnn_nodes=1 if nd["it1"] else 0,
                phases=ph, ph_fath=int(bool(fath)), ph_stall=int(bool(st)), ph_limit=int(bool(lim)),
                ph_giveup=int(bool(gu)))


def load(arm, use_all=False):
    """-> {c: dict(top=[nodes], dives=[[nodes]], runs=[names])}"""
    S = {}
    for d in sorted(glob.glob(f"{HERE}/runs/{arm}_s*_*")):
        m = re.match(rf"{arm}_s(\d+)_(\w+)$", os.path.basename(d))
        if not m or not os.path.exists(d + "/trace.txt"):
            continue
        done = "exit 0" in open(d + "/cmd.txt").read()
        if not done and not use_all:
            continue
        top, dives = parse(d + "/trace.txt")
        if not done:
            dives = dives[:-1]
        c = int(m.group(1))
        e = S.setdefault(c, dict(top=top, dives=[], runs=[]))
        e["dives"] += dives
        e["runs"].append(os.path.basename(d) + ("" if done else "*"))
    return S


def dive_tot(dv, key, depth=None):
    return sum(n["factor"] * q(n)[key] for n in dv if depth is None or n["depth"] == depth)


def top_tot(top, key, depth=None):
    return sum(q(n)[key] for n in top if depth is None or n["depth"] == depth)


def estimate(S, key, rng=random.Random(11), boots=4000):
    """tree total (without the root), se, 80% band, per-stratum (mean, se, n)"""
    per, tot, var = {}, 0.0, 0.0
    for c, e in S.items():
        x = [dive_tot(d, key) for d in e["dives"]]
        t = top_tot(e["top"], key)
        mu = statistics.fmean(x) if x else float("nan")
        se = statistics.stdev(x) / math.sqrt(len(x)) if len(x) > 1 else float("nan")
        per[c] = (t + mu, se, len(x), x, t)
        tot += t + mu
        var += se * se
    bs = sorted(sum(per[c][4] + statistics.fmean(rng.choices(per[c][3], k=len(per[c][3]))) for c in per)
                for _ in range(boots)) if all(per[c][2] for c in per) else [float("nan")] * boots
    return tot, math.sqrt(var), (bs[boots // 10], bs[boots * 9 // 10]), per


def fmt(x):
    if x != x:
        return "-"
    return f"{x:,.0f}" if abs(x) >= 100 else f"{x:.3g}"


def main():
    use_all = "--all" in sys.argv
    arms = ["cut06", "hand", "base"]
    if "--arms" in sys.argv:
        arms = sys.argv[sys.argv.index("--arms") + 1].split(",")
    R = {}
    for arm in arms:
        S = load(arm, use_all)
        if not S:
            continue
        R[arm] = S
        missing = [c for c in STRATA if c not in S or not S[c]["dives"]]
        nd = sum(len(e["dives"]) for e in S.values())
        print(f"\n## arm {arm}: {nd} dives in {len(S)} strata" + (f"; NO DIVES YET in strata {missing}: totals are partial" if missing else ""))
        print("\n| quantity | tree (root aside) | standard error | 80% band (bootstrap) | largest single dive's share |")
        print("|---|---|---|---|---|")
        for key in ("nodes", "cpu", "it1", "itc", "phases", "ph_fath", "ph_stall", "ph_limit", "ph_giveup"):
            tot, se, band, per = estimate(S, key)
            share = max((max(p[3]) / p[2] if p[3] else 0) for p in per.values()) / tot if tot else 0
            print(f"| {key} | {fmt(tot)} | {fmt(se)} | {fmt(band[0])} to {fmt(band[1])} | {share:.2f} |")
        print("\n| stratum | dives | nodes: mean (se) | CPU s: mean (se) | first-solve it | phase it | min-max dive CPU s | process CPU s per dive (unweighted) | runs |")
        print("|---|---|---|---|---|---|---|---|---|")
        _, _, _, pn = estimate(S, "nodes")
        _, _, _, pc = estimate(S, "cpu")
        _, _, _, p1 = estimate(S, "it1")
        _, _, _, p2 = estimate(S, "itc")
        for c in STRATA:
            if c not in S:
                continue
            raw = [sum(n["cpu"] for n in d) for d in S[c]["dives"]]
            x = pc[c][3]
            print(f"| 26:{c} | {pn[c][2]} | {fmt(pn[c][0])} ({fmt(pn[c][1])}) | {fmt(pc[c][0])} ({fmt(pc[c][1])}) | {fmt(p1[c][0])} | {fmt(p2[c][0])} | "
                  f"{fmt(min(x)) if x else '-'} to {fmt(max(x)) if x else '-'} | {fmt(statistics.fmean(raw)) if raw else '-'} | {' '.join(S[c]['runs'])} |")
        print("\n| depth | nodes | CPU s | share of CPU | first-solve it | phase it | phases | fathomed | stalled | at the limit | gave up | nodes seen in dives (unweighted) | phases seen | of them fathomed / stalled / limit |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        depths = sorted({n["depth"] for e in S.values() for d in e["dives"] for n in d} | {1})
        tot_cpu = estimate(S, "cpu")[0]
        for dp in depths:
            row = {}
            for key in QUANT:
                row[key] = sum(top_tot(e["top"], key, dp) + (statistics.fmean(dive_tot(d, key, dp) for d in e["dives"]) if e["dives"] else 0)
                               for e in S.values())
            seen = [n for e in S.values() for d in e["dives"] for n in d if n["depth"] == dp] + [n for e in S.values() for n in e["top"] if n["depth"] == dp]
            ph = [q(n) for n in seen if n["cut"]]
            print(f"| {dp} | {fmt(row['nodes'])} | {fmt(row['cpu'])} | {row['cpu'] / tot_cpu:.1%} | {fmt(row['it1'])} | {fmt(row['itc'])} | {fmt(row['phases'])} | "
                  f"{fmt(row['ph_fath'])} | {fmt(row['ph_stall'])} | {fmt(row['ph_limit'])} | {fmt(row['ph_giveup'])} | {len(seen)} | {len(ph)} | "
                  f"{sum(p['ph_fath'] for p in ph)} / {sum(p['ph_stall'] for p in ph)} / {sum(p['ph_limit'] for p in ph)} |")
        # seconds per iteration under load, first solves and phases by m
        print("\n| depth (m) | first solves seen | s per first-solve it | phases seen | s per phase it |")
        print("|---|---|---|---|---|")
        for dp in depths:
            seen = [n for e in S.values() for d in e["dives"] for n in d if n["depth"] == dp] + [n for e in S.values() for n in e["top"] if n["depth"] == dp]
            a = [n for n in seen if n["it1"] >= 50]
            b = [n for n in seen if n["cut"] and n["itc"] >= 50]
            if a:
                print(f"| {dp} ({40 - dp}) | {len(a)} | {sum(n['c1'] for n in a) / sum(n['it1'] for n in a):.4f} | {len(b)} | "
                      + (f"{sum(n['cc'] for n in b) / sum(n['itc'] for n in b):.4f}" if b else "-") + " |")
    if "cut06" in R and "base" in R:
        print("\n## cut06 / base on matched strata (this box, same estimator)\n")
        print("| stratum | base dives | cut dives | base nodes | cut nodes | nodes ratio | base CPU s | cut CPU s | CPU ratio cut / base |")
        print("|---|---|---|---|---|---|---|---|---|")
        en = {a: estimate(R[a], "nodes")[3] for a in R}
        ec = {a: estimate(R[a], "cpu")[3] for a in R}
        sb = sc = nb = nc = 0
        for c in STRATA:
            if c in en["base"] and c in en["cut06"] and en["base"][c][2] and en["cut06"][c][2]:
                b, k = ec["base"][c][0], ec["cut06"][c][0]
                sb += b; sc += k; nb += en["base"][c][0]; nc += en["cut06"][c][0]
                print(f"| 26:{c} | {en['base'][c][2]} | {en['cut06'][c][2]} | {fmt(en['base'][c][0])} | {fmt(en['cut06'][c][0])} | {en['cut06'][c][0] / en['base'][c][0]:.4f} | {fmt(b)} | {fmt(k)} | {k / b:.2f} |")
        if sb:
            print(f"| all matched | | | {fmt(nb)} | {fmt(nc)} | {nc / nb:.4f} | {fmt(sb)} | {fmt(sc)} | {sc / sb:.2f} |")


if __name__ == "__main__":
    main()
