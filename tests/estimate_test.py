#!/usr/bin/env python3
"""Check the tree-size estimator (bnbdrv -e, latest/estimate.C; thesis chapter 4).

  python3 tests/estimate_test.py [--solver PATH]

nug15 is solved once to count its tree. Then the estimator runs 3 x 400 dives
each way a dive can weight children: by their predicted bounds (the default)
and by their own QPB bounds (QAP_EST_CHILD_QPB=1). Knuth's estimate is unbiased
whatever the weights, so both pooled means must come within 15% of the count,
and so must their CPU estimates of the run's CPU (a weighting that changed the
search itself, or counted its own bounds as the tree's time, fails here: the
first try at QAP_EST_CHILD_QPB read 0.16 of nug20's nodes). With the setting the
first draw must say it scored its nodes.

Then the estimator resumed from a state file (-e -r, as dives from a saved frontier are run): its
first step must bound the state's nodes exactly as a plain resume of the same state does (the
same bounds, strategy and branching line at every depth-1 node). Before Nov 16 it bounded them
with root bound 0, so with the deepest level's strategy. Last, QAP_EST_SURV's dives from the depth-2
frontier, three seeds of 1,000 dives each: survivor dives (b = 0 and 0.4), Knuth's in strata (k1.5),
and Knuth's with 30% of the strata sampled by QPB, each within 5% of the count on the three's mean.
About 40 seconds in all.
"""

import argparse
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCE, UB = "nug15", 1150
SEEDS, DIVES, TOLERANCE = (1, 2, 3), 400, 0.15


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    a = ap.parse_args()
    lat = os.path.join(ROOT, "latest")
    common = ["-p", os.path.join(lat, "PARAM", f"{INSTANCE}.par"), "-s", os.path.join(lat, "SYM", f"{INSTANCE}.sym"),
              "-u", str(UB), os.path.join(ROOT, "qaplib", f"{INSTANCE}.dat")]
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")
    full = subprocess.run([a.solver] + common, cwd=lat, env=env, capture_output=True, text=True, check=True).stdout
    nodes = int(re.search(r"Nodes:\s*(\d+)", full).group(1))
    cpu = float(re.search(r"Total Time:\s*(\S+)", full).group(1))
    print(f"{INSTANCE}: {nodes:,} nodes, {cpu:.2f} CPU s")
    failed = 0
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "estimate.par"), "w") as f:   # start level, dives, print every, print all every, q, ignore depth
            f.write(f"2 {DIVES} 1 100000 2 1000\n")
        for arm in ("0", "1"):
            n_est, t_est = [], []
            for seed in SEEDS:
                hist = os.path.join(d, f"dives_{arm}_{seed}.txt")
                out = subprocess.run([a.solver, "-e"] + common, cwd=d, capture_output=True, text=True, check=True,
                                     env=dict(env, QAP_EST_CHILD_QPB=arm, QAP_EST_SEED=str(seed), QAP_EST_HIST=hist)).stdout
                if arm == "1" and "QAP_EST_CHILD_QPB: the first draw's" not in out:
                    print(f"FAIL: QAP_EST_CHILD_QPB=1, seed {seed}: the first draw did not score its nodes")
                    failed += 1
                rows = [l.split() for l in open(hist) if l.strip()]
                if len(rows) != DIVES:
                    print(f"FAIL: QAP_EST_CHILD_QPB={arm}, seed {seed}: {len(rows)} dives, not {DIVES}")
                    failed += 1
                n_est += [float(r[1]) for r in rows]
                t_est += [float(r[2]) for r in rows]
            rn, rt = statistics.fmean(n_est) / nodes, statistics.fmean(t_est) / cpu
            ok = abs(rn - 1) <= TOLERANCE and abs(rt - 1) <= TOLERANCE
            failed += not ok
            print(f"{'ok  ' if ok else 'FAIL'} QAP_EST_CHILD_QPB={arm}: {len(n_est)} dives, nodes {rn:.3f}, CPU {rt:.3f} of the run's")
        failed += resume_check(a.solver, common, env, d)
        failed += survivor_check(a.solver, common, env, d, nodes, full)
        failed += spill_check(a.solver, env, d)
    print("estimator: all checks passed" if not failed else f"estimator: {failed} checks failed")
    return 1 if failed else 0


def spill_check(solver, env, d):
    """the strata past QAP_EST_WARM_KEEP (estimate.C, spill_save): nug24 under the nug30 recipe (the DNN
    at every depth), stopped at 12 CPU s with its open nodes saved, then Knuth's dives from that state,
    one seed three ways: every stratum's state kept in memory, every one spilled to a file and read
    back around its dives, and dropped as before Nov 19 (QAP_EST_SPILL_MB=0: the first step cold,
    without the parent's pool). The spilled run must give the in-memory run's node estimate and a
    first-step CPU within 30% of it, and say it spilled; the cold one costs more (tho40-est2-261119:
    2.1-4.6x). About a minute. (1 failure or 0)"""
    lat = os.path.join(ROOT, "latest")
    par, dat = os.path.join(ROOT, "infra", "phase1", "recipes", "nug30.par"), os.path.join(ROOT, "qaplib", "nug24.dat")
    st = os.path.join(d, "spill.state")
    r = subprocess.run([solver, "-p", par, "-u", "3488", "-o", st, "-t", "12", dat], cwd=lat, env=env,
                       capture_output=True, text=True)
    if r.returncode not in (0, 3) or not os.path.exists(st):
        print(f"FAIL spill: the nug24 state was not written (exit {r.returncode})")
        return 1
    res = {}
    for arm, extra in (("mem", {"QAP_EST_WARM_KEEP": "1000"}), ("spill", {"QAP_EST_WARM_KEEP": "0"}),
                       ("cold", {"QAP_EST_WARM_KEEP": "0", "QAP_EST_SPILL_MB": "0"})):
        dv = os.path.join(d, f"spill_{arm}.jsonl")
        out = subprocess.run([solver, "-p", par, "-u", "3488", "-e", "-r", st, dat], cwd=lat, capture_output=True, text=True,
                             env=dict(env, QAP_EST_SURV="k1.5", QAP_EST_BUDGET="1", QAP_EST_MIN_DIVES="4",
                                      QAP_EST_DIVES=dv, QAP_EST_SEED="7", **extra)).stdout
        lines = [json.loads(l) for l in open(dv)]
        strata = {x["key"]: x for x in lines if x["type"] == "stratum"}
        dives = [x for x in lines if x["type"] == "dive"]
        first = []
        for x in dives:
            k = str(strata[x["key"]]["depth"] + 1)
            if x["by_depth"].get(k):
                first.append(x["cpu_by_depth"].get(k, 0.0) / x["by_depth"][k])
        nodes = sum(s["exact_nodes"] for s in strata.values()) + \
            sum(statistics.fmean(x["nodes"] for x in dives if x["key"] == k) for k in {x["key"] for x in dives})
        res[arm] = (nodes, statistics.fmean(first) if first else 0.0, sum(s.get("warm") == "s" for s in strata.values()),
                    "strata's states spilled" in out)
    (mn, mc, _, _), (sn, sc, ss, said), (cn, cc, _, _) = res["mem"], res["spill"], res["cold"]
    ok = abs(sn - mn) <= 0.02 * mn and mc > 0 and abs(sc / mc - 1) <= 0.3 and ss > 0 and said and not \
        any(f.startswith("qap-est-spill-") for f in os.listdir(d))
    print(f"{'ok  ' if ok else 'FAIL'} spill: Knuth's dives from a saved nug24 frontier (nug30 recipe): nodes "
          f"{mn:.1f} in memory, {sn:.1f} spilled ({ss} strata), {cn:.1f} cold; the first child's CPU {mc:.2f} s, "
          f"{sc:.2f} s, {cc:.2f} s (x{cc / mc if mc else 0:.1f} cold); spill files removed")
    return 0 if ok else 1


def resume_check(solver, common, env, d):
    """-e -r: step 1 bounds a state's nodes as a plain -r does (1 failure or 0)"""
    flags, inst = common[:-1], common[-1]
    s0, s1 = os.path.join(d, "s0.state"), os.path.join(d, "s1.state")
    plain_t, est_t = os.path.join(d, "plain.nodes"), os.path.join(d, "est.nodes")
    run = lambda extra, **e: subprocess.run([solver] + flags + extra + [inst], cwd=d, capture_output=True, text=True,
                                            check=True, env=dict(env, **e))
    run(["-d", "0", "-o", s0])                                       # the root, its children saved
    run(["-r", s0, "-d", "1", "-o", s1], QAP_NODE_TRACE=plain_t)     # those children, as a resume bounds them
    with open(os.path.join(d, "estimate.par"), "w") as f:            # dives start at level 2: step 1 bounds level 1
        f.write("2 20 1 100000 2 1000\n")
    run(["-e", "-r", s0], QAP_NODE_TRACE=est_t, QAP_EST_SEED="1", QAP_EST_HIST=os.path.join(d, "resume_dives.txt"))
    level1 = lambda path: sorted(" ".join(l.split()[:-1]) for l in open(path) if l.split()[:2] == ["node", "1"])
    a, b = level1(plain_t), level1(est_t)   # the last field is the node's CPU time
    ok = bool(a) and a == b
    print(f"{'ok  ' if ok else 'FAIL'} -e -r: step 1 bounds the state's {len(a)} depth-1 nodes as a plain resume "
          f"does ({sum(x == y for x, y in zip(a, b))} of {len(a)} lines identical)")
    return 0 if ok else 1


def survivor_check(solver, common, env, d, nodes, full):
    """QAP_EST_SURV (survivor dives, estimate.C) from the depth-2 frontier, 3 seeds x 1,000 dives each,
    b = 0 and 0.4: the frontier's nodes plus each stratum's exact part and dive mean must come within 5%
    of the count (on nug15 every seed and weight is within 2%, Nov 18); failures counted. Since Nov 21 the
    lines also carry fathoms, locations considered and eliminated by depth (the RESULTS table's %fathom and
    %c.elim in coord/dives.py's results_estimate.txt): the estimated fathoms over the estimated nodes must
    come within 0.03 of the run's share (an unweighted count reads far below), and with survivor dives
    %c.elim within 0.05 of the run's at every level below the frontier with 5% of the tree's nodes. Knuth's
    dives are printed, not held to it: ten k1.5 seeds read the share at 0.787 on average against 0.786, but
    %c.elim at levels 6-7 from 0.58 to 0.88 seed by seed (0.63 and 0.66 in the run), Nov 21"""
    import json
    true_f = int(re.search(r"Nodes:\s*\d+\s*/\s*(\d+)", full).group(1))
    lv = {}
    table = full.split("%c.elim", 1)[1].split("\n\n", 1)[0]     # the RESULTS table, not the DNN one below it
    for m in re.finditer(r"^\s*(\d+)\s+(\d+)\t([0-9.]+)\t(N/A|[0-9.]+)\t", table, re.M):
        lv[int(m.group(1))] = (int(m.group(2)), None if m.group(4) == "N/A" else float(m.group(4)))
    flags, inst = common[:-1], common[-1]
    fr = os.path.join(d, "f2.state")
    top = subprocess.run([solver] + flags + ["-d", "2", "-o", fr, inst], cwd=d, capture_output=True, text=True,
                         check=True, env=env).stdout
    top_nodes = int(re.search(r"Nodes:\s*(\d+)", top).group(1))
    top_f = int(re.search(r"Nodes:\s*\d+\s*/\s*(\d+)", top).group(1))
    with open(os.path.join(d, "estimate.par"), "w") as f:
        f.write("2 1000 100 100000 1.5 1000\n")
    bad = 0
    for b, extra in (("0", {}), ("0.4", {}), ("k1.5", {}),
                     ("k1.5", {"QAP_EST_SAMPLE": "0.3", "QAP_EST_SAMPLE_W": "qpb"})):
        ests, fests, locs, elim = [], [], {}, {}
        for seed in (1, 2, 3):
            jl = os.path.join(d, f"surv_{b}_{extra.get('QAP_EST_SAMPLE', '1')}_{seed}.jsonl")
            if os.path.exists(jl):
                os.remove(jl)                    # the estimator appends
            subprocess.run([solver, "-e", "-r", fr] + flags + [inst], cwd=d, capture_output=True, text=True, check=True,
                           env=dict(env, QAP_EST_SURV=b, QAP_EST_DIVES=jl, QAP_EST_SEED=str(seed), **extra))
            strata, dives, dl = {}, {}, {}
            for line in open(jl):
                x = json.loads(line)
                if x["type"] == "stratum":
                    strata[x["key"]] = x
                elif x["type"] == "dive":
                    dives.setdefault(x["key"], []).append(x["nodes"])
                    dl.setdefault(x["key"], []).append(x)
            pi = {k: s.get("pi", 1.0) for k, s in strata.items()}
            ests.append(top_nodes + sum(s["exact_nodes"] / pi[k] for k, s in strata.items()) +
                        sum(statistics.fmean(v) / pi[k] for k, v in dives.items()))
            # fathoms, locations and eliminated by depth, weighted as the nodes are (coord/dives.py, combine)
            byd = lambda m: {int(k): float(v) for k, v in (m or {}).items()}
            f = top_f
            for k, st in strata.items():
                f += sum(byd(st.get("exact_fathom_by_depth")).values()) / pi[k]
                for acc, name in ((locs, "exact_locs_by_depth"), (elim, "exact_elim_by_depth")):
                    for dd, v in byd(st.get(name)).items():
                        acc[dd] = acc.get(dd, 0.0) + v / pi[k] / 3
            for k, xs in dl.items():
                for x in xs:
                    f += sum(byd(x.get("fathom_by_depth")).values()) / len(xs) / pi[k]
                    for acc, name in ((locs, "locs_by_depth"), (elim, "elim_by_depth")):
                        for dd, v in byd(x.get(name)).items():
                            acc[dd] = acc.get(dd, 0.0) + v / len(xs) / pi[k] / 3
            fests.append(f)
        # 1,000 dives a seed: one seed's estimate moves by 5-9% (Knuth's tails, a sample of the strata),
        # so the mean of the three is held to 5%
        sampled = bool(extra)
        ok = abs(statistics.fmean(ests) / nodes - 1) <= 0.05
        bad += not ok
        what = f"QAP_EST_SURV={b}" + (f" QAP_EST_SAMPLE={extra['QAP_EST_SAMPLE']} ({extra['QAP_EST_SAMPLE_W']})"
                                      if sampled else "")
        print(f"{'ok  ' if ok else 'FAIL'} {what}: dives from the depth-2 frontier, 3 seeds: "
              f"{', '.join(f'{e / nodes:.3f}' for e in ests)} of the nodes, mean {statistics.fmean(ests) / nodes:.3f}")
        share = sum(fests) / sum(ests)
        okf = abs(share - true_f / nodes) <= 0.03
        rows = [(dd, elim.get(dd, 0.0) / locs[dd], lv[dd][1]) for dd in sorted(locs)
                if dd >= 2 and locs[dd] > 0 and dd in lv and lv[dd][1] is not None and lv[dd][0] >= 0.05 * nodes]
        oke = bool(rows) and (b.startswith("k") or all(abs(e - t) <= 0.05 for _, e, t in rows))
        bad += (not okf) + (not oke)
        print(f"{'ok  ' if okf and oke else 'FAIL'} {what}: fathomed share {share:.3f} (the run's {true_f / nodes:.3f}); "
              f"%c.elim by level, estimated / run{' (not checked)' if b.startswith('k') else ''}: "
              + ", ".join(f"{dd}: {e:.3f}/{t:.3f}" for dd, e, t in rows))
    return bad


if __name__ == "__main__":
    sys.exit(main())
