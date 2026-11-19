#!/usr/bin/env python3
"""The estimate's decision point (coord/dives.py, decision()) replayed on a collected dives run.

    python3 experiments/estimates_nov19/decide_sim.py infra/phase1/runs/tai50b-est-261119 [--frac 0.15]

From the run's dives export (dives/units.tsv, dives/dives.jsonl.gz: `collect.sh`), the units are
taken in their merge order; after each, the part-way estimate of the rest (coord/dives.py's ratio
imputation, depth by depth) and the estimator CPU the dives run will have spent (the merged units'
CPU, the others at the merged units' mean). The running units' time so far is not in the export,
so the projection is the low one: the replay decides no earlier than the live check would.
A unit's open nodes by depth are not in the export either: they are its strata's, each weighted
1/pi, scaled to its open count (the sample line), which is exact where every open node is a stratum.
"""
import argparse
import csv
import gzip
import json
import math
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from coord.dives import DECIDE, DECIDE_MIN, _parts, read_lines, rest_from, unit_rest  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--frac", type=float, default=DECIDE)
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(a.run_dir, "dives", "units.tsv")), delimiter="\t"))
    lines = gzip.open(os.path.join(a.run_dir, "dives", "dives.jsonl.gz"), "rt").read().splitlines()
    units, keys, strata, dives, open_by = [], {}, {}, {}, {}
    for r in rows:
        if r["has_dives"] != "yes":
            continue
        text = "\n".join(lines[int(r["line_from"]):int(r["line_to"])])
        samples = []
        s1, d1 = read_lines(text, samples)
        strata.update(s1)
        for k, v in d1.items():
            dives.setdefault(k, []).extend(v)
        keys[r["uid"]] = list(s1)
        w = {}
        for st in s1.values():
            w[int(st["depth"])] = w.get(int(st["depth"]), 0.0) + 1.0 / (float(st.get("pi", 1)) or 1.0)
        n_open = sum(x.get("open", 0) for x in samples) or int(r["nodes"])
        tw = sum(w.values())
        if tw:
            ob = {d: n_open * v / tw for d, v in w.items()}
        else:       # no stratum sampled: its open nodes at its shallowest depth
            ob = {int(r["depth_min"]): float(n_open)}
        open_by[r["uid"]] = ob
        units.append({"uid": r["uid"], "seq": int(r["merged_seq"]), "cpu": float(r["cpu_s"] or 0),
                      "by_depth": {str(d): v for d, v in ob.items()}})
    units.sort(key=lambda u: u["seq"])
    est = json.load(open(os.path.join(a.run_dir, "estimate.json")))
    parts = _parts(strata, dives)
    n = len(units)
    need = max(math.ceil(DECIDE_MIN * n), min(10, math.ceil(n / 2)))
    accs = dict(zip([u["uid"] for u in units], unit_rest(units, keys, strata, parts)))
    total_cpu = sum(u["cpu"] for u in units)
    print(f"{os.path.basename(a.run_dir.rstrip('/'))}: {n} dives units, {total_cpu / 3600:,.1f} estimator CPU h; "
          f"the estimate's rest {est['rest']['cpu'] / 3600:,.1f} CPU h, so the whole dives run cost "
          f"{total_cpu / est['rest']['cpu']:.1%} of it")
    # (1) the first units to merge, as a check on merged units would see them: biased to the small ones
    for k in (need,):
        m = units[:k]
        pend = {}
        for u in units[k:]:
            for d, v in open_by[u["uid"]].items():
                pend[d] = pend.get(d, 0.0) + v
        _, rc, *_ = rest_from([accs[u["uid"]] for u in m], pend)
        spent = sum(u["cpu"] for u in m)
        print(f"  the first {k} to merge: rest {rc / 3600:,.1f} CPU h, projected estimator {spent / k * n / 3600:,.1f} h: "
              f"share {spent / k * n / rc:.1%}")
    # (2) a pilot: a random {need} units fixed in advance, all of them merged before the decision
    rng = random.Random(a.seed)
    shares, costs, rests = [], [], []
    for _ in range(a.reps):
        pilot = rng.sample(units, need)
        ids = {u["uid"] for u in pilot}
        others = {}
        for u in units:
            if u["uid"] not in ids:
                for d, v in open_by[u["uid"]].items():
                    others[d] = others.get(d, 0.0) + v
        _, rc, *_ = rest_from([accs[u["uid"]] for u in pilot], others)
        pc = sum(u["cpu"] for u in pilot)
        proj = pc / need * n
        shares.append(proj / rc if rc > 0 else float("inf"))
        costs.append(pc)
        rests.append(rc)
    shares.sort()
    q = lambda xs, f: sorted(xs)[min(int(f * len(xs)), len(xs) - 1)]
    solve = sum(1 for x in shares if x > a.frac) / len(shares)
    print(f"  a random pilot of {need} units, {a.reps} draws: decides to solve in {solve:.0%} of them at {a.frac:.0%}; "
          f"share 10-50-90%: {q(shares, .1):.1%} {q(shares, .5):.1%} {q(shares, .9):.1%}; "
          f"the pilot's rest 10-50-90%: {q(rests, .1) / 3600:,.1f} {q(rests, .5) / 3600:,.1f} {q(rests, .9) / 3600:,.1f} CPU h; "
          f"its CPU 10-50-90%: {q(costs, .1) / 3600:,.1f} {q(costs, .5) / 3600:,.1f} {q(costs, .9) / 3600:,.1f} h")


if __name__ == "__main__":
    main()
