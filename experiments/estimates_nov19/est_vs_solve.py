#!/usr/bin/env python3
"""A dives estimate against a solve of the same tree: every sampled open node (stratum) of the estimate's dives
run, its estimate (exact part + its dives' mean) against its true subtree in the solve's node trace, weighted
by 1/pi and summed by the open node's depth; and the stratum node's own bound in each.

    python3 experiments/estimates_nov19/est_vs_solve.py infra/phase1/runs/sko42-est-261121b infra/phase1/runs/sko42-261120b
"""
import collections
import gzip
import json
import os
import statistics as st
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "experiments", "estimates", "sko42", "zone_nov19"))
from common import load  # noqa: E402


def main():
    est_dir, solve_dir = sys.argv[1:3]
    T = load(os.path.join(solve_dir, "nodes.txt.gz"))
    n, cpu, size = T["nodes"], T["cpu"], T["size"]
    strata, dives = {}, collections.defaultdict(list)
    for line in gzip.open(os.path.join(est_dir, "dives", "dives.jsonl.gz"), "rt"):
        x = json.loads(line)
        if x["type"] == "stratum":
            strata[x["key"]] = x
        elif x["type"] == "dive":
            dives[x["key"]].append(x)
    own = [(s["exact_cpu_by_depth"].get(str(s["depth"]), 0.0), n[k]["cpu"], n[k]["out"]) for k, s in strata.items() if k in n]
    print(f"{len(strata):,} strata, {sum(1 for k in strata if k not in n)} not in the solve's tree")
    print(f"their own bound: {sum(a for a, _, _ in own) / 3600:.1f} CPU h in the estimate, {sum(b for _, b, _ in own) / 3600:.1f} "
          f"in the solve ({sum(a for a, _, _ in own) / sum(b for _, b, _ in own):.2f}); by the solve's outcome: " +
          ", ".join(f"{o} {sum(a for a, _, x in own if x == o) / sum(b for _, b, x in own if x == o):.2f}"
                    for o in "DCB" if any(x == o for _, _, x in own)))
    agg = collections.defaultdict(lambda: [0.0] * 6)
    for k, s in strata.items():
        if k not in n:
            continue
        ds = dives.get(k, [])
        w = 1.0 / (float(s.get("pi", 1)) or 1.0)
        a = agg[s["depth"]]
        a[0] += w * (s["exact_cpu"] + (st.fmean(d["cpu"] for d in ds) if ds else 0.0))
        a[1] += w * cpu[k]
        a[2] += w * (s["exact_nodes"] + (st.fmean(d["nodes"] for d in ds) if ds else 0.0))
        a[3] += w * size[k]
        a[4] += 1
        a[5] += bool(ds)
    print("\n| open nodes at depth | strata (dived) | estimate CPU h | true CPU h | ratio | estimate nodes | true nodes | ratio |")
    print("|---|---|---|---|---|---|---|---|")
    for d in sorted(agg):
        a = agg[d]
        print(f"| {d} | {int(a[4]):,} ({int(a[5]):,}) | {a[0] / 3600:,.1f} | {a[1] / 3600:,.1f} | {a[0] / a[1]:.2f} | {a[2]:,.0f} | {a[3]:,.0f} | {a[2] / a[3]:.2f} |")
    t = [sum(a[i] for a in agg.values()) for i in range(4)]
    print(f"| all | | {t[0] / 3600:,.1f} | {t[1] / 3600:,.1f} | {t[0] / t[1]:.2f} | {t[2]:,.0f} | {t[3]:,.0f} | {t[2] / t[3]:.2f} |")


if __name__ == "__main__":
    main()
