#!/usr/bin/env python3
"""The cost and precision of a sampled A/B's draw (README.md, "Tuning on paired subtrees"), on a
collected dives run's strata.

    python3 experiments/estimates_nov19/sample_design.py infra/phase1/runs/tho40-est2-261119 [--depths 3,4]

The population at each depth: the dives run's strata there whose first step was warm (the first 8 of
their unit: see "Knuth dives from cold strata"), each sized by its estimate (exact part + dives' mean).
For a draw with pi ~ size^a (a = 0 uniform, 0.5, 1; coord/sample.py's inclusion, mix 0.1 for a > 0) and n
nodes: the expected CPU of the sample (sum of pi x size), the Horvitz-Thompson total's relative standard
error, and a ratio's between arms whose effect varies by a share sigma from node to node (about sigma
times the total's).
"""
import argparse
import csv
import gzip
import json
import math
import os
import statistics as stx
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from coord.sample import inclusion  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--depths", default="3,4")
    a = ap.parse_args()
    depths = [int(x) for x in a.depths.split(",")]
    rows = list(csv.DictReader(open(os.path.join(a.run_dir, "dives", "units.tsv")), delimiter="\t"))
    lines = gzip.open(os.path.join(a.run_dir, "dives", "dives.jsonl.gz"), "rt").read().splitlines()
    sizes = {d: [] for d in depths}
    for r in rows:
        if r["has_dives"] != "yes":
            continue
        n, cur = 0, {}
        for l in lines[int(r["line_from"]):int(r["line_to"])]:
            x = json.loads(l)
            if x["type"] == "stratum":
                cur[x["key"]] = dict(depth=x["depth"], warm=n < 8, ec=x["exact_cpu"], dives=[])
                n += 1
            elif x["type"] == "dive":
                cur[x["key"]]["dives"].append(x["cpu"])
        for s in cur.values():
            if s["depth"] in sizes and s["warm"]:
                sizes[s["depth"]].append(s["ec"] + (stx.fmean(s["dives"]) if s["dives"] else 0.0))
    for d in depths:
        v = sizes[d]
        T = sum(v)
        print(f"depth {d}: {len(v)} warm strata, {T / 3600:,.0f} CPU h; median {stx.median(v) / 3600:.3f} h, "
              f"mean {T / len(v) / 3600:.2f} h")
        for power in (0.0, 0.5, 1.0):
            for n in (25, 50, 100, 150, 200, 300):
                if n > len(v):
                    continue
                pi = inclusion([s ** power for s in v], n, 0.0 if power == 0 else 0.1)
                cost = sum(p * s for p, s in zip(pi, v)) / 3600
                rel = math.sqrt(sum((1 - p) / p * s * s for p, s in zip(pi, v) if p > 0)) / T
                print(f"  pi ~ size^{power:.1f}, n {n:3d}: {cost:6,.0f} CPU h an arm; total {rel:5.1%}; "
                      f"ratio at sigma 20% {0.2 * rel:5.1%}, 40% {0.4 * rel:5.1%}")


if __name__ == "__main__":
    main()
