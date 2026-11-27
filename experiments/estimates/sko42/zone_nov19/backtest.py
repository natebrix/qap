#!/usr/bin/env python3
"""Backtest rest.py's imputation: cut a node trace at a share of its lines (the trace is in merge order, about
time order), impute the rest of each depth-1 subtree from the cut, and compare with the depth-1 subtrees that
are finished, or nearly (imputed rest under 5% of their CPU), by the end of the whole trace. Favourable to the
imputation: the subtrees that finish are the ones without a long tail.

    python3 experiments/estimates/sko42/zone_nov19/backtest.py [TRACE ...]
"""
import sys

from common import ARMS, build, read_lines
from rest import impute


def main():
    for path in sys.argv[1:] or [ARMS["z25"], ARMS["z10"]]:
        lines = read_lines(path)
        final, _, _ = impute(build(lines))
        truth = {k: g["done"] + g["sib"] for k, g in final.items() if k and g["sib"] < 0.05 * g["done"]}
        tot = sum(truth.values())
        print(f"{path.split('/')[-2]}: {len(truth)} depth-1 subtrees finished or nearly by the end, "
              f"{tot / 3600:,.0f} CPU h")
        for share in (0.3, 0.5, 0.7, 0.85):
            e, _, _ = impute(build(lines[:int(share * len(lines))]))
            done = sum(e[k]["done"] for k in truth)
            rd, rs = sum(e[k]["dep"] for k in truth), sum(e[k]["sib"] for k in truth)
            print(f"  cut at {share:.0%}: done {done / 3600:,.0f} CPU h; of the true rest, the depth means give "
                  f"{rd / (tot - done):.2f}, the siblings' {rs / (tot - done):.2f}")


if __name__ == "__main__":
    main()
