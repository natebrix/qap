#!/usr/bin/env python3
"""The subtrees both arms closed, paired by their root's fixes: nodes and CPU in each arm, by the root's depth
(a maximal pair: closed in both, its parent not), with a bootstrap over the pairs for the CPU ratio.

    python3 experiments/estimates/sko42/zone_nov19/closed_pairs.py [A_TRACE B_TRACE]
"""
import random
import sys

from common import ARMS, load, parent


def main():
    pa, pb = sys.argv[1:3] if len(sys.argv) > 2 else (ARMS["z10"], ARMS["z25"])
    A, B = load(pa), load(pb)
    both = {k for k in A["nodes"] if k in B["nodes"] and A["closed"][k] and B["closed"][k]}
    top = [k for k in both if k and parent(k) not in both]
    print(f"closed in both: {len(both):,} nodes; maximal such subtrees: {len(top):,}\n")
    print("| root depth | subtrees | A nodes | B nodes | ratio | A CPU h | B CPU h | B/A CPU (5-95%) |")
    print("|---|---|---|---|---|---|---|---|")
    rng = random.Random(1)
    groups = [(d, [k for k in top if A["nodes"][k]["depth"] == d]) for d in sorted({A["nodes"][k]["depth"] for k in top})]
    for d, ks in groups + [("all", top)]:
        na, nb = sum(A["size"][k] for k in ks), sum(B["size"][k] for k in ks)
        ca, cb = sum(A["cpu"][k] for k in ks), sum(B["cpu"][k] for k in ks)
        boots = sorted(sum(B["cpu"][k] for k in s) / sum(A["cpu"][k] for k in s)
                       for s in ([ks[rng.randrange(len(ks))] for _ in ks] for _ in range(2000)))
        print(f"| {d} | {len(ks):,} | {na:,} | {nb:,} | {nb / na:.2f} | {ca / 3600:,.1f} | {cb / 3600:,.1f} | "
              f"{cb / ca:.2f} ({boots[100]:.2f}-{boots[1900]:.2f}) |")
    ca = sum(A["cpu"][k] for k in top)
    cb = sum(B["cpu"][k] for k in top)
    print(f"\nshare of each arm's CPU in these subtrees: A {ca / A['cpu']['']:.1%}, B {cb / B['cpu']['']:.1%}")
    same = [k for k in top if A["size"][k] == B["size"][k] and A["nodes"][k]["out"] == "B"]
    print(f"branched roots whose subtrees have the same node count in both: {len(same):,}; B/A CPU on them "
          f"{sum(B['cpu'][k] for k in same) / sum(A['cpu'][k] for k in same):.2f}")


if __name__ == "__main__":
    main()
