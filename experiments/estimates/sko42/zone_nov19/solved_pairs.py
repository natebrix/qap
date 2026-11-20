#!/usr/bin/env python3
"""sko42 solved (sko42-261120b) paired with the earlier runs on the subtrees both closed, leaving out any
subtree that holds, or lies under, a failed pooled solve in either (README.md, "Failed pooled solves").

    python3 experiments/estimates/sko42/zone_nov19/solved_pairs.py
"""
import os
import random

from common import RUNS, load, parent

def failed_roots(T):
    n = T["nodes"]
    lb = {k: max([v for v in (x["qpb"], x["dnn"], x["cut"]) if v is not None] or [float("-inf")]) for k, x in n.items()}
    f = [k for k, x in n.items() if k and x["dnn"] is not None and parent(k) in n and x["out"] != "D" and x["dnn"] < lb[parent(k)] - 158.14]
    return [k for k in f if not any(k.startswith(g + ",") for g in f if g != k)]
def pair(an, bn):
    A = load(os.path.join(RUNS, an, "nodes.txt.gz")); B = load(os.path.join(RUNS, bn, "nodes.txt.gz"))
    fa, fb = failed_roots(A), failed_roots(B)
    bad = fa + fb
    touches = lambda k: any(k == f or k.startswith(f + ",") or f.startswith(k + ",") for f in bad)
    both = {k for k in A["nodes"] if k in B["nodes"] and A["closed"][k] and B["closed"][k] and not touches(k)}
    top = [k for k in both if k and parent(k) not in both]
    na, nb = sum(A["size"][k] for k in top), sum(B["size"][k] for k in top)
    ca, cb = sum(A["cpu"][k] for k in top), sum(B["cpu"][k] for k in top)
    rng = random.Random(1)
    boots = sorted(sum(B["cpu"][k] for k in s) / sum(A["cpu"][k] for k in s) for s in ([top[rng.randrange(len(top))] for _ in top] for _ in range(1000)))
    same = [k for k in top if A["size"][k] == B["size"][k]]
    print(f"{bn} against {an}: failed roots {len(fa)} / {len(fb)} left out; {len(top):,} paired subtrees, {ca/A['cpu']['']:.0%} of {an}'s CPU and {cb/B['cpu']['']:.0%} of {bn}'s")
    print(f"   nodes {na:,} / {nb:,} ({nb/na:.2f}); CPU h {ca/3600:,.0f} / {cb/3600:,.0f}: {cb/ca:.2f} ({boots[50]:.2f}-{boots[950]:.2f}); "
          f"same node count in {len(same):,} ({sum(A['size'][k] for k in same):,} nodes), CPU on those {sum(B['cpu'][k] for k in same)/sum(A['cpu'][k] for k in same):.2f}")
    for d in (1, 2, 3, 4):
        ks = [k for k in top if A["nodes"][k]["depth"] == d]
        if ks:
            print(f"   root depth {d}: {len(ks):,} subtrees, nodes {sum(B['size'][k] for k in ks)/sum(A['size'][k] for k in ks):.2f}, CPU {sum(B['cpu'][k] for k in ks)/sum(A['cpu'][k] for k in ks):.2f}")
if __name__ == "__main__":
    for a in ("sko42-zone-261119-z25", "sko42-261120", "sko42-zone-261119-z10"):
        pair(a, "sko42-261120b")
