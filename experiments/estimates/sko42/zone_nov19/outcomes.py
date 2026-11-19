#!/usr/bin/env python3
"""The nodes both arms bounded: how each arm decided them (B branched, D the DNN bound fathomed, C the cut
phase fathomed, Q the QPB), and the CPU a node where the outcome is the same.

    python3 experiments/estimates/sko42/zone_nov19/outcomes.py [A_TRACE B_TRACE]
"""
import collections
import statistics as st
import sys

from common import ARMS, load


def main():
    pa, pb = sys.argv[1:3] if len(sys.argv) > 2 else (ARMS["z10"], ARMS["z25"])
    A, B = load(pa)["nodes"], load(pb)["nodes"]
    for name, X in (("A", A), ("B", B)):
        c = sum(1 for x in X.values() if x["cut"] is not None)
        print(f"{name}: {len(X):,} nodes, the cut phase run at {c:,} ({c / len(X):.1%}); "
              f"outcomes {dict(collections.Counter(x['out'] for x in X.values()))}")
    com = [k for k in A if k in B]
    print(f"nodes in both: {len(com):,}")
    print("outcome pairs (A, B): " + ", ".join(f"{a}/{b} {n:,}" for (a, b), n in
                                              collections.Counter((A[k]["out"], B[k]["out"]) for k in com).most_common(8)))
    for a, b in (("D", "D"), ("C", "D"), ("C", "C"), ("B", "B")):
        ks = [k for k in com if A[k]["out"] == a and B[k]["out"] == b]
        if ks:
            print(f"  {a}/{b}: {len(ks):,} nodes, B/A CPU {sum(B[k]['cpu'] for k in ks) / sum(A[k]['cpu'] for k in ks):.2f}; "
                  f"median CPU s A {st.median(A[k]['cpu'] for k in ks):.1f}, B {st.median(B[k]['cpu'] for k in ks):.1f}")


if __name__ == "__main__":
    main()
