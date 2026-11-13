#!/usr/bin/env python3
"""Two arms' tree estimates from paired anchors, and their ratio with a paired interval:

  paired.py [estimate.py's options] A_DIR,B_DIR=SHORTFALL ...

Each argument is one anchor node probed under both arms (the same directory twice
for an anchor whose subtree the arms treat alike). It prints both arms' iso and fit
totals (as estimate_se.py does) and B/A. The interval resamples stage B's sample
and the anchor nodes with replacement, the same draws for both arms, so the
scatter that both arms share (which nodes are dear) mostly cancels in the ratio."""
import os, random, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estimate_se as T  # noqa: E402
E = T.E


def main():
    args = sys.argv[1:]
    table, top, top_cpu, front = None, None, None, None
    while args and args[0].startswith("--"):
        o = args.pop(0)
        if o == "--table": table = args.pop(0)
        elif o == "--top": top = int(args.pop(0))
        elif o == "--top-cpu": top_cpu = float(args.pop(0))
        elif o == "--front": front = int(args.pop(0))
        else: raise SystemExit(f"unknown option {o}")
    if None in (table, top, top_cpu, front):
        raise SystemExit("--table, --top, --top-cpu and --front are all needed")
    P = []
    for a in args:
        ds, _, s = a.rpartition("=")
        da, db = ds.split(",")
        P.append((float(s), E.anchor(da)[1], E.anchor(db)[1]))
    S = E.sample(table)
    rng = random.Random(1)
    B = 2000
    print(f"{len(P)} paired anchors, {len(S)} sampled frontier nodes, {B} draws")
    print("| model | A core-hours | B core-hours | B/A (95%, paired) | draws with B < A |")
    print("|---|---|---|---|---|")
    for mi, name in ((1, "iso"), (2, "fit")):
        tot = [T.total(S, T.models([(s, x[k]) for s, *x in P])[mi], 3, top_cpu, front) for k in (0, 1)]
        r = []
        for _ in range(B):
            Sb = [S[rng.randrange(len(S))] for _ in S]
            while True:
                Pb = [P[rng.randrange(len(P))] for _ in P]
                if len({x[0] for x in Pb}) >= 3:
                    break
            a, b = (T.total(Sb, T.models([(s, x[k]) for s, *x in Pb])[mi], 3, top_cpu, front) for k in (0, 1))
            r.append(b / a)
        r.sort()
        print(f"| {name} | {tot[0] / 3600:,.0f} | {tot[1] / 3600:,.0f} | {tot[1] / tot[0]:.2f}"
              f" ({r[int(0.025 * B)]:.2f}-{r[int(0.975 * B) - 1]:.2f}) | {sum(x < 1 for x in r) / B:.0%} |")


if __name__ == "__main__":
    main()
