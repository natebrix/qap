"""Pick a tuning set from a complete node trace (docs/tho40b_results_plan.md,
(a) 4): subtrees whose size is known exactly, each about the length of one
tuning run, spread over the range of gaps, none inside another.

A subtree's CPU is the sum of its nodes' CPU in the trace. Its wall on one
instance is CPU / (slots x busy), busy being the share of slots a subtree run
keeps busy (tho40ab-32: 74%). Candidates are branched nodes, not the trace's
root, whose subtree is complete and whose wall falls within --tolerance of
--minutes. They are sorted by gap and cut into --n bins of equal count; each bin
gives the candidate closest to the target wall that is not inside, or above, one
already picked.

usage: python3 tuning_set.py TRACE[.gz] [--minutes 30] [--slots 64] [--busy 0.74] [--n 8] [--tolerance 0.5]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ntrace import load  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace")
    ap.add_argument("--minutes", type=float, default=30)
    ap.add_argument("--slots", type=int, default=64)
    ap.add_argument("--busy", type=float, default=0.74)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--tolerance", type=float, default=0.5, help="accept walls within this fraction of the target")
    a = ap.parse_args()
    t = load(a.trace)
    roots = set(t.roots)
    rate = a.slots * a.busy                      # CPU s per wall s
    target = a.minutes * 60 * rate
    lo, hi = target * (1 - a.tolerance), target * (1 + a.tolerance)
    cand = sorted((i for i in range(len(t)) if t.outcome[i] == "B" and i not in roots and t.complete[i]
                   and lo <= t.sub_cpu[i] <= hi), key=t.gap_pct)
    print(f"{a.trace}: {len(t)} nodes; target {a.minutes:g} min on {a.slots} slots at {a.busy:.0%} busy "
          f"= {target:,.0f} CPU s (accepting {lo:,.0f} to {hi:,.0f}); {len(cand)} candidates")
    if not cand:
        return
    picked = []

    def clashes(i):
        anc = set(t.ancestors(i))
        return any(p in anc or i in set(t.ancestors(p)) for p in picked)
    n = min(a.n, len(cand))
    for k in range(n):
        bin_ = cand[k * len(cand) // n:(k + 1) * len(cand) // n]
        for i in sorted(bin_, key=lambda i: abs(t.sub_cpu[i] - target)):
            if not clashes(i):
                picked.append(i)
                break
    print(f"{'depth':>5} {'gap %':>7} {'nodes':>8} {'CPU s':>9} {'wall min':>9} {'max depth':>9}  -F")
    for i in sorted(picked, key=t.gap_pct):
        print(f"{t.depth[i]:>5} {t.gap_pct(i):>7.3f} {t.size[i]:>8,} {t.sub_cpu[i]:>9,.0f} "
              f"{t.sub_cpu[i] / rate / 60:>9.1f} {t.max_depth[i]:>9}  {t.fixes_flag(i)}")


if __name__ == "__main__":
    main()
