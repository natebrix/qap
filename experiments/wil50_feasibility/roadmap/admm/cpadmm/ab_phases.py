#!/usr/bin/env python3
"""The node-10 A/B of the CP cut phase in the search (README.md, "In the search"):
per cut phase and per pooled child, from two arms' per-node traces
(pipeline_ab/run_pair.sh, runs/ARM/trace.txt). The totals and the kill line are
pipeline_ab/summarize.py's; this adds:

- every node with a phase in either arm: its first bound, the cutoff, and each
  arm's phase (iterations, separations or rounds, final bound, status), the CP
  phase's stop (fathomed, the stall rule, the budget, no cuts);
- the phase nodes' children (the nodes one fix below): count, fathomed, their
  first solves' iterations, in each arm;
- the pooled children (first solve with an inherited pool) on the nodes both
  arms bounded: first-solve iterations and CPU, and the outcome flips.

  python3 ab_phases.py RUNS [REF CP] [--budget 3400]
"""
import argparse
import os
import re


def num(pat, line, conv=float):
    m = re.search(pat, line)
    return conv(m.group(1)) if m else None


def read_trace(path):
    rows = {}
    for line in open(path):
        if " fixes " not in line:
            continue
        fx = frozenset(re.search(r"fixes ([\d: ]*?) \|", line).group(1).split())
        cut = num(r"\| cut (\d)", line, int) == 1
        rows[fx] = dict(
            fixes=fx, depth=num(r"depth (\d+)", line, int), lb=num(r"first lb (\S+)", line),
            status=num(r"\| first lb \S+ primal \S+ status (\S+)", line, str),
            iters=num(r"\| first lb \S+ primal \S+ status \S+ iters (\d+)", line, int),
            cpu=num(r"\| first lb \S+ primal \S+ status \S+ iters \d+ cpu (\S+)", line),
            pool=num(r" pool (\d) kept", line, int), pcp=" pcp 1 " in line,
            cut=cut, best=num(r"\| cut 1 best (\S+)", line) if cut else None,
            rounds=num(r"\| cut 1 best \S+ rounds (\d+)", line, int) if cut else 0,
            cut_iters=num(r"\| cut 1 best \S+ rounds \d+ iters (\d+)", line, int) if cut else 0,
            cut_cpu=num(r"\| cut 1 best \S+ rounds \d+ iters \d+ cpu (\S+)", line) if cut else 0.0,
            cut_status=num(r"\| cut 1 best \S+ rounds \d+ iters \d+ cpu \S+ status (\S+)", line, str) if cut else None,
            cp=num(r" cp (\d) active", line, int), active=num(r" cp \d active (\d+)", line, int),
            stalled=num(r" stalled (\d)", line, int),
            cutoff=num(r"cutoff (\S+)", line), fathomed=num(r"fathomed (\d) ", line, int),
            pool_out=num(r"pool_out (\d+)", line, int))
    return rows


def children(rows, fx):
    return [r for k, r in rows.items() if len(k) == len(fx) + 1 and fx < k]


def stop_reason(r, budget):
    if r["cut_status"] == "fathomed" or r["fathomed"] == 1:
        return "fathomed"
    if r["stalled"] == 1:
        return "stall rule"
    if r["cut_iters"] == 0:
        return "no cuts"
    if r["cut_iters"] >= budget:
        return "budget"
    return r["cut_status"] or "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("ref", nargs="?", default="bestcore")
    ap.add_argument("cp", nargs="?", default="cpcut")
    ap.add_argument("--budget", type=int, default=3400)
    a = ap.parse_args()
    X = read_trace(os.path.join(a.runs, a.ref, "trace.txt"))
    Y = read_trace(os.path.join(a.runs, a.cp, "trace.txt"))
    print(f"# {a.cp} against {a.ref}: cut phases and pooled children\n")
    ph = sorted({k for k, r in X.items() if r["cut"]} | {k for k, r in Y.items() if r["cut"]},
                key=lambda k: ((X.get(k) or Y.get(k))["depth"], sorted(k)))
    print(f"## The {len(ph)} nodes with a cut phase\n")
    print(f"| depth | node | first lb | cutoff | {a.ref}: iterations / rounds, best, status | "
          f"{a.cp}: iterations / separations, best, stop, active | best {a.cp} - {a.ref} |")
    print("|---|---|---|---|---|---|---|")
    for k in ph:
        x, y = X.get(k), Y.get(k)
        r = x or y
        dx = (f"{x['cut_iters']:,} / {x['rounds']}, {x['best']:.1f}, {x['cut_status']}" if x and x["cut"]
              else ("no phase" if x else "not bounded"))
        dy = (f"{y['cut_iters']:,} / {y['rounds']}, {y['best']:.1f}, {stop_reason(y, a.budget)}"
              + (f", {y['active']}" if y["active"] is not None else "") + ("" if y["cp"] != 0 else " (as rounds)")
              if y and y["cut"] else ("no phase" if y else "not bounded"))
        dd = f"{y['best'] - x['best']:+.1f}" if x and y and x["cut"] and y["cut"] else ""
        print(f"| {r['depth']} | {' '.join(sorted(k, key=lambda s: int(s.split(':')[0])))} | {r['lb']:.1f} "
              f"| {r['cutoff']:.1f} | {dx} | {dy} | {dd} |")
    both = [k for k in ph if k in X and k in Y and X[k]["cut"] and Y[k]["cut"]]
    if both:
        ix = sum(X[k]["cut_iters"] for k in both); iy = sum(Y[k]["cut_iters"] for k in both)
        cx = sum(X[k]["cut_cpu"] for k in both); cy = sum(Y[k]["cut_cpu"] for k in both)
        print(f"\n- the {len(both)} phases in both arms: iterations {ix:,} -> {iy:,}, CPU {cx:,.0f} -> {cy:,.0f} s "
              f"({1000 * cx / max(1, ix):.0f} -> {1000 * cy / max(1, iy):.0f} ms per iteration); fathomed "
              f"{sum(X[k]['fathomed'] == 1 for k in both)} -> {sum(Y[k]['fathomed'] == 1 for k in both)}; "
              f"mean best - first lb {sum(X[k]['best'] - X[k]['lb'] for k in both) / len(both):+.1f} -> "
              f"{sum(Y[k]['best'] - Y[k]['lb'] for k in both) / len(both):+.1f}; handed on a pool "
              f"{sum(X[k]['pool_out'] > 0 for k in both)} -> {sum(Y[k]['pool_out'] > 0 for k in both)}")
    cps = [r for r in Y.values() if r["cut"]]
    if cps:
        reasons = {}
        for r in cps:
            reasons.setdefault(stop_reason(r, a.budget), []).append(r)
        print(f"- {a.cp}'s {len(cps)} phases by stop: " + "; ".join(
            f"{s} {len(v)} ({sum(r['cut_iters'] for r in v):,} iterations)" for s, v in sorted(reasons.items())))
        bud = reasons.get("budget", [])
        if bud:
            print(f"  - at the budget ({a.budget:,}) unfathomed: " + "; ".join(
                f"[{' '.join(sorted(r['fixes']))}] {r['best']:.1f} of {r['cutoff']:.1f} (short {r['cutoff'] - r['best']:.1f})"
                for r in bud))
    print(f"\n## The phase nodes' children (one fix below)\n")
    print(f"| node | {a.ref}: children (fathomed), first-solve iterations | {a.cp}: children (fathomed), first-solve iterations |")
    print("|---|---|---|")
    for k in ph:
        cx, cy = children(X, k), children(Y, k)
        f = lambda c: (f"{len(c)} ({sum(r['fathomed'] == 1 for r in c)}), {sum(r['iters'] for r in c):,}" if c else "-")
        print(f"| {' '.join(sorted(k, key=lambda s: int(s.split(':')[0])))} | {f(cx)} | {f(cy)} |")
    print(f"\n## Pooled children (first solve with an inherited pool)\n")
    for nm, R in ((a.ref, X), (a.cp, Y)):
        p = [r for r in R.values() if r["pool"] == 1]
        print(f"- {nm}: {len(p)} pooled nodes ({sum(r['pcp'] for r in p)} in CP form), first solves "
              f"{sum(r['iters'] for r in p):,} iterations, {sum(r['cpu'] for r in p):,.0f} s, "
              f"{sum(r['status'] == 'fathomed' for r in p)} fathomed by it")
    sh = [k for k in X if k in Y and X[k]["pool"] == 1 and Y[k]["pool"] == 1]
    if sh:
        ix = sum(X[k]["iters"] for k in sh); iy = sum(Y[k]["iters"] for k in sh)
        cx = sum(X[k]["cpu"] for k in sh); cy = sum(Y[k]["cpu"] for k in sh)
        print(f"- the {len(sh)} nodes pooled in both: first-solve iterations {ix:,} -> {iy:,} "
              f"({ix / iy if iy else 0:.3f}x), CPU {cx:,.0f} -> {cy:,.0f} s ({cx / cy if cy else 0:.3f}x); "
              f"mean first lb {a.cp} - {a.ref} {sum(Y[k]['lb'] - X[k]['lb'] for k in sh) / len(sh):+.2f}")
        fl = [k for k in sh if X[k]["fathomed"] != Y[k]["fathomed"]]
        for k in fl:
            print(f"  - flip [{' '.join(sorted(k))}]: {a.ref} {'fathomed' if X[k]['fathomed'] else 'branched'} "
                  f"(lb {X[k]['lb']:.1f}), {a.cp} {'fathomed' if Y[k]['fathomed'] else 'branched'} (lb {Y[k]['lb']:.1f})")
    only = [k for k in Y if Y[k]["pool"] == 1 and (k not in X or X[k]["pool"] != 1)]
    if only:
        print(f"- {len(only)} nodes pooled in {a.cp} only (not bounded, or not pooled, in {a.ref}): first solves "
              f"{sum(Y[k]['iters'] for k in only):,} iterations")


if __name__ == "__main__":
    main()
