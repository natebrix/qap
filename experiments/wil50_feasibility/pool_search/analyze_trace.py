#!/usr/bin/env python3
"""Per-node summary of QAP_DNN_TRACE lines (bnb.C): nodes with and without
an inherited pool, fathoms by the first solve and by the cut phase, DNN
iterations and CPU, the pools handed on.
  python3 analyze_trace.py trace.txt"""
import collections, re, sys

rows = []
for line in open(sys.argv[1] if len(sys.argv) > 1 else "trace.txt"):
    g = lambda pat, conv=float: (lambda m: conv(m.group(1)) if m else None)(re.search(pat, line))
    rows.append(dict(
        depth=g(r"depth (\d+)", int),
        lb=g(r"first lb (\S+)"), primal=g(r"primal (\S+)"), status=g(r"status (\S+)", str),
        iters=g(r"iters (\d+)", int), cpu=g(r"cpu (\S+)"), warm=g(r"warm (\d)", int),
        pool=g(r"pool (\d)", int), kept=g(r"kept (\d+)", int), cut=g(r"\| cut (\d)", int),
        cut_best=g(r"\| cut 1 best (\S+)"), cut_iters=g(r"rounds \d+ iters (\d+)", int),
        cut_cpu=g(r"rounds \d+ iters \d+ cpu (\S+)"), fathomed=g(r"fathomed (\d) ", int),
        pool_out=g(r"pool_out (\d+)", int), cutoff=g(r"cutoff (\S+)")))

print(f"{len(rows)} DNN nodes, {sum(r['fathomed'] for r in rows)} fathomed")
by = collections.defaultdict(list)
for r in rows:
    by[(r["pool"], r["depth"])].append(r)
print("\n| inherited pool | depth | nodes | fathomed by the first solve | cut phases | fathomed by the phase | branched | first-solve iterations (mean) | DNN CPU s (first + phase) |")
print("|---|---|---|---|---|---|---|---|---|")
for (p, d) in sorted(by):
    rs = by[(p, d)]
    f1 = sum(r["status"] == "fathomed" for r in rs)
    cp = sum(r["cut"] == 1 for r in rs)
    fc = sum(r["cut"] == 1 and r["fathomed"] == 1 for r in rs)
    br = sum(r["fathomed"] == 0 for r in rs)
    it = sum(r["iters"] for r in rs) / len(rs)
    cpu = sum(r["cpu"] + (r["cut_cpu"] or 0) for r in rs)
    print(f"| {'yes' if p else 'no'} | {d} | {len(rs)} | {f1} | {cp} | {fc} | {br} | {it:.0f} | {cpu:,.0f} |")
tot_first = sum(r["cpu"] for r in rows)
tot_cut = sum(r["cut_cpu"] or 0 for r in rows)
print(f"\nDNN CPU: first solves {tot_first:,.0f} s, cut phases {tot_cut:,.0f} s, total {tot_first + tot_cut:,.0f} s")
print(f"first-solve iterations {sum(r['iters'] for r in rows):,}, cut-phase iterations {sum(r['cut_iters'] or 0 for r in rows):,}")
pooled = [r for r in rows if r["pool"] == 1]
if pooled:
    print(f"pooled nodes: {len(pooled)}, kept cuts median {sorted(r['kept'] for r in pooled)[len(pooled)//2]}, "
          f"fathomed by the first solve {sum(r['status'] == 'fathomed' for r in pooled)}, mean first-solve iterations "
          f"{sum(r['iters'] for r in pooled)/len(pooled):.0f}")
print(f"pools handed on: {sum(r['pool_out'] > 0 for r in rows)} (from a cut phase {sum(r['pool_out'] > 0 and r['cut'] == 1 for r in rows)})")
print("\ncut phases:")
for r in rows:
    if r["cut"] == 1:
        print(f"  depth {r['depth']} pool {r['pool']} first {r['lb']:.2f} best {r['cut_best']:.2f} (+{r['cut_best'] - r['lb']:.2f}) "
              f"iters {r['cut_iters']} cpu {r['cut_cpu']:.0f} fathomed {r['fathomed']} pool_out {r['pool_out']}")
