#!/usr/bin/env python3
# Table of a directory of dnntest -v logs (node<k>.log, specs in nodes.txt): certified bound at
# 400 iterations (the search's give-up point) and at the end, short of the cutoff, status, CPU.
import re, sys, os, statistics
d = sys.argv[1]; cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 23386 - 1.95
specs = [l.split()[-1] for l in open(os.path.join(d, "nodes.txt")) if l.strip()]
rows = []
for k, spec in enumerate(specs, 1):
    f = os.path.join(d, f"node{k}.log")
    if not os.path.exists(f): continue
    txt = open(f).read()
    its = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"^it +(\d+) lb ([\d.]+)", txt, re.M)}
    lb400 = its.get(400, max((v for i, v in its.items() if i <= 400), default=float("nan")))
    m = re.search(r"RESULT node=\d+ m=(\d+) lb=([\d.]+) primal=([\d.]+) status=(\w+) iters=(\d+) cpu=([\d.]+)", txt)
    if not m:
        last = max(its) if its else 0
        rows.append((k, spec, lb400, its.get(last, float("nan")), "running", last, float("nan"), float("nan"))); continue
    rows.append((k, spec, lb400, float(m.group(2)), m.group(4), int(m.group(5)), float(m.group(6)), float(m.group(3))))
print(f"| node | spec | lb at 400 | short at 400 | final lb | short | primal | status | iters | CPU s |")
print("|---|---|---|---|---|---|---|---|---|---|")
for k, spec, a, b, st, it, cpu, pr in rows:
    print(f"| {k} | `{spec}` | {a:,.1f} | {cutoff-a:,.1f} | {b:,.1f} | {cutoff-b:,.1f} | {pr:,.1f} | {st} | {it} | {cpu:,.0f} |")
done = [r for r in rows if r[4] != "running"]
if done:
    sh = sorted(cutoff - r[3] for r in done)
    print(f"\n{len(done)} done: fathomed {sum(1 for r in done if r[4]=='fathomed')}, short values {', '.join(f'{x:.0f}' for x in sh)}, median {statistics.median(sh):.0f}")
