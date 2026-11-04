#!/usr/bin/env python3
"""Join child_gl.py's per-child bounds with the search trace (question 3).

For each child on the branching line: the parent's free bounds (outer-fixed,
exact per-child), and what the child cost in the search (first-solve
iterations and CPU, outcome). A child whose free bound exceeds the cutoff
would have needed no DNN solve at all.
Usage: join_children.py <childgl log> <trace> <parent fixes f:l,...> <row:F | col:L> [n=50]
"""
import re, sys

gl_log, trace, fixes, line = sys.argv[1:5]
n = int(sys.argv[5]) if len(sys.argv) > 5 else 50
fixed = [tuple(map(int, p.split(":"))) for p in fixes.split(",")]
free_f = [a for a in range(n) if a not in {f for f, _ in fixed}]
free_l = [b for b in range(n) if b not in {l for _, l in fixed}]
kind, v = line.split(":"); v = int(v)

bounds = {}
cutoff = None
for s in open(gl_log):
    m = re.match(r"parent .* cutoff (\S+)", s)
    if m:
        cutoff = float(m.group(1)); print(s.strip())
    m = re.match(r"child (\d+):(\d+) outer-fixed (\S+) exact-child (\S+)", s)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        bounds[(free_f[a], free_l[b])] = (float(m.group(3)), float(m.group(4)))

par = {f"{f}:{l}" for f, l in fixed}
rows = []
for s in open(trace):
    fx = re.search(r"fixes ([^|]*)\|", s).group(1).split()
    if len(fx) != len(fixed) + 1 or not par.issubset(fx):
        continue
    new = [x for x in fx if x not in par][0]
    f, l = map(int, new.split(":"))
    if (kind == "row" and f != v) or (kind == "col" and l != v):
        continue
    it = int(re.search(r"iters (\d+)", s).group(1)); cpu = float(re.search(r"cpu ([\d.]+)", s).group(1))
    st = re.search(r"status (\w+)", s).group(1); fat = re.search(r"fathomed (\d)", s).group(1)
    ob, eb = bounds.get((f, l), (float("nan"), float("nan")))
    rows.append((it, new, st, fat, cpu, ob, eb))
rows.sort()
saved_cpu = saved_it = 0.0
nf = 0
for it, new, st, fat, cpu, ob, eb in rows:
    hit = eb > cutoff
    if hit:
        nf += 1; saved_cpu += cpu; saved_it += it
    print(f"{new:6s} search: {st:9s} it {it:4d} cpu {cpu:7.1f} | free bounds outer-fixed {cutoff - ob:8.2f} short, "
          f"exact-child {cutoff - eb:8.2f} short {'<- FREE FATHOM' if hit else ''}")
tot = sum(r[4] for r in rows)
print(f"children {len(rows)}, free fathoms {nf}, their search CPU {saved_cpu:.1f} s of {tot:.1f} s "
      f"({saved_it:.0f} iterations)")
