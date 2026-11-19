#!/usr/bin/env python3
"""Candidate subtrees for the tai60b / tai50b check: for every trace node N, its
subtree's size and CPU in the trace, and the nodes the pair test would drop
inside it (under the highest dominated nodes strictly below N), from
domcheck's per-node output. Prints the cheapest subtrees with drops.
  python3 subtree_pick.py INSTANCE TRACE.gz [min_saved]"""
import gzip, subprocess, sys, os
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
name, trace = sys.argv[1:3]; mins = int(sys.argv[3]) if len(sys.argv) > 3 else 5
out = subprocess.run(f"zcat {trace} | {HERE}/bin/domcheck /home/user/wt-dombuild/qaplib/{name}.dat | grep -v -e '^objective' -e '^symm'",
                     shell=True, capture_output=True, text=True).stdout.split("\n")
tr = [l.split() for l in gzip.open(trace, "rt") if l.startswith("node ")]
dom = {}
for t, l in zip(tr, out):
    f = l.split()
    dom[t[2]] = len(f) == 4 and f[3] == "1"
size = defaultdict(int); cpu = defaultdict(float); saved = defaultdict(int); savedcpu = defaultdict(float)
for t in tr:
    k = t[2]; parts = [] if k == "-" else k.split(",")
    pres = ["-"] + [",".join(parts[:i]) for i in range(1, len(parts) + 1)]
    # highest dominated prefix of k
    top = next((i for i, p in enumerate(pres) if dom.get(p)), None)
    for i, p in enumerate(pres):
        size[p] += 1; cpu[p] += float(t[11])
        if top is not None and top > i:
            saved[p] += 1; savedcpu[p] += float(t[11])
outc = {t[2]: t[7] for t in tr}
c = [(cpu[p], p) for p in size if saved[p] >= mins and not dom.get(p)]
c = [x for x in c if (0 if x[1]=="-" else len(x[1].split(","))) <= int(os.environ.get("MAXD","99")) and x[0] <= 3600*float(os.environ.get("MAXH","1e9"))]
c.sort(key=lambda x: -saved[x[1]] if os.environ.get("BYSAVED") else x[0])
for cp, p in c[:int(os.environ.get("TOP", "25"))]:
    print(f"depth {0 if p=='-' else len(p.split(','))} subtree {size[p]} cpu {cp/3600:.2f} h  saved {saved[p]} nodes {savedcpu[p]/3600:.2f} h  outcome {outc.get(p)}  {p}")
