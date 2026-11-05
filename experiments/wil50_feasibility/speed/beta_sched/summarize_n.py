#!/usr/bin/env python3
"""Per-node comparison of dnntest logs holding several bounds (-N node files,
-K children): summarize_n.py base.log other.log ...
Each RESULT line closes one bound; prints status, iterations, the certified
lb at iterations 100 and 200 and the final lb, then totals over the bounds
that fathom in every log."""
import re, sys


def parse(f):
    out, tr = [], {}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m: tr[int(m.group(1))] = float(m.group(2))
        if line.startswith("RESULT"):
            r = dict(kv.split("=") for kv in line.split()[1:])
            r["trace"] = tr; out.append(r); tr = {}
    return out


logs = [parse(f) for f in sys.argv[1:]]
n = min(len(x) for x in logs)
names = [f.split("/")[-1].replace(".log", "") for f in sys.argv[1:]]
print("bound " + " | ".join("%-38s" % x for x in names))
for k in range(n):
    row = []
    for L in logs:
        r = L[k]
        row.append("%-38s" % ("%s %4s lb100 %8.1f lb200 %8.1f end %8.1f" % (
            "F" if r["status"] == "fathomed" else "m", r["iters"], r["trace"].get(100, float("nan")),
            r["trace"].get(200, float("nan")), float(r["lb"]))))
    print("%5d " % k + " | ".join(row))
common = [k for k in range(n) if all(L[k]["status"] == "fathomed" for L in logs)]
for name, L in zip(names, logs):
    it = sum(int(L[k]["iters"]) for k in common)
    cpu = sum(float(L[k]["cpu"]) for k in common)
    fath = sum(1 for k in range(n) if L[k]["status"] == "fathomed")
    print("%-10s fathomed %d of %d; on the %d fathomed by all: %d iterations, %.1f s; all bounds: %d iterations, %.1f s" % (
        name, fath, n, len(common), it, cpu, sum(int(L[k]["iters"]) for k in range(n)),
        sum(float(L[k]["cpu"]) for k in range(n))))
