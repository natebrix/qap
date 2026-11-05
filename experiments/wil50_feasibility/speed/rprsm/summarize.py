#!/usr/bin/env python3
"""Summarize runs/ of run_all32.sh (adapted from solvers/accel/summarize.py):
per configuration, the nodes fathomed, ms per iteration, median and total
iterations and CPU to fathom over the nodes every configuration fathoms, the
CPU on the misses, and for the baseline's misses the certified bound at
iteration 200 relative to the first configuration. Per-node lines with -v.
   usage: summarize.py [-v] base name [name ...]"""
import glob, os, re, statistics, sys
D = os.path.dirname(os.path.abspath(__file__))


def parse(f):
    r = {"trace": {}}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m:
            r["trace"][int(m.group(1))] = float(m.group(2))
        if line.startswith("RESULT"):
            for kv in line.split()[1:]:
                k, v = kv.split("=")
                r[k] = v
    return r if "status" in r else None


verbose = "-v" in sys.argv
names = [a for a in sys.argv[1:] if a != "-v"]
runs = {}
for f in glob.glob(os.path.join(D, "runs", "*.log")):
    name, k = re.match(r"(.+)_(\d+)\.log", os.path.basename(f)).groups()
    r = parse(f)
    if r: runs.setdefault(name, {})[int(k)] = r
names = [n for n in names if n in runs]
nodes = sorted(set.intersection(*[set(runs[n]) for n in names]))
fath = {n: {k for k in nodes if runs[n][k]["fathomed"] == "1"} for n in names}
common = set.intersection(*fath.values())
print("%d nodes done in every configuration (%s); %d fathomed by all" % (
    len(nodes), " ".join(map(str, nodes)), len(common)))
print("%-6s %5s %7s %8s %8s %8s %9s %8s  %s" % ("config", "fath", "ms/it", "med_it", "sum_it", "it/base",
                                               "sum_cpu", "miss_cpu", "lb@200 - base (base misses)"))
base = runs[names[0]]
bsum = sum(int(base[k]["iters"]) for k in common)
for n in names:
    R = runs[n]
    it = [int(R[k]["iters"]) for k in common]
    cpu = [float(R[k]["cpu"]) for k in common]
    misses = [k for k in nodes if k not in fath[n]]
    mcpu = sum(float(R[k]["cpu"]) for k in misses)
    d200 = [R[k]["trace"].get(200, float("nan")) - base[k]["trace"].get(200, float("nan"))
            for k in nodes if k not in fath[names[0]]]
    d200 = [x for x in d200 if x == x]
    msit = sum(float(R[k]["cpu"]) for k in nodes) / sum(int(R[k]["iters"]) for k in nodes) * 1e3
    print("%-6s %5d %7.0f %8.0f %8d %8.3f %9.0f %8.0f  median %+.2f  min %+.2f  max %+.2f" % (
        n, len(fath[n]), msit, statistics.median(it) if it else 0, sum(it), sum(it) / bsum if bsum else 0,
        sum(cpu), mcpu, statistics.median(d200) if d200 else 0,
        min(d200) if d200 else 0, max(d200) if d200 else 0))
    lost = sorted(fath[names[0]] - fath[n]); won = sorted(fath[n] - fath[names[0]])
    if lost or won: print("       vs %s: lost %s, gained %s" % (names[0], lost, won))
if verbose:
    print("node  " + "  ".join("%-26s" % n for n in names))
    for k in nodes:
        print("%-5d " % k + "  ".join("%-9s %4s %6.1f %7.1f" % (
            runs[n][k]["status"][:9], runs[n][k]["iters"], float(runs[n][k]["cpu"]),
            runs[n][k]["trace"].get(200, float("nan"))) for n in names))
