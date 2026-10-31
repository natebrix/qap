#!/usr/bin/env python3
"""Summarize runs/ of run_all32.sh: per configuration, the nodes fathomed, the
median and total iterations and CPU to fathom (over the nodes every
configuration fathoms), give-ups, and for the misses the certified bound at
iteration 200 relative to d4."""
import glob, os, re, statistics, sys
D = os.path.dirname(os.path.abspath(__file__))
CUT = 48816 - 1.95


def parse(f):
    r = {"lb": {}}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m:
            r["lb"][int(m.group(1))] = float(m.group(2))
        if line.startswith("RESULT"):
            for kv in line.split()[1:]:
                k, v = kv.split("=")
                r[k] = v
    return r if "status" in r else None


runs = {}
for f in glob.glob(os.path.join(D, "runs", "*.log")):
    name, k = re.match(r"(.+)_(\d+)\.log", os.path.basename(f)).groups()
    r = parse(f)
    if r: runs.setdefault(name, {})[int(k)] = r
names = [n for n in ["d4", "d3", "d3L2", "d4L2", "d4W5", "d4A25"] if n in runs]
nodes = sorted(set.intersection(*[set(runs[n]) for n in names]))
fath = {n: {k for k in nodes if runs[n][k]["fathomed"] == "1"} for n in names}
common = set.intersection(*fath.values())
print("%d nodes done in every configuration; %d fathomed by all" % (len(nodes), len(common)))
print("%-6s %5s %8s %8s %9s %9s %8s  %s" % ("config", "fath", "med_it", "med_cpu", "sum_it",
                                             "sum_cpu", "miss_cpu", "lb@200 - d4 (misses)"))
base = runs[names[0]]
for n in names:
    R = runs[n]
    it = [int(R[k]["iters"]) for k in common]
    cpu = [float(R[k]["cpu"]) for k in common]
    misses = [k for k in nodes if k not in fath[n]]
    mcpu = sum(float(R[k]["cpu"]) for k in misses)
    d200 = [R[k]["lb"].get(200, float("nan")) - base[k]["lb"].get(200, float("nan"))
            for k in nodes if k not in fath[names[0]]]
    d200 = [x for x in d200 if x == x]
    print("%-6s %5d %8.0f %8.1f %9d %9.0f %8.0f  median %+.2f  min %+.2f  max %+.2f" % (
        n, len(fath[n]), statistics.median(it) if it else 0, statistics.median(cpu) if cpu else 0,
        sum(it), sum(cpu), mcpu, statistics.median(d200) if d200 else 0,
        min(d200) if d200 else 0, max(d200) if d200 else 0))
    lost = sorted(fath[names[0]] - fath[n]); won = sorted(fath[n] - fath[names[0]])
    if lost or won: print("       vs %s: lost %s, gained %s" % (names[0], lost, won))
