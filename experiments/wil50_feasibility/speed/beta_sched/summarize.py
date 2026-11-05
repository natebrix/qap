#!/usr/bin/env python3
"""Summarize runs/ of run_all32.sh against the first configuration named:
  summarize.py base two05 ... [--nodes 1,2,9]
Per configuration: the nodes fathomed (and lost or gained against the first),
iterations to fathom on the nodes both fathom (median, total, per-node
difference), the certified lb at iterations 100 and 200 on the first
configuration's misses relative to it, ms per iteration and total CPU over
the nodes every configuration finished. --table also prints every node."""
import glob, os, re, statistics, sys

D = os.path.dirname(os.path.abspath(__file__))


def parse(f):
    r = {"trace": {}, "beta": []}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m:
            r["trace"][int(m.group(1))] = float(m.group(2))
        m = re.match(r"bsched it (\d+) beta \S+ -> (\S+)", line)
        if m:
            r["beta"].append((int(m.group(1)), float(m.group(2))))
        if line.startswith("RESULT"):
            for kv in line.split()[1:]:
                k, v = kv.split("=")
                r[k] = v
    return r if "status" in r else None


args = sys.argv[1:]
table = "--table" in args
if table: args.remove("--table")
only = None
if "--nodes" in args:
    k = args.index("--nodes"); only = {int(x) for x in args[k + 1].split(",")}; del args[k:k + 2]
runs = {}
for f in glob.glob(os.path.join(D, "runs", "*.log")):
    name, k = re.match(r"(.+)_(\d+)\.log", os.path.basename(f)).groups()
    r = parse(f)
    if r: runs.setdefault(name, {})[int(k)] = r
names = [n for n in args if n in runs]
nodes = sorted(set.intersection(*[set(runs[n]) for n in names]))
if only: nodes = [k for k in nodes if k in only]
fath = {n: {k for k in nodes if runs[n][k]["fathomed"] == "1"} for n in names}
b = names[0]
B = runs[b]
print("%d nodes done in every configuration: %s" % (len(nodes), ",".join(map(str, nodes))))
print("%-8s %4s %6s %7s %7s %8s %9s %9s  %-28s %-28s" % (
    "config", "fath", "ms/it", "med_it", "sum_it", "d_it_med", "early>=10", "late>=10",
    "lb@100 - %s (misses)" % b, "lb@200 - %s (misses)" % b))
for n in names:
    R = runs[n]
    common = sorted(fath[n] & fath[b])
    it = [int(R[k]["iters"]) for k in common]
    dit = [int(R[k]["iters"]) - int(B[k]["iters"]) for k in common]
    misses = [k for k in nodes if k not in fath[b]]

    def dl(at):
        d = [R[k]["trace"].get(at, float("nan")) - B[k]["trace"].get(at, float("nan")) for k in misses]
        d = [x for x in d if x == x]
        return "med %+7.2f [%+7.2f, %+7.2f]" % (statistics.median(d), min(d), max(d)) if d else "-"
    msit = sum(float(R[k]["cpu"]) for k in nodes) / sum(int(R[k]["iters"]) for k in nodes) * 1e3
    print("%-8s %4d %6.0f %7.0f %7d %+8.0f %9d %9d  %-28s %-28s" % (
        n, len(fath[n]), msit, statistics.median(it) if it else 0, sum(it),
        statistics.median(dit) if dit else 0, sum(1 for x in dit if x <= -10),
        sum(1 for x in dit if x >= 10), dl(100), dl(200)))
    lost = sorted(fath[b] - fath[n]); won = sorted(fath[n] - fath[b])
    if lost or won: print("         vs %s: lost %s, gained %s" % (b, lost, won))
    print("         total CPU %.0f s (fathoms %.0f, misses %.0f)" % (
        sum(float(R[k]["cpu"]) for k in nodes), sum(float(R[k]["cpu"]) for k in nodes if k in fath[n]),
        sum(float(R[k]["cpu"]) for k in nodes if k not in fath[n])))
if table:
    print("\nnode " + " ".join("%-30s" % n for n in names))
    for k in nodes:
        row = []
        for n in names:
            r = runs[n][k]
            row.append("%-30s" % ("%s %4s lb100 %.1f lb200 %.1f" % (
                "F" if r["fathomed"] == "1" else "m", r["iters"], r["trace"].get(100, float("nan")),
                r["trace"].get(200, float("nan")))))
        print("%4d " % k + " ".join(row))
