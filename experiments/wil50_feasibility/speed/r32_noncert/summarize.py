#!/usr/bin/env python3
"""Summarize runs/ of run_all32.sh (base vs nc = DNN_R32_NONCERT=1): fathom
set, iterations to fathom, ms per iteration, the misses' certified bound at
iterations 200 and 300, the certificate margin on certificate iterations, and
the form-R split of the -v -v profile. Usage: summarize.py [base nc]."""
import glob, os, re, statistics, sys
D = os.path.dirname(os.path.abspath(__file__))


def parse(f):
    r = {"trace": {}, "margins": []}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb (\S+) \(cert (\S+), margin (\S+?)[,)]", line)
        if m:
            r["trace"][int(m.group(1))] = float(m.group(2))
            r["margins"].append(float(m.group(4)))
        m = re.match(r"form R \(ms per iteration\): F (\S+)\s+syrk (\S+)\s+symmetrize (\S+)", line)
        if m:
            r["formR"] = [float(x) for x in m.groups()]
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
names = [n for n in (sys.argv[1:] or ["base", "nc"]) if n in runs]
nodes = sorted(set.intersection(*[set(runs[n]) for n in names]))
fath = {n: {k for k in nodes if runs[n][k]["fathomed"] == "1"} for n in names}
common = set.intersection(*fath.values())
base = runs[names[0]]
print("%d nodes done in every configuration; %d fathomed by all" % (len(nodes), len(common)))
print("%-5s %5s %7s %8s %8s %9s %9s %9s %10s %10s" % (
    "cfg", "fath", "ms/it", "med_it", "sum_it", "cpu_fath", "cpu_miss", "cpu_all", "margin_med", "margin_max"))
for n in names:
    R = runs[n]
    it = [int(R[k]["iters"]) for k in common]
    cpu = [float(R[k]["cpu"]) for k in common]
    misses = [k for k in nodes if k not in fath[n]]
    mcpu = sum(float(R[k]["cpu"]) for k in misses)
    acpu = sum(float(R[k]["cpu"]) for k in nodes)
    msit = acpu / sum(int(R[k]["iters"]) for k in nodes) * 1e3
    mg = [x for k in nodes for x in R[k]["margins"]]
    print("%-5s %5d %7.1f %8.0f %8d %9.1f %9.1f %9.1f %10.2e %10.2e" % (
        n, len(fath[n]), msit, statistics.median(it) if it else 0, sum(it), sum(cpu), mcpu, acpu,
        statistics.median(mg), max(mg)))
    lost = sorted(fath[names[0]] - fath[n]); won = sorted(fath[n] - fath[names[0]])
    if lost or won: print("      vs %s: lost %s, gained %s" % (names[0], lost, won))
for n in names[1:]:
    R = runs[n]
    ratio = [(float(R[k]["cpu"]) / int(R[k]["iters"])) / (float(base[k]["cpu"]) / int(base[k]["iters"]))
             for k in nodes]
    print("%s/%s ms per iteration, per node: median %.4f  min %.4f  max %.4f" % (
        n, names[0], statistics.median(ratio), min(ratio), max(ratio)))
    same = sum(1 for k in common if R[k]["iters"] == base[k]["iters"])
    d_it = [int(R[k]["iters"]) - int(base[k]["iters"]) for k in sorted(common)]
    print("iterations to fathom, %s - %s: same on %d of %d; differences %s" % (
        n, names[0], same, len(common), [d for d in d_it if d]))
    pairs = [(a, b) for k in nodes for a, b in zip(R[k]["margins"], base[k]["margins"])]
    rel = [abs(a - b) / b for a, b in pairs]
    print("certificate margins, %s vs %s, paired by iteration: %d pairs, max |ratio - 1| %.3f, "
          "%d differ by more than 1%%; %s max %.2e" % (
              n, names[0], len(pairs), max(rel), sum(1 for x in rel if x > 0.01), n, max(a for a, _ in pairs)))
    for t in (200, 300):
        d = [R[k]["trace"][t] - base[k]["trace"][t] for k in nodes
             if k not in fath[names[0]] and t in R[k]["trace"] and t in base[k]["trace"]]
        if d:
            print("misses' lb at %d, %s - %s: median %+.3f  min %+.3f  max %+.3f  (%d nodes)" % (
                t, n, names[0], statistics.median(d), min(d), max(d), len(d)))
print("form R, ms per iteration (F, syrk, symmetrize), mean over nodes:")
for n in names:
    fr = [runs[n][k]["formR"] for k in nodes if "formR" in runs[n][k]]
    if fr:
        print("  %-5s F %.1f  syrk %.1f  symmetrize %.1f  total %.1f" % (
            n, *[statistics.mean(x[i] for x in fr) for i in range(3)],
            statistics.mean(sum(x) for x in fr)))
print("per node (id: outcome iters ms/it, per configuration):")
for k in nodes:
    print("  %2d " % k + "  ".join("%s %s %4s %6.1f" % (
        n, "F" if runs[n][k]["fathomed"] == "1" else "-", runs[n][k]["iters"],
        float(runs[n][k]["cpu"]) / int(runs[n][k]["iters"]) * 1e3) for n in names))
