#!/usr/bin/env python3
"""Summarize warm/ of run_warmset.sh against the first configuration named:
  summarize_warm.py prod cross w05 ... [--table] [--families f9,f6]
Per configuration, over the families every configuration finished: the parents'
iterations and final lb; the children fathomed (lost/gained against the first);
children's iterations to fathom on the common fathoms (total, per-child
difference, how many earlier/later by >= 10); the certified lb at iterations 50,
100 and 200 on the first configuration's missed children, relative to it; CPU
(children, parents) and the extra certificates (DNN_ADMM_CROSS)."""
import glob, os, re, statistics, sys

D = os.path.dirname(os.path.abspath(__file__))


def parse(f):
    """-> list of solves: {node, iters, fathomed, lb, cpu, trace{it: cert}, extra}"""
    solves, cur = [], None
    for line in open(f, errors="replace"):
        if line.startswith("node ") or line.startswith("CHILD-PASS"):
            cur = {"trace": {}, "extra": 0}
        m = re.match(r"it\s+(\d+) lb (\S+) \(cert ([0-9.eE+-]+)", line)
        if m and cur is not None:
            if "extra" in line: cur["extra"] += 1
            else: cur["trace"][int(m.group(1))] = float(m.group(2))   # running best
        m = re.match(r"admm cross: (\d+) extra", line)
        if m and cur is not None: cur["extra"] = int(m.group(1))
        if line.startswith("RESULT") and cur is not None:
            kv = dict(x.split("=", 1) for x in line.split()[1:] if "=" in x)
            cur.update(node=int(kv["node"]), iters=int(kv["iters"]), lb=float(kv["lb"]), cpu=float(kv["cpu"]),
                       fathomed=kv.get("fathomed", "1" if kv["status"] == "fathomed" else "0") == "1")
            solves.append(cur); cur = None
    return solves if len(solves) == 5 else None


args = sys.argv[1:]
table = "--table" in args
if table: args.remove("--table")
fams = None
if "--families" in args:
    k = args.index("--families"); fams = set(args[k + 1].split(",")); del args[k:k + 2]
runs = {}
for f in glob.glob(os.path.join(D, "warm", "*.log")):
    name, fam = re.match(r"(.+)_([a-z0-9]+)\.log", os.path.basename(f)).groups()
    s = parse(f)
    if s: runs.setdefault(name, {})[fam] = s
names = [n for n in args if n in runs]
common = sorted(set.intersection(*[set(runs[n]) for n in names]))
if fams: common = [f for f in common if f in fams]
print("families done in every configuration: %s" % ",".join(common))
b = names[0]
keys = [(fam, k) for fam in common for k in range(1, 5)]
fath = {n: {(fam, k) for fam, k in keys if runs[n][fam][k]["fathomed"]} for n in names}
print("%-8s %-22s %4s %7s %8s %6s %6s  %-26s %-26s %-26s %9s %7s %5s" % (
    "config", "parents: it (lb-%s)" % b, "fath", "sum_it", "d_it_med", "e>=10", "l>=10",
    "lb@50 - %s (misses)" % b, "lb@100", "lb@200", "CPU kids", "CPU par", "extra"))
for n in names:
    R = runs[n]; B = runs[b]
    par_it = sum(R[fam][0]["iters"] for fam in common)
    par_lb = [R[fam][0]["lb"] - B[fam][0]["lb"] for fam in common]
    cf = sorted(fath[n] & fath[b])
    it = [R[fam][k]["iters"] for fam, k in cf]
    dit = [R[fam][k]["iters"] - B[fam][k]["iters"] for fam, k in cf]
    misses = [(fam, k) for fam, k in keys if (fam, k) not in fath[b]]

    def dl(at):
        d = [R[fam][k]["trace"].get(at, float("nan")) - B[fam][k]["trace"].get(at, float("nan")) for fam, k in misses]
        d = [x for x in d if x == x]
        return "med %+7.2f [%+6.1f,%+6.1f]" % (statistics.median(d), min(d), max(d)) if d else "-"
    print("%-8s %5d (%+6.2f..%+6.2f) %4d %7d %+8.1f %6d %6d  %-26s %-26s %-26s %9.0f %7.0f %5d" % (
        n, par_it, min(par_lb), max(par_lb), len(fath[n]), sum(it), statistics.median(dit) if dit else 0,
        sum(1 for x in dit if x <= -10), sum(1 for x in dit if x >= 10), dl(50), dl(100), dl(200),
        sum(R[fam][k]["cpu"] for fam, k in keys), sum(R[fam][0]["cpu"] for fam in common),
        sum(R[fam][k]["extra"] for fam, k in keys)))
    lost = sorted(fath[b] - fath[n]); won = sorted(fath[n] - fath[b])
    if lost or won: print("         vs %s: lost %s, gained %s" % (b, lost, won))
if table:
    print("\nchild      " + " ".join("%-34s" % n for n in names))
    for fam, k in keys:
        row = []
        for n in names:
            r = runs[n][fam][k]
            row.append("%-34s" % ("%s %4d it lb100 %.1f lb200 %.1f x%d" % (
                "F" if r["fathomed"] else "m", r["iters"], r["trace"].get(100, float("nan")),
                r["trace"].get(200, float("nan")), r["extra"])))
        print("%-4s %d     " % (fam, k) + " ".join(row))
