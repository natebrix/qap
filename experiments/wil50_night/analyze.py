#!/usr/bin/env python3
"""wil50_night/analyze.py: every finished run under runs/ against the base arm of its subtree: nodes, DNN
iterations (first solves + cut phases), CPU, and the split of the base arm's work by kind of node."""
import os, re, sys, collections
D = os.path.dirname(os.path.abspath(__file__))


def load(d):
    rows = []
    p = os.path.join(d, "dnn.txt")
    if not os.path.exists(p): return rows
    for ln in open(p):
        m = re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.\-]+) primal [\d.\-]+ status (\w+) iters (\d+) cpu ([\d.]+)', ln)
        if not m: continue
        ph = re.search(r'\| cut 1 best ([\d.\-]+) rounds \d+ iters (\d+) cpu ([\d.]+)', ln)
        rows.append(dict(key=" ".join(sorted(m.group(2).split())), depth=int(m.group(1)), st=m.group(4), it=int(m.group(5)), cpu=float(m.group(6)),
                         pit=int(ph.group(2)) if ph else 0, pcpu=float(ph.group(3)) if ph else 0.0, ph=bool(ph), fath=' fathomed 1' in ln))
    return rows


def done(d):
    p = os.path.join(d, "log.txt")
    if not os.path.exists(p): return None
    t = open(p).read()
    m = re.search(r"Total Time:\s*([\d.]+)", t); n = re.search(r"Nodes:\s*(\d+)", t)
    return (int(n.group(1)), float(m.group(1)), "New objective" in t) if m and n and "exit 0" in t else None


runs = collections.defaultdict(dict)
for r in sorted(os.listdir(os.path.join(D, "runs"))):
    name, arm = r.split("_", 1)
    runs[name][arm] = os.path.join(D, "runs", r)
for name, arms in runs.items():
    if "base" not in arms or not done(arms["base"]):
        print(f"{name}: base not finished"); continue
    b = load(arms["base"]); bn, bcpu, _ = done(arms["base"])
    bit = sum(r["it"] + r["pit"] for r in b)
    kinds = collections.OrderedDict()
    for r in b:
        k = "first solve fathoms" if r["fath"] and not r["ph"] else "phase fathoms" if r["fath"] else "phase misses" if r["ph"] else "no phase, branched"
        e = kinds.setdefault(k, [0, 0, 0.0, 0, 0.0]); e[0] += 1; e[1] += r["it"]; e[2] += r["cpu"]; e[3] += r["pit"]; e[4] += r["pcpu"]
    tot = sum(e[2] + e[4] for e in kinds.values())
    print(f"## {name}: base {bn} nodes, {bit:,} DNN iterations, {bcpu:,.0f} CPU s\n")
    print("| kind of node (base) | nodes | first-solve iterations, CPU s | phase iterations, CPU s | share of DNN CPU |\n|---|---|---|---|---|")
    for k, e in kinds.items():
        print(f"| {k} | {e[0]} | {e[1]:,}, {e[2]:,.0f} | {e[3]:,}, {e[4]:,.0f} | {(e[2] + e[4]) / tot:.0%} |")
    print(f"\nFirst solves {sum(e[2] for e in kinds.values()) / tot:.0%} of the DNN CPU, phases {sum(e[4] for e in kinds.values()) / tot:.0%}; "
          f"DNN {tot / bcpu:.0%} of the run's CPU.\n")
    print("| arm | nodes | DNN iterations (first + phase) | CPU s | iterations / base | CPU / base | nodes both bounded: iterations / base |\n|---|---|---|---|---|---|---|")
    bk = {r["key"]: r for r in b}
    for arm, d in arms.items():
        dn = done(d)
        if not dn: print(f"| {arm} | running | | | | | |"); continue
        a = load(d); f = sum(r["it"] for r in a); p = sum(r["pit"] for r in a)
        common = [r for r in a if r["key"] in bk]
        ci = sum(r["it"] + r["pit"] for r in common); cb = sum(bk[r["key"]]["it"] + bk[r["key"]]["pit"] for r in common)
        print(f"| {arm} | {dn[0]} | {f + p:,} ({f:,} + {p:,}) | {dn[1]:,.0f} | {(f + p) / bit:.3f} | {dn[1] / bcpu:.3f} | "
              f"{len(common)} nodes: {ci / cb:.3f} |" + (" NEW OBJECTIVE" if dn[2] else ""))
    print()
