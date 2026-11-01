#!/usr/bin/env python3
"""Table of the depth-4 sample (nodes.txt; uncut_<id>.log from dnntest -v,
cuts_<id>.log from dnncuts from the give-up point): per node the parent, the
uncut certified bound at 400 and 1,500 iterations, the shortfall at
convergence, the best bound with cuts, its gain over the 400-iteration bound,
whether the node fathomed (uncut, or with cuts), rounds, iterations and CPU
of the cut run; then the counts by class.

  python3 summarize_depth4.py [dir] [cutoff]
"""
import os, re, sys

d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 48814.05
rows = []
for line in open(os.path.join(d, "nodes.txt")):
    k, parent, spec = line.split()
    k = int(k)
    uf = os.path.join(d, "uncut_%d.log" % k)
    cf = os.path.join(d, "cuts_%d.log" % k)
    if not (os.path.exists(uf) and os.path.exists(cf)):
        continue
    lb400 = lb1500 = None
    ustatus = "?"
    for l in open(uf):
        m = re.match(r"it\s+400 lb (\S+)", l)
        if m: lb400 = float(m.group(1))
        m = re.match(r"RESULT .*lb=(\S+) .*status=(\S+) iters=(\d+)", l)
        if m: lb1500 = float(m.group(1)); ustatus = m.group(2); uiters = int(m.group(3))
    rounds = []
    for l in open(cf):
        m = re.match(r"round\s+(\d+)\s+lb (\S+) .* status (\S+)\s+iters (\d+)\s+cpu (\S+)", l)
        if m: rounds.append((int(m.group(1)), float(m.group(2)), m.group(3), int(m.group(4)), float(m.group(5))))
    if not rounds: continue
    first = rounds[0][1]
    best = max(lb for _, lb, *_ in rounds)
    cfath = any(st == "fathomed" for _, _, st, _, _ in rounds)
    ufath = ustatus == "fathomed"
    rows.append(dict(k=k, parent=int(parent), lb400=lb400, lb1500=lb1500, short=cutoff - lb1500,
                     first=first, best=best, gain=best - first, ufath=ufath, uiters=uiters,
                     cfath=cfath, rounds=rounds[-1][0], its=sum(i for *_, i, _ in rounds),
                     cpu=sum(c for *_, c in rounds)))
print("node parent  uncut@400  uncut@1500  short   uncut fathom  | cut best   gain  fathom  rounds  iters  cpu_s")
for r in rows:
    print("%4d %6d  %9.1f  %10.1f  %6.1f  %-13s | %9.1f  %+6.1f  %-6s  %6d  %5d  %5.0f" % (
        r["k"], r["parent"], r["lb400"] if r["lb400"] else r["lb1500"], r["lb1500"], r["short"],
        ("yes (%d it)" % r["uiters"]) if r["ufath"] else "no", r["best"], r["gain"],
        "yes" if r["cfath"] else "no", r["rounds"], r["its"], r["cpu"]))
n = len(rows)
uf = [r for r in rows if r["ufath"]]
near = [r for r in rows if not r["ufath"]]
cf = [r for r in near if r["cfath"]]
print("\n%d children: %d fathom uncut (%d of them within 400 iterations); %d near misses, %.0f to %.0f short at convergence, of which %d fathom with cuts from the give-up point" % (
    n, len(uf), sum(1 for r in uf if r["uiters"] <= 400), len(near),
    min(r["short"] for r in near), max(r["short"] for r in near), len(cf)))
if cf:
    print("cut fathoms: short at convergence %s; mean %.0f iterations, %.0f CPU-s" % (
        ", ".join("%.0f" % r["short"] for r in cf), sum(r["its"] for r in cf) / len(cf), sum(r["cpu"] for r in cf) / len(cf)))
nf = [r for r in near if not r["cfath"]]
if nf:
    print("not fathomed: short at convergence %s; gains %s; %.0f CPU-s each" % (
        ", ".join("%.0f" % r["short"] for r in nf), ", ".join("%+.0f" % r["gain"] for r in nf), sum(r["cpu"] for r in nf) / len(nf)))
