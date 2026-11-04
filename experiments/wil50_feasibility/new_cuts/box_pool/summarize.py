#!/usr/bin/env python3
"""Table and gate of run_box_pool.sh: per node and arm the round-0 bound, the
best Polyak round, the quiet round, the retained certificate (FINAL best_lb),
fathom, iterations and process CPU; the gate is the median over nodes 10, 13
and d2 of (box - tri) retained certificates, KILL below +15.
  python3 summarize.py [runs/]"""
import os, re, statistics, sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
CUTOFF = 48814.05


def parse(path):
    if not os.path.exists(path):
        return None
    t = open(path).read()
    r = {"host": (re.search(r"^host\s+(.*)$", t, re.M) or [None, "?"])[1]}
    rounds = [(int(m.group(1)), float(m.group(2)), int(m.group(3)))
              for m in re.finditer(r"^round\s+(\d+)\s+lb (\S+) .*?iters (\d+)", t, re.M)]
    r["first"] = rounds[0][1] if rounds else None
    r["best_round"] = max((lb for k, lb, _ in rounds if k > 0), default=None)
    q = re.search(r"^quiet round\s+lb (\S+) .*?iters (\d+)", t, re.M)
    r["quiet"] = float(q.group(1)) if q else None
    f = re.search(r"^FINAL best_lb (\S+)", t, re.M)
    r["final"] = float(f.group(1)) if f else None
    r["iters"] = sum(i for _, _, i in rounds) + (int(q.group(2)) if q else 0)
    c = re.search(r"PROCESS_CPU user (\S+) sys (\S+)", t)
    r["cpu"] = float(c.group(1)) + float(c.group(2)) if c else None
    b = [int(m.group(1)) for m in re.finditer(r"addedB (\d+)", t)]
    r["boxes"] = sum(b)
    return r


rows, gains = [], {}
print("| node | arm | round 0 | best Polyak round | quiet round | retained | fathomed | iterations | process CPU s | box cuts added | host |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for k in ["10", "13", "d2", "14"]:
    res = {}
    for arm in ["tri", "box"]:
        r = parse(os.path.join(D, f"{arm}_{k}.log"))
        res[arm] = r
        if r is None:
            print(f"| {k} | {arm} | (not done) |||||||||")
            continue
        fmt = lambda v: "-" if v is None else f"{v:,.2f}"
        print(f"| {k} | {arm} | {fmt(r['first'])} | {fmt(r['best_round'])} | {fmt(r['quiet'])} | **{fmt(r['final'])}** | "
              f"{'yes' if r['final'] and r['final'] > CUTOFF else 'no'} | {r['iters']:,} | {fmt(r['cpu'])} | {r['boxes']} | {r['host']} |")
    if res["tri"] and res["box"] and res["tri"]["final"] and res["box"]["final"]:
        gains[k] = res["box"]["final"] - res["tri"]["final"]
print()
for k, g in gains.items():
    print(f"node {k}: box - tri = {g:+.2f}")
gate = [gains[k] for k in ["10", "13", "d2"] if k in gains]
if len(gate) == 3:
    med = statistics.median(gate)
    print(f"gate: median over 10, 13, d2 = {med:+.2f} -> {'PASS' if med >= 15 else 'KILL'} (threshold +15)")
else:
    print("gate: incomplete")
