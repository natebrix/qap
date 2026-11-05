#!/usr/bin/env python3
"""Summary of the pipeline A/B arms (run_pair.sh): per arm the tree, the DNN
iterations and CPU from the node traces (QAP_DNN_TRACE lines, bnb.C) and the
coordinator's totals; then every arm against base (and warmy against zone) on
the whole tree and on the nodes both trees share.

  python3 summarize.py RUNS_DIR        (default: runs/ next to this file)

Iterations do not depend on the machine's speed; CPU does, so CPU ratios are
meaningful only between arms run on the same machine, best at the same time
(hosts.txt in each arm's directory)."""
import os
import re
import sys

ARMS = ["base", "zone", "warmy", "certx", "ship", "full", "zone2", "base2", "ship2", "fast"]
# kill lines (README.md): the ratio of DNN iterations an arm must reach. For
# warmy and certx the end-to-end ratio is below the run-to-run noise of one
# configuration (README, base vs zone), so they are judged on the matched
# nodes (shared, same parent-Y warm start and pool status) and, for warmy, on
# the phases both arms ran
GATES = {("zone", "base"): 1.15, ("warmy", "zone"): 1.05, ("certx", "zone"): 1.03, ("ship", "zone"): 1.1,
         ("full", "ship"): 1.1, ("full", "warmy"): 1.1}
# the pairs compared, (arm, reference); the first reference that finished is used
PAIRS = [("zone", ["base"]), ("warmy", ["zone"]), ("certx", ["zone", "base"]), ("ship", ["zone"]), ("full", ["ship"]), ("full", ["warmy"]),
         ("full", ["zone"]), ("full", ["base"]), ("zone2", ["zone"]), ("ship2", ["ship"]),
         ("base2", ["base"]), ("fast", ["full"])]


def num(pat, line, conv=float):
    m = re.search(pat, line)
    return conv(m.group(1)) if m else None


def read_trace(path):
    rows = {}
    for line in open(path):
        if " fixes " not in line:
            continue
        fixes = frozenset(re.search(r"fixes ([\d: ]*?) \|", line).group(1).split())
        cut = num(r"\| cut (\d)", line, int) == 1
        rows[fixes] = dict(
            depth=num(r"depth (\d+)", line, int), lb=num(r"first lb (\S+)", line),
            status=num(r"status (\S+)", line, str), iters=num(r"iters (\d+)", line, int),
            cpu=num(r"cpu (\S+)", line), warm=num(r"warm (\d)", line, int), pool=num(r"pool (\d)", line, int),
            cut=cut, cut_best=num(r"\| cut 1 best (\S+)", line) if cut else None,
            cut_iters=num(r"rounds \d+ iters (\d+)", line, int) if cut else 0,
            cut_cpu=num(r"rounds \d+ iters \d+ cpu (\S+)", line) if cut else 0.0,
            cut_status=num(r"rounds \d+ iters \d+ cpu \S+ status (\S+)", line, str) if cut else None,
            fathomed=num(r"fathomed (\d) ", line, int), cutoff=num(r"cutoff (\S+)", line))
    return rows


def read_arm(d):
    arm = os.path.basename(d)
    out = {"arm": arm, "dir": d, "done": os.path.exists(os.path.join(d, "results.txt"))}
    if not out["done"]:
        return out
    out["rows"] = rows = read_trace(os.path.join(d, "trace.txt"))
    res = open(os.path.join(d, "results.txt")).read()
    m = re.search(r"coordinator: (\d+) slots, (\d+) units, wall (\S+) s, CPU (\S+) s .*slots busy (\S+)%, "
                  r"best objective (\S+)", res)
    if m:
        out.update(slots=int(m.group(1)), units=int(m.group(2)), wall=float(m.group(3)), cpu=float(m.group(4)),
                   busy=float(m.group(5)), best=m.group(6))
    m = re.search(r"Nodes: (\d+) / (\d+)", res)
    if m:
        out.update(nodes=int(m.group(1)), fathomed=int(m.group(2)))
    out["first_it"] = sum(r["iters"] for r in rows.values())
    out["cut_it"] = sum(r["cut_iters"] for r in rows.values())
    out["dnn_it"] = out["first_it"] + out["cut_it"]
    out["dnn_cpu"] = sum(r["cpu"] + r["cut_cpu"] for r in rows.values())
    out["phases"] = sum(r["cut"] for r in rows.values())
    out["phase_fathoms"] = sum(r["cut"] and r["fathomed"] == 1 for r in rows.values())
    out["dnn_fathomed"] = sum(r["fathomed"] == 1 for r in rows.values())
    out["warm"] = sum(r["warm"] == 1 for r in rows.values())
    out["pooled"] = sum(r["pool"] == 1 for r in rows.values())
    out["hosts"] = open(os.path.join(d, "hosts.txt")).read().strip().splitlines() \
        if os.path.exists(os.path.join(d, "hosts.txt")) else []
    return out


def ratio(a, b):
    return f"{a / b:.3f}x" if a and b else "-"


def compare(x, y):
    """x against y (y the reference)"""
    print(f"\n## {x['arm']} against {y['arm']}")
    print(f"- tree: {len(x['rows'])} DNN nodes against {len(y['rows'])}; "
          f"cut phases {x['phases']} ({x['phase_fathoms']} fathomed) against {y['phases']} ({y['phase_fathoms']})")
    print(f"- DNN iterations {y['dnn_it']:,} -> {x['dnn_it']:,}: **{ratio(y['dnn_it'], x['dnn_it'])}** "
          f"(first solves {y['first_it']:,} -> {x['first_it']:,}, cut phases {y['cut_it']:,} -> {x['cut_it']:,})")
    print(f"- DNN CPU {y['dnn_cpu']:,.0f} -> {x['dnn_cpu']:,.0f} s: {ratio(y['dnn_cpu'], x['dnn_cpu'])}; "
          f"coordinator CPU {y.get('cpu', 0):,.0f} -> {x.get('cpu', 0):,.0f} s: {ratio(y.get('cpu'), x.get('cpu'))}; "
          f"wall {y.get('wall', 0):,.0f} -> {x.get('wall', 0):,.0f} s: {ratio(y.get('wall'), x.get('wall'))}")
    shared = set(x["rows"]) & set(y["rows"])
    if shared:
        ix = sum(x["rows"][k]["iters"] + x["rows"][k]["cut_iters"] for k in shared)
        iy = sum(y["rows"][k]["iters"] + y["rows"][k]["cut_iters"] for k in shared)
        print(f"- the {len(shared)} shared nodes: {iy:,} -> {ix:,} DNN iterations, {ratio(iy, ix)}; "
              f"{len(set(x['rows']) - shared)} nodes only in {x['arm']}, {len(set(y['rows']) - shared)} only in {y['arm']}")
        flips = [k for k in shared if x["rows"][k]["fathomed"] != y["rows"][k]["fathomed"]]
        if flips:
            print(f"- outcome changed on {len(flips)} shared nodes:")
            for k in sorted(flips, key=lambda k: (x["rows"][k]["depth"], sorted(k))):
                rx, ry = x["rows"][k], y["rows"][k]
                desc = lambda r: (f"{'fathomed' if r['fathomed'] else 'branched'} (first {r['iters']} it, lb {r['lb']:.1f}"
                                  + (f"; phase {r['cut_iters']} it, best {r['cut_best']:.1f}" if r["cut"] else "") + ")")
                print(f"  - depth {rx['depth']} [{' '.join(sorted(k))}]: {y['arm']} {desc(ry)}, {x['arm']} {desc(rx)}")
    it = lambda r: r["iters"] + r["cut_iters"]
    matched = [k for k in shared if x["rows"][k]["warm"] == y["rows"][k]["warm"]
               and x["rows"][k]["pool"] == y["rows"][k]["pool"]]
    rm = None
    if matched:
        mx, my = sum(it(x["rows"][k]) for k in matched), sum(it(y["rows"][k]) for k in matched)
        rm = my / mx if mx else 0
        print(f"- the {len(matched)} matched nodes (shared, same warm start and pool status): {my:,} -> {mx:,} "
              f"DNN iterations, {ratio(my, mx)}")
    both = [k for k in shared if x["rows"][k]["cut"] and y["rows"][k]["cut"]]
    if both:
        px, py = sum(x["rows"][k]["cut_iters"] for k in both), sum(y["rows"][k]["cut_iters"] for k in both)
        fx = sum(x["rows"][k]["fathomed"] == 1 for k in both); fy = sum(y["rows"][k]["fathomed"] == 1 for k in both)
        gx = sum(x["rows"][k]["cut_best"] - x["rows"][k]["lb"] for k in both) / len(both)
        gy = sum(y["rows"][k]["cut_best"] - y["rows"][k]["lb"] for k in both) / len(both)
        print(f"- the {len(both)} nodes with a phase in both: phase iterations {py:,} -> {px:,} ({ratio(py, px)}), "
              f"fathomed by the phase {fy} -> {fx}, mean gain over the first solve {gy:+.1f} -> {gx:+.1f}")
    depths = sorted({r["depth"] for r in list(x["rows"].values()) + list(y["rows"].values())})
    print(f"- by depth (nodes, DNN iterations), {y['arm']} -> {x['arm']}: " + "; ".join(
        f"d{d} {sum(r['depth'] == d for r in y['rows'].values())}/{sum(it(r) for r in y['rows'].values() if r['depth'] == d):,}"
        f" -> {sum(r['depth'] == d for r in x['rows'].values())}/{sum(it(r) for r in x['rows'].values() if r['depth'] == d):,}"
        for d in depths))
    gate = GATES.get((x["arm"], y["arm"]))
    if gate:
        r = y["dnn_it"] / x["dnn_it"] if x["dnn_it"] else 0
        if x["arm"] in ("zone", "ship", "full"):
            print(f"- kill line: {gate:.2f}x in DNN iterations end to end -> {'PASS' if r >= gate else 'FAIL'} ({r:.3f}x)")
        elif rm is not None:
            print(f"- kill line: {gate:.2f}x in DNN iterations on the matched nodes -> "
                  f"{'PASS' if rm >= gate else 'FAIL'} ({rm:.3f}x; end to end {r:.3f}x)")


def main():
    runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
    arms = {}
    for a in ARMS:
        d = os.path.join(runs, a)
        if os.path.isdir(d):
            arms[a] = read_arm(d)
    print("# pipeline A/B on the node-10 subtree\n")
    print("| arm | nodes / fathomed | DNN nodes / DNN-fathomed | cut phases (fathomed) | warm / pooled nodes "
          "| DNN iterations (first + phase) | DNN CPU s | CPU s | wall s | busy | best |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for a, x in arms.items():
        if not x["done"]:
            print(f"| {a} | not finished | | | | | | | | | |")
            continue
        print(f"| {a} | {x.get('nodes', '?')} / {x.get('fathomed', '?')} | {len(x['rows'])} / {x['dnn_fathomed']} "
              f"| {x['phases']} ({x['phase_fathoms']}) | {x['warm']} / {x['pooled']} "
              f"| {x['dnn_it']:,} ({x['first_it']:,} + {x['cut_it']:,}) | {x['dnn_cpu']:,.0f} | {x.get('cpu', 0):,.0f} "
              f"| {x.get('wall', 0):,.0f} | {x.get('busy', 0):.1f}% | {x.get('best', '?')} |")
    for a, x in arms.items():
        if x["done"]:
            for h in x["hosts"]:
                print(f"  {a}: {h}")
    done = {a: x for a, x in arms.items() if x["done"]}
    for a, refs in PAIRS:
        ref = next((r for r in refs if r in done), None)
        if a in done and ref:
            compare(done[a], done[ref])
    # the -F node itself: did it run a phase (the equip phase), and what did it hand on
    for a, x in done.items():
        top = min(x["rows"].values(), key=lambda r: r["depth"]) if x["rows"] else None
        if top:
            print(f"\n{a}: the subtree's root (depth {top['depth']}): first lb {top['lb']:.2f} after {top['iters']} it"
                  + (f", phase best {top['cut_best']:.2f} in {top['cut_iters']} it ({top['cut_status']})" if top["cut"] else ", no phase"))


main()
