#!/usr/bin/env python3
"""A layered, sampled probe of one depth-3 wil50 subtree (README.md).

The subtrees of the depth-3 nodes 329-405 short (node 1 of tuning/lines,
34:40,28:22,14:29, 455 short at 300 iterations, 405 converged) carry 85-91% of
the wil50 estimate (roadmap/tree_size.md) and are far too large to finish on a
laptop. So the probe bounds the node and every one of its children, then a
seeded sample of the next levels, and estimates each level from the sample:

  level 3   the node itself                       (-F SPEC -d 3)
  level 4   all its children                      (-r task -d 4, one task each)
  level 5   a sample of --sample5 grandchildren   (-r task -d 5)
  level 6   a sample of --sample6 of their children
  deeper    estimated as open nodes, not bounded

Every task is one node in its own state file, run by its own solver process,
so children start from their parent's DNN state only when the arm carries it
(QAP_DNN_SHIP): that is how a fleet run splits work, and how the arms differ.

  python3 tail_probe.py ARM [--slots N] [--sample5 30] [--sample6 20] [--seed 1]
  ARM: base (the Nov 4 production setting; every node starts cold, since each
  is its own task), baseship (base + the shipped parent state: every node
  warm, as in a serial run), full (+ cut zone 0.005, warm-Y rounds with
  damping, shipped state), fullwide (full with cut zone 0.010, so that node 1
  itself runs a cut phase and its children inherit a pool)

Resumable: a finished task (out.state and trace present) is not rerun.
Results: runs/ARM/ (tasks, traces) and runs/ARM/summary.txt.
"""
import argparse
import concurrent.futures as cf
import os
import platform
import random
import re
import subprocess
import sys
import time

D = os.path.dirname(os.path.abspath(__file__))
R = os.path.abspath(os.path.join(D, "../../.."))
sys.path.insert(0, os.path.join(R, "bench"))
from parallel import State  # noqa: E402

BASE_ENV = {
    "QAP_DNN_GAMMA1": "1.0", "QAP_DNN_GAMMA": "1.0", "QAP_DNN_EXT": "20",
    "QAP_DNN_CUT_ZONE": "0.003", "QAP_DNN_CUT_GRID": "5:10", "QAP_DNN_CUT_SIDE": "0", "QAP_DNN_CUT_FREEZE": "1500",
    "DNN_AMX": "1", "DNN_CHOLQR_PASSES": "1", "QAP_DNN_CHEB": "3", "QAP_DNN_TIGHT_LAM": "2", "QAP_DNN_REFRESH": "40",
    "QAP_DNN_WARM_Y": "2", "QAP_DNN_POOL": "1", "QAP_DNN_POOL_CHECK": "20",
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
FULL = {"QAP_DNN_CUT_ZONE": "0.005", "QAP_DNN_CUT_WARMY": "1", "QAP_DNN_CUT_ROUNDS": "12", "QAP_DNN_CUT_ITERS": "200",
        "QAP_DNN_CUT_FREEZE": "1000", "QAP_DNN_CUT_DAMP1": "0.3", "QAP_DNN_SHIP": "16"}
# bf16 shipping (16), not the production fp32 (32): every task file carries its
# parent's state, about 13 MB at m = 46 instead of 26, and on 8 children of
# node 10 bf16, fp32 and exact gave the same iterations (../ship/README.md)
ARMS = {"base": {}, "baseship": {"QAP_DNN_SHIP": "16"}, "full": FULL,
        "fullwide": dict(FULL, QAP_DNN_CUT_ZONE="0.010"),
        # equip3: a phase starting more than 120 short runs at most 3 rounds (README, "Equip phases")
        "equip3": dict(FULL, QAP_DNN_CUT_EQUIP="120:3"),
        # the two level-4 fathoms equip3 lost started 127-128 short: cap at 150
        "equip150": dict(FULL, QAP_DNN_CUT_EQUIP="150:3"),
        "equip150r2": dict(FULL, QAP_DNN_CUT_EQUIP="150:2")}


def run_task(args, env, path_in, level, tag, wd, spec=None):
    """bound the node(s) of path_in (or the -F node) down to `level`; children open in out"""
    out = os.path.join(wd, tag + ".out.state")
    trace = os.path.join(wd, tag + ".trace")
    log = os.path.join(wd, tag + ".log")
    if os.path.exists(out) and os.path.exists(trace):
        return out, trace
    e = dict(os.environ)
    e.update(env)
    e["QAP_DNN_TRACE"] = trace + ".tmp"
    if os.path.exists(trace + ".tmp"):
        os.remove(trace + ".tmp")
    cmd = [args.solver, "-p", args.params] + args.extra.split() + ["-u", str(args.ub), "-d", str(level), "-o", out]
    cmd += (["-F", spec] if spec else ["-r", path_in]) + [os.path.join(R, args.instance)]
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, cwd=os.path.join(R, "latest"), env=e, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0 or not os.path.exists(out):
        raise RuntimeError(f"{tag}: exit {rc} (see {log})")
    if not os.path.exists(trace + ".tmp"):
        open(trace + ".tmp", "w").close()      # a node QPB fathomed: no DNN line
    os.replace(trace + ".tmp", trace)
    if path_in and os.path.exists(path_in):
        os.remove(path_in)                     # disk: the input's warm line is the parent's state
    return out, trace


def open_count(path):
    """the number of open nodes in a state file, without parsing the warm lines
    (with shipped DNN state a file holds tens of MB)"""
    with open(path) as f:
        for line in f:
            if line.startswith("nodes "):
                return int(line.split()[1])
    return 0


def write_tasks(picks, wd, lv):
    """one state file per picked (src, i): node i of src alone, statistics zeroed,
    with its warm start line; each src parsed once"""
    out = {}
    by = {}
    for j, (src, i) in enumerate(picks):
        by.setdefault(src, []).append((j, i))
    for src, items in by.items():
        todo = [(j, i) for j, i in items if not os.path.exists(os.path.join(wd, f"L{lv}_{j}.in.state"))]
        st = None
        if todo:
            st = State(src)
            st.zero_stats()
        for j, i in items:
            p = os.path.join(wd, f"L{lv}_{j}.in.state")
            if not os.path.exists(p):
                st.write(p, [st.nodes[i]])
            out[j] = p
        del st
    return [(out[j], f"L{lv}_{j}") for j in range(len(picks))]


def trace_rows(paths):
    rows = []
    for p in paths:
        for line in open(p):
            if " fixes " not in line:
                continue
            g = lambda pat, conv=float: (lambda m: conv(m.group(1)) if m else None)(re.search(pat, line))
            cut = g(r"\| cut (\d)", int) == 1
            rows.append(dict(depth=g(r"depth (\d+)", int), lb=g(r"first lb (\S+)"), iters=g(r"iters (\d+)", int),
                             cpu=g(r"cpu (\S+)"), status=g(r"status (\S+)", str), warm=g(r"warm (\d)", int),
                             pool=g(r"pool (\d)", int), cut=cut,
                             cut_best=g(r"\| cut 1 best (\S+)") if cut else None,
                             cut_iters=g(r"rounds \d+ iters (\d+)", int) if cut else 0,
                             cut_cpu=g(r"rounds \d+ iters \d+ cpu (\S+)") if cut else 0.0,
                             fathomed=g(r"fathomed (\d) ", int), cutoff=g(r"cutoff (\S+)")))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arm", choices=sorted(ARMS))
    ap.add_argument("--node", default="34:40,28:22,14:29", help="the -F spec of the depth-3 node (default node 1)")
    ap.add_argument("--level", type=int, default=3, help="its depth")
    ap.add_argument("--instance", default="qaplib/wil50.dat")
    ap.add_argument("--ub", type=float, default=48816)
    ap.add_argument("--params", default="PARAM/nug_deep_xbar9_all.par")
    ap.add_argument("--extra", default="-x 0:4:1:300:1")
    ap.add_argument("--grid", default="5:10")
    ap.add_argument("--slots", type=int, default=None)
    ap.add_argument("--sample5", type=int, default=30)
    ap.add_argument("--sample6", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--runs", default=os.path.join(D, "runs"))
    ap.add_argument("--solver", default=os.path.join(R, "latest", "bnbdrv"))
    args = ap.parse_args()
    if args.slots is None:
        n = os.cpu_count() or 2
        if platform.system() == "Darwin":
            try:
                n = int(subprocess.run(["sysctl", "-n", "hw.perflevel0.physicalcpu"], capture_output=True,
                                       text=True).stdout.strip())
            except ValueError:
                pass
        args.slots = max(2, n // 2)
    env = dict(BASE_ENV, **ARMS[args.arm])
    env["QAP_DNN_CUT_GRID"] = args.grid
    wd = os.path.join(args.runs, args.arm)
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "env.txt"), "w") as f:
        f.write(" ".join(f"{k}={v}" for k, v in sorted(env.items())) + "\n" + " ".join(sys.argv) + "\n")
    log = lambda m: print(f"{time.strftime('%H:%M:%S')} {args.arm}: {m}", flush=True)
    rng = random.Random(args.seed)
    L0 = args.level
    levels = {}                      # level -> dict(n_est, sampled traces, open children per sampled node)
    t0 = time.time()

    # the node itself
    out, tr = run_task(args, env, None, L0, "L%d_0" % L0, wd, spec=args.node)
    nk = open_count(out)
    levels[L0] = dict(n_est=1.0, traces=[tr], children=[nk], outs=[out])
    log(f"level {L0}: 1 node, {nk} children open")

    # level L0+1: every child; then samples
    frontier = [(out, nk)]
    n_parent_est = 1.0
    for lv, k in ((L0 + 1, None), (L0 + 2, args.sample5), (L0 + 3, args.sample6)):
        pool = [(src, i) for src, n in frontier for i in range(n)]
        if not pool:
            log(f"level {lv}: nothing open; the subtree is finished above it")
            break
        if k is not None and len(pool) > k:
            # stratified: round robin over parents, seeded
            by = {}
            for src, i in pool:
                by.setdefault(src, []).append((src, i))
            for v in by.values():
                rng.shuffle(v)
            order = [x for grp in zip(*[v + [None] * (max(map(len, by.values())) - len(v)) for v in by.values()])
                     for x in grp if x is not None]
            pick = order[:k]
        else:
            pick = pool
        tasks = write_tasks(pick, wd, lv)
        log(f"level {lv}: {len(pool)} nodes open, bounding {len(tasks)} on {args.slots} slots")
        res = {}
        with cf.ThreadPoolExecutor(max_workers=args.slots) as ex:
            futs = {ex.submit(run_task, args, env, p, lv, tag, wd): tag for p, tag in tasks}
            for f in cf.as_completed(futs):
                res[futs[f]] = f.result()
        outs = [res[tag][0] for _, tag in tasks]
        traces = [res[tag][1] for _, tag in tasks]
        children = [open_count(o) for o in outs]
        n_est = n_parent_est * (sum(levels[lv - 1]["children"]) / len(levels[lv - 1]["children"]))
        levels[lv] = dict(n_est=n_est, traces=traces, children=children, outs=outs, open=len(pool), bounded=len(tasks))
        n_parent_est = n_est
        frontier = [(o, c) for o, c in zip(outs, children)]
        log(f"level {lv}: {sum(c > 0 for c in children)} of {len(tasks)} branched, {sum(children)} children")

    # summary
    lines = [f"# tail probe, arm {args.arm}: node {args.node} (depth {L0}), {args.instance}, ub {args.ub:g}",
             f"# {platform.node()} {platform.machine()}, {args.slots} slots, wall {time.time() - t0:.0f} s (this session)",
             "", "| level | nodes (est.) | bounded | fathomed | branched | DNN iterations per node (first + phase) | "
             "CPU s per node | est. iterations at the level | est. CPU s at the level | warm | pooled | phases |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    tot_it = tot_cpu = 0.0
    last = max(levels)
    for lv in sorted(levels):
        L = levels[lv]
        rows = trace_rows(L["traces"])
        b = len(L["traces"])
        it = sum(r["iters"] + r["cut_iters"] for r in rows) / b
        cpu = sum(r["cpu"] + r["cut_cpu"] for r in rows) / b
        fi = sum(r["iters"] for r in rows) / b
        ci = sum(r["cut_iters"] for r in rows) / b
        br = sum(c > 0 for c in L["children"])
        tot_it += L["n_est"] * it
        tot_cpu += L["n_est"] * cpu
        lines.append(f"| {lv} | {L['n_est']:,.0f} | {b} | {b - br} | {br} | {it:,.0f} ({fi:,.0f} + {ci:,.0f}) | {cpu:,.1f} | "
                     f"{L['n_est'] * it:,.0f} | {L['n_est'] * cpu:,.0f} | {sum(r['warm'] == 1 for r in rows)} | "
                     f"{sum(r['pool'] == 1 for r in rows)} | {sum(r['cut'] for r in rows)} |")
    below = levels[last]["n_est"] * sum(levels[last]["children"]) / len(levels[last]["children"])
    lines += ["", f"Estimated through level {last}: {tot_it:,.0f} DNN iterations, {tot_cpu:,.0f} CPU s. "
              f"Open below level {last} (estimated): {below:,.0f} nodes at level {last + 1}, not bounded.", "",
              "Shortfalls of the branched nodes (cutoff minus the best certified value):"]
    for lv in sorted(levels):
        rows = trace_rows(levels[lv]["traces"])
        sh = sorted(round(r["cutoff"] - (r["cut_best"] if r["cut"] else r["lb"]), 1) for r in rows if r["fathomed"] == 0)
        lines.append(f"  level {lv}: {sh}")
    txt = "\n".join(lines) + "\n"
    open(os.path.join(wd, "summary.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
