#!/usr/bin/env python3
"""Validation of the stage-A launcher and the stage-B restart path on a small
instance (README.md, "Validation"):

1. the parallel stage-A frontier (stage_a.py, one task per node) is the same
   node set, with the same inherited and predicted bounds, as one serial
   `bnbdrv -d D -o` run with the same binary and environment, and the per-level
   statistics and the per-node DNN results agree;
2. a frontier node restarted the stage-B way (a one-node state file, -r) gets
   the same first DNN bound as the same node started with -F, when the restart
   carries no warm start (-F has none); with the parent's warm start, as stage
   B runs it, the first bound is printed next to it.

  validate.py --instance qaplib/nug20.dat --sym latest/SYM/nug20.sym --ub 2570 \\
      --grid 4:5 --depth 2 --ship 64 --tag nug20_ship64 [--restart-nodes 3]

--ship 64 ships the parent's DNN state as exact doubles (a resumed run is then
the same run); --ship 32 is the production setting.
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (D, PARAMS, R, XSPEC, State, fullrun_env, node_fields, parse_trace, private_binary,  # noqa: E402
                    run_bnbdrv, solver_env, spec, write_one)


def frontier_of(path):
    st = State(path)
    out = {}
    for text, w in st.nodes:
        d, b, pb, f = node_fields(text, st.n)
        out[spec(f)] = (d, b, pb)
    return st, out


def compare_levels(a, b, name):
    bad = []
    for k in ("n_level", "f_level", "l_level", "e_level"):
        if a.arrays[k] != b.arrays[k]:
            bad.append(f"{k}: serial {a.arrays[k][:4]} parallel {b.arrays[k][:4]}")
    if a.arrays["dnn_counts"] != b.arrays["dnn_counts"]:
        bad.append(f"dnn_counts: serial {a.arrays['dnn_counts']} parallel {b.arrays['dnn_counts']}")
    if a.counters[:3] != b.counters[:3]:
        bad.append(f"counters: serial {a.counters} parallel {b.counters}")
    print(f"{name}: per-level statistics " + ("agree" if not bad else "DIFFER:\n  " + "\n  ".join(bad)))
    return not bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instance", required=True)
    ap.add_argument("--sym", required=True)
    ap.add_argument("--ub", type=float, required=True)
    ap.add_argument("--grid", required=True, help="QAP_DNN_CUT_GRID of the instance (side 1)")
    ap.add_argument("--depth", type=int, default=2)
    ap.add_argument("--ship", default="32")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--restart-nodes", type=int, default=3)
    ap.add_argument("--env", action="append", default=[], help="K=V on top (diagnostics)")
    a = ap.parse_args()
    inst, sym = os.path.abspath(a.instance), os.path.abspath(a.sym)
    vdir = os.path.join(D, "runs", "val_" + a.tag)
    os.makedirs(vdir, exist_ok=True)
    overrides = {"QAP_DNN_CUT_GRID": a.grid, "QAP_DNN_SHIP": a.ship}
    for kv in a.env:
        k, v = kv.split("=", 1)
        overrides[k] = v
    report = {"tag": a.tag, "instance": a.instance, "depth": a.depth, "ship": a.ship}

    # 1. the parallel stage A
    run = os.path.join(vdir, "stageA")
    cmd = [sys.executable, os.path.join(D, "stage_a.py"), "--run", run, "--instance", inst, "--sym", sym,
           "--ub", repr(a.ub), "--depth", str(a.depth), "--slots", str(a.slots), "--tau", "1.0"]
    for k, v in overrides.items():
        cmd += ["--env", f"{k}={v}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"stage_a.py failed ({r.returncode}):\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    binary, commit = private_binary(os.path.join(R, "latest", "bnbdrv"), os.path.join(run, "bin"))
    settings, _ = fullrun_env(binary)
    env = solver_env(settings, overrides)

    # 2. the serial run, same binary and environment
    ser = os.path.join(vdir, "serial.state")
    strace = os.path.join(vdir, "serial.trace")
    if not os.path.exists(ser):
        rc, wall = run_bnbdrv(binary, ["-p", PARAMS, "-x", XSPEC, "-s", sym, "-u", repr(a.ub), "-d", str(a.depth),
                                       "-o", ser, inst], env, os.path.join(vdir, "serial.log"), strace, cpu=0)
        if rc != 0:
            sys.exit(f"serial run failed: {rc}")
    sst, sfr = frontier_of(ser)
    pst, pfr = frontier_of(os.path.join(run, "frontier.state"))
    same_set = set(sfr) == set(pfr)
    print(f"{a.tag}: frontier serial {len(sfr)} nodes, parallel {len(pfr)} nodes; same node set: {same_set}")
    if not same_set:
        print("  only serial:", sorted(set(sfr) - set(pfr))[:10])
        print("  only parallel:", sorted(set(pfr) - set(sfr))[:10])
    common_ = set(sfr) & set(pfr)
    db = max((abs(sfr[k][1] - pfr[k][1]) for k in common_), default=0.0)
    dp = max((abs(sfr[k][2] - pfr[k][2]) for k in common_), default=0.0)
    print(f"  inherited bounds differ by at most {db:.3g}, predicted bounds by at most {dp:.3g}")
    levels_ok = compare_levels(sst, pst, a.tag)
    # per-node DNN results
    srows = {x["fixes"]: x for x in parse_trace(strace)}
    prows = {}
    for tdir in os.listdir(os.path.join(run, "tasks")):
        for x in parse_trace(os.path.join(run, "tasks", tdir, "trace.txt")):
            prows[x["fixes"]] = x
    keys = set(srows) & set(prows)
    dlb = max((abs(srows[k]["lb"] - prows[k]["lb"]) for k in keys), default=0.0)
    dit = sum(srows[k]["iters"] != prows[k]["iters"] for k in keys)
    dst = sum(srows[k]["status"] != prows[k]["status"] for k in keys)
    warm_s = sum(x["warm"] for x in srows.values())
    warm_p = sum(x["warm"] for x in prows.values())
    print(f"  DNN nodes: serial {len(srows)}, parallel {len(prows)}, shared {len(keys)}; first lb differs by at most "
          f"{dlb:.3g}; iterations differ on {dit}, status on {dst}; warm-started serial {warm_s}, parallel {warm_p}")
    report.update(frontier_serial=len(sfr), frontier_parallel=len(pfr), same_set=same_set, max_bound_diff=db,
                  max_pred_diff=dp, levels_agree=levels_ok, dnn_nodes_serial=len(srows), dnn_nodes_parallel=len(prows),
                  max_first_lb_diff=dlb, iteration_mismatches=dit, status_mismatches=dst, warm_serial=warm_s,
                  warm_parallel=warm_p, cpu_serial=sst.times[0], cpu_parallel=pst.times[0])

    # 3. the stage-B restart against -F, on a few frontier nodes (spread over parents)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(run, "frontier.tsv"))][1:]
    rows.sort(key=lambda x: float(x[5]))
    picked = rows[:a.restart_nodes]
    report["restart"] = []
    for parent, idx, sp, *_ in picked:
        res = {"spec": sp}
        tdir = os.path.join(run, "tasks", parent)
        src = os.path.join(tdir, "out.full.state")
        if not os.path.exists(src):
            src = os.path.join(tdir, "out.slim.state")
        for mode in ("none", "keep"):
            tag = f"rs_{sp.replace(':', '_').replace(',', '-')}_{mode}"
            inp = os.path.join(vdir, tag + ".in.state")
            write_one(src, int(idx), inp, warm=mode)
            rc, _ = run_bnbdrv(binary, ["-p", PARAMS, "-x", XSPEC, "-s", sym, "-u", repr(a.ub),
                                        "-d", str(a.depth + 1), "-r", inp, "-o", os.path.join(vdir, tag + ".out.state"),
                                        inst], env, os.path.join(vdir, tag + ".log"),
                               os.path.join(vdir, tag + ".trace"), cpu=0)
            tr = parse_trace(os.path.join(vdir, tag + ".trace"))
            res[mode] = tr[0] if tr else None
        tag = f"rs_{sp.replace(':', '_').replace(',', '-')}_F"
        rc, _ = run_bnbdrv(binary, ["-p", PARAMS, "-x", XSPEC, "-s", sym, "-u", repr(a.ub), "-d", str(a.depth + 1),
                                    "-F", sp, "-o", os.path.join(vdir, tag + ".out.state"), inst], env,
                           os.path.join(vdir, tag + ".log"), os.path.join(vdir, tag + ".trace"), cpu=0)
        tr = parse_trace(os.path.join(vdir, tag + ".trace"))
        res["F"] = tr[0] if tr else None
        f = lambda x: "QPB-fathomed" if x is None else f"{x['lb']:.6f} ({x['status']}, {x['iters']} it, warm {x['warm']})"
        same = (res["none"] is None and res["F"] is None) or (
            res["none"] is not None and res["F"] is not None and res["none"]["lb"] == res["F"]["lb"]
            and res["none"]["iters"] == res["F"]["iters"])
        res["none_equals_F"] = same
        print(f"  restart {sp}: -F {f(res['F'])}; stage-B path without warm start {f(res['none'])} "
              f"[{'same' if same else 'DIFFERENT'}]; with the parent's warm start {f(res['keep'])}")
        report["restart"].append(res)
    report["commit"] = commit
    with open(os.path.join(vdir, "report.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    ok = same_set and levels_ok and all(x["none_equals_F"] for x in report["restart"])
    print(f"{a.tag}: {'PASS' if ok else 'CHECK'} (report {os.path.join(vdir, 'report.json')})")


if __name__ == "__main__":
    main()
