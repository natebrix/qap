#!/usr/bin/env python3
"""Phase 2 gate, cloud box: wil50 stage-A tasks replayed from their kept parent
states, the stage-A binary (a667727, stage A's own environment and command
line) against the new binary (HEAD, fullrun.par through the sampling code path:
fullrun_env + params_with + write_one + run_bnbdrv, as stage_a.py runs a task).

  gate.py run JOB[,JOB...]     (JOB = side:task, side new|old)
"""
import filecmp
import json
import os
import queue
import sys
import threading
import time

W = "/home/user/qap/.claude/worktrees/agent-aebd90129de6196e4"
SA = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA"
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
sys.path.insert(0, os.path.join(W, "experiments", "wil50_feasibility", "sampling"))
from common import (PARAMS, XSPEC, dnn_flags, fullrun_env, params_with, private_binary, run_bnbdrv,  # noqa: E402
                    solver_env, write_one)

import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("parallel_a667727", os.path.join(G, "old_py", "parallel.py"))
old_parallel = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(old_parallel)


def old_write_one(src, index, out):
    """a667727's write_one (sampling/common.py at a667727) with a667727's State"""
    st = old_parallel.State(src)
    assert st.magic[1] == "6", st.magic
    st.zero_stats()
    st.write(out, [st.nodes[index]])


def task_info(tid):
    d = json.load(open(os.path.join(SA, "tasks", tid, "done.json")))
    return d


def prep(side, tid, wd):
    """the task's one-node input state, as stage A built it (parent's out.state, node index)"""
    d = task_info(tid)
    parent, index = d["parent"], d["index"]
    src = os.path.join(SA, "tasks", parent, "out.state")
    inp = os.path.join(wd, "in.state")
    write_one(src, index, inp, warm="keep")
    chk = os.path.join(wd, "in.a667727.state")
    old_write_one(src, index, chk)
    same = filecmp.cmp(inp, chk, shallow=False)
    os.remove(chk)
    return d, inp, same


def main():
    jobs = [tuple(j.split(":")) for j in sys.argv[2].split(",")]
    cpus = [int(c) for c in os.environ.get("GATE_CPUS", "0,1,2,3").split(",")]
    # the new side: the sampling launchers' path (stage_a.py's StageA.__init__ and _run_task)
    nbin, ncommit = private_binary(os.path.join(W, "latest", "bnbdrv"), os.path.join(G, "bin_new"))
    nsettings, info = fullrun_env(nbin)
    npar = os.path.join(G, "params.par")
    nenv = solver_env(nsettings, params_with(PARAMS, {}, npar, nbin))
    new = {"binary": nbin, "env": nenv, "args": ["-p", npar] + dnn_flags(XSPEC), "sym": info["FULLRUN_SYM"],
           "ub": float(info["FULLRUN_UB"]), "instance": info["FULLRUN_INSTANCE"]}
    # the old side: stage A's own record (config.json: the settings environment, the flags)
    cfg = json.load(open(os.path.join(SA, "config.json")))
    oenv = solver_env(cfg["settings_env"], cfg["env_overrides"])
    old = {"binary": os.path.join(SA, "bin", "bnbdrv"), "env": oenv,
           "args": ["-p", cfg["params"]] + dnn_flags(cfg["xspec"]), "sym": cfg["sym"], "ub": cfg["ub"],
           "instance": cfg["instance"]}
    sides = {"new": new, "old": old}
    # attribution: 18ebd58's parent (347ff84) and 18ebd58 itself, run exactly as the old side
    for name, d in (("pre", "b_pre"), ("fix", "b_fix")):
        sides[name] = dict(old, binary=os.path.join(G, d, "latest", "bnbdrv"))
    with open(os.path.join(G, "sides.json"), "w") as f:
        json.dump({k: {kk: vv for kk, vv in v.items()} for k, v in sides.items()} | {"new_commit": ncommit}, f,
                  indent=1, sort_keys=True)
    q = queue.Queue()
    for j in jobs:
        q.put(j)
    lock = threading.Lock()

    def slot(cpu):
        while True:
            try:
                side, tid = q.get_nowait()
            except queue.Empty:
                return
            s = sides[side]
            wd = os.path.join(G, "runs", side, tid)
            os.makedirs(wd, exist_ok=True)
            out = os.path.join(wd, "out.state")
            if tid == "r":      # the root task: no input state (stage_a.py: fixes [], no -r)
                d, inp, same = task_info(tid), None, None
                resume = []
            else:
                d, inp, same = prep(side, tid, wd)
                resume = ["-r", inp]
            args = s["args"] + ["-s", s["sym"], "-u", repr(s["ub"]), "-d", str(d["depth"]), "-o", out] + resume + [
                s["instance"]]
            with lock:
                print(f"{time.strftime('%F %T')} start {side} {tid} cpu {cpu} (in.state same as a667727's writer: "
                      f"{same})", flush=True)
            rc, wall = run_bnbdrv(s["binary"], args, s["env"], os.path.join(wd, "solver.log"),
                                  os.path.join(wd, "trace.txt"), cpu=cpu)
            if inp:
                os.remove(inp)
            with open(os.path.join(wd, "done.json"), "w") as f:
                json.dump({"side": side, "task": tid, "rc": rc, "wall": wall, "cpu_slot": cpu,
                           "in_state_same_as_a667727_writer": same, "finished": time.strftime("%F %T")}, f)
            with lock:
                print(f"{time.strftime('%F %T')} done {side} {tid} rc {rc} wall {wall:.0f} s", flush=True)

    th = [threading.Thread(target=slot, args=(c,)) for c in cpus[:min(len(cpus), len(jobs))]]
    for t in th:
        t.start()
    for t in th:
        t.join()
    print(f"{time.strftime('%F %T')} all done", flush=True)


if __name__ == "__main__":
    main()
