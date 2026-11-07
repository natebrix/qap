#!/usr/bin/env python3
"""Stage A of the wil50 sampling plan (README.md; roadmap/tree_size.md
section 3): bound every node at depths 0..D (default D = 2) of the full run's
search with the full run's settings (fullrun/settings.sh) and save the
depth-(D+1) frontier.

Every node is its own task: one bnbdrv process, pinned to one CPU, resumes a
state file holding the node alone with its parent's warm start (the QPB point,
the cut pool, and under QAP_DNN_SHIP=32 the parent's DNN state) and bounds it
with -d <its depth>, so its children stay open in the task's out.state. That
is how a fleet run hands a node to another process. It also keeps memory and
disk bounded: a depth-2 parent's DNN state is 28 MB in a state file and about
85 MB in memory, so one process to depth 2 (or a coordinator unit holding
many depth-2 nodes) would hold 635 of them, 18 GB on disk and 54 GB in memory.

  stage_a.py [--run DIR] [--slots 4] [--max-new N] [--env K=V ...]

Defaults: the instance, symmetry file and incumbent of fullrun/settings.sh,
PARAM/nug_deep_xbar9_all.par, -x 1:4:1:300:1 (common.py). Order: shallower
tasks first; at the last depth round-robin over parents, so the first tasks
are a spread sample for the cost projection, printed after every task.

Resumable: a task with a done.json is never rerun; run the same command again.
SIGTERM or SIGINT drains (no new task starts, running ones finish); a second
signal kills the running solvers, whose tasks rerun on the next start.

Frontier files (in --run):
  tasks/<id>/out.slim.state  a last-depth task's children, without DNN state
  tasks/<id>/out.full.state  the same with the parent's DNN state, kept only
                             when a child is a stage-B candidate (hash < tau)
  frontier.state, frontier.tsv   all frontier nodes, at the end
The stage-B sample is the frontier nodes with the smallest unit_hash(seed, .):
tau is set so that about 430 of the 29,580 expected nodes fall below it.
"""
import argparse
import json
import os
import signal
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (D, PARAMS, R, XSPEC, State, disk_free_gb, fullrun_env, log_checks, log_line,  # noqa: E402
                    node_fields, parse_trace, private_binary, run_bnbdrv, slim_copy, solver_env, spec,
                    stats_of, unit_hash, write_json, write_one)

SEED = 20261107
TAU = 0.01454        # 430 / 29,580


def child_line(own, kids):
    """the branching line from the children's assignments: 'row i' (facility i
    placed at each child's location), 'col j' (location j), or '?' with one child"""
    new = [sorted(set(k) - set(own)) for k in kids]
    new = [x[0] for x in new if len(x) == 1]
    if len(new) < 2:
        return "?" if new else "-"
    if len({i for i, _ in new}) == 1:
        return f"row {new[0][0]}"
    if len({j for _, j in new}) == 1:
        return f"col {new[0][1]}"
    return "mixed"


class StageA:
    def __init__(self, a):
        self.a = a
        self.run = os.path.abspath(a.run)
        os.makedirs(os.path.join(self.run, "tasks"), exist_ok=True)
        self.log_path = os.path.join(self.run, "launcher.log")
        cfg_path = os.path.join(self.run, "config.json")
        self.binary, self.commit = private_binary(a.binary, os.path.join(self.run, "bin"))
        settings, info = fullrun_env(self.binary)
        cfg = {"instance": a.instance or info["FULLRUN_INSTANCE"], "sym": a.sym or info["FULLRUN_SYM"],
               "ub": a.ub if a.ub is not None else float(info["FULLRUN_UB"]), "params": a.params,
               "xspec": a.xspec, "depth": a.depth, "seed": a.seed, "tau": a.tau,
               "env_overrides": dict(kv.split("=", 1) for kv in a.env), "settings_env": settings,
               "commit": self.commit}
        if os.path.exists(cfg_path):
            old = json.load(open(cfg_path))
            for k in ("instance", "sym", "ub", "params", "xspec", "depth", "seed", "tau", "env_overrides",
                      "settings_env"):
                if old[k] != cfg[k]:
                    raise SystemExit(f"{cfg_path}: {k} was {old[k]!r}, now {cfg[k]!r}; use another --run")
            cfg = old
        else:
            write_json(cfg_path, cfg)
        self.cfg = cfg
        self.env = solver_env(cfg["settings_env"], cfg["env_overrides"])
        write_json(os.path.join(self.run, "env.json"),
                   {k: v for k, v in sorted(self.env.items()) if k.startswith(("QAP_", "DNN_", "OMP_", "OPENBLAS_",
                                                                              "VECLIB_"))})
        self.cpus = [int(x) for x in a.cpus.split(",")][:a.slots]
        self.cond = threading.Condition()
        self.tasks = {}          # id -> {id, depth, parent, index}
        self.done = {}           # id -> done.json
        self.ready = []
        self.running = {}        # id -> cpu
        self.started_kids = {}   # parent id -> children started
        self.procs = {}
        self.new_started = 0
        self.draining = False
        self.killing = False
        self.failed = []
        self._load()

    # ------------------------------------------------------------ the task graph
    def tdir(self, tid):
        return os.path.join(self.run, "tasks", tid)

    def _add(self, tid, depth, parent, index):
        self.tasks[tid] = {"id": tid, "depth": depth, "parent": parent, "index": index}
        p = os.path.join(self.tdir(tid), "done.json")
        if os.path.exists(p):
            self.done[tid] = json.load(open(p))
            if depth < self.cfg["depth"]:
                for k in range(self.done[tid]["n_children"]):
                    self._add(f"{tid}.{k}", depth + 1, tid, k)
        else:
            self.ready.append(tid)

    def _load(self):
        self._add("r", 0, None, None)
        log_line(self.log_path, f"stage A: {len(self.done)} tasks done, {len(self.ready)} ready "
                                f"(run {self.run}, binary {self.commit['commit'][:10]})")

    def _pick(self):
        if not self.ready:
            return None

        def key(tid):
            t = self.tasks[tid]
            return (t["depth"], self.started_kids.get(t["parent"], 0), t["parent"] or "", t["index"] or 0)
        tid = min(self.ready, key=key)
        self.ready.remove(tid)
        t = self.tasks[tid]
        self.started_kids[t["parent"]] = self.started_kids.get(t["parent"], 0) + 1
        return tid

    # ------------------------------------------------------------ one task
    def _run_task(self, tid, cpu):
        t = self.tasks[tid]
        wd = self.tdir(tid)
        os.makedirs(wd, exist_ok=True)
        inp, out = os.path.join(wd, "in.state"), os.path.join(wd, "out.state")
        cfg = self.cfg
        args = ["-p", cfg["params"], "-x", cfg["xspec"], "-s", cfg["sym"], "-u", repr(cfg["ub"]),
                "-d", str(t["depth"]), "-o", out]
        if t["parent"] is None:
            fixes = []
        else:
            src = os.path.join(self.tdir(t["parent"]), "out.state")
            _, _, _, fixes = write_one(src, t["index"], inp, warm="keep")
            args += ["-r", inp]
        args.append(cfg["instance"])
        trace = os.path.join(wd, "trace.txt")
        logf = os.path.join(wd, "solver.log")
        rc, wall = run_bnbdrv(self.binary, args, self.env, logf, trace, cpu=cpu, procs=self.procs, tag=tid)
        if rc != 0 or not os.path.exists(out):
            raise RuntimeError(f"task {tid}: bnbdrv exit {rc} (see {logf})")
        st = State(out)
        kids = []
        for text, _ in st.nodes:
            d, b, pb, f = node_fields(text, st.n)
            kids.append({"depth": d, "bound": b, "pred": pb, "fixes": f})
        doc = {"id": tid, "depth": t["depth"], "parent": t["parent"], "index": t["index"], "fixes": fixes,
               "spec": spec(fixes), "cpu_slot": cpu, "wall": wall, "rc": rc, "finished": time.strftime("%F %T"),
               "stats": stats_of(st), "checks": log_checks(logf), "trace": parse_trace(trace),
               "n_children": len(kids), "line": child_line(fixes, [k["fixes"] for k in kids])}
        if any(k["depth"] != t["depth"] + 1 for k in kids):
            raise RuntimeError(f"task {tid}: children not all at depth {t['depth'] + 1}")
        if t["depth"] >= cfg["depth"]:
            # frontier nodes: a slim copy, and the full one only for stage-B candidates
            for k in kids:
                k["hash"] = unit_hash(cfg["seed"], k["fixes"])
            doc["frontier"] = [{"spec": spec(k["fixes"]), "bound": k["bound"], "pred": k["pred"], "hash": k["hash"]}
                               for k in kids]
            slim_copy(out, os.path.join(wd, "out.slim.state"))
            doc["retained"] = any(k["hash"] < cfg["tau"] for k in kids)
            if doc["retained"]:
                os.replace(out, os.path.join(wd, "out.full.state"))
            else:
                os.remove(out)
        else:
            doc["children"] = [{"spec": spec(k["fixes"]), "bound": k["bound"], "pred": k["pred"]} for k in kids]
        if os.path.exists(inp):
            os.remove(inp)
        write_json(os.path.join(wd, "done.json"), doc)
        return doc

    # ------------------------------------------------------------ the slots
    def _slot(self, cpu):
        while True:
            with self.cond:
                while True:
                    stop = self.draining or (self.a.max_new is not None and self.new_started >= self.a.max_new)
                    if not stop and self.ready and disk_free_gb(self.run) < self.a.min_free_gb:
                        log_line(self.log_path, f"stage A: only {disk_free_gb(self.run):.1f} GB free; draining")
                        self.draining = stop = True
                    if stop or (not self.ready and not self.running):
                        self.cond.notify_all()
                        return
                    if self.ready:
                        tid = self._pick()
                        self.running[tid] = cpu
                        self.new_started += 1
                        break
                    self.cond.wait(timeout=5)
            try:
                doc = self._run_task(tid, cpu)
                err = None
            except Exception as e:      # noqa: BLE001  (a failed task is reported and rerun next start)
                doc, err = None, e
            with self.cond:
                del self.running[tid]
                if doc is not None:
                    self.done[tid] = doc
                    t = self.tasks[tid]
                    if t["depth"] < self.cfg["depth"]:
                        for k in range(doc["n_children"]):
                            cid = f"{tid}.{k}"
                            self.tasks[cid] = {"id": cid, "depth": t["depth"] + 1, "parent": tid, "index": k}
                            self.ready.append(cid)
                    self._report(tid, doc)
                elif not self.killing:
                    self.failed.append(tid)
                    log_line(self.log_path, f"stage A: task {tid} failed: {err}")
                self.cond.notify_all()

    def _report(self, tid, doc):
        tr = doc["trace"][0] if doc["trace"] else None
        what = (f"{tr['status']} lb {tr['lb']:.1f} ({tr['short']:.1f} short) {tr['iters']} it warm {tr['warm']}"
                + (f" cut {tr.get('cut_status')} best {tr.get('cut_best', 0):.1f}" if tr['cut'] else "")
                if tr else "no DNN")
        log_line(self.log_path, f"task {tid} d{doc['depth']} [{doc['spec']}] cpu {doc['stats']['times'][0]:.1f} s "
                                f"wall {doc['wall']:.0f} s: {what}; {doc['n_children']} children ({doc['line']})"
                                f"{' retained' if doc.get('retained') else ''}")
        p = self.projection()
        write_json(os.path.join(self.run, "status.json"), p)
        if p.get("total_cpu_est"):
            log_line(self.log_path, f"projection: {p['done']} tasks done, CPU {p['done_cpu']:.0f} s; expected nodes "
                                    f"{p['expected']}; total CPU {p['total_cpu_est']:.0f} s, remaining wall "
                                    f"{p['remaining_wall_h']:.2f} h on {len(self.cpus)} slots; disk free "
                                    f"{disk_free_gb(self.run):.1f} GB")

    def projection(self):
        """expected tasks per depth (children per done node at each depth) and the
        CPU they cost at each depth's measured mean"""
        D = self.cfg["depth"]
        by = {d: [x for x in self.done.values() if x["depth"] == d] for d in range(D + 1)}
        expected, total = {}, 0.0
        e = 1.0
        for d in range(D + 1):
            expected[d] = e
            done = by[d]
            if done:
                mean_cpu = sum(x["stats"]["times"][0] for x in done) / len(done)
                total += e * mean_cpu
                mean_kids = sum(x["n_children"] for x in done) / len(done)
                e = sum(x["n_children"] for x in done) + max(0.0, e - len(done)) * mean_kids
            else:
                total = None
                break
        done_cpu = sum(x["stats"]["times"][0] for x in self.done.values())
        out = {"updated": time.strftime("%F %T"), "done": len(self.done), "running": len(self.running),
               "ready": len(self.ready), "done_cpu": done_cpu,
               "expected": {d: round(v, 1) for d, v in expected.items()}, "frontier_expected": round(e, 1)}
        if total is not None:
            out["total_cpu_est"] = total
            out["remaining_wall_h"] = max(0.0, total - done_cpu) / len(self.cpus) / 3600
        return out

    # ------------------------------------------------------------ run
    def signal(self, *_):
        with self.cond:
            if not self.draining:
                self.draining = True
                log_line(self.log_path, "stage A: draining (running tasks finish; signal again to kill them)")
            else:
                self.killing = True
                log_line(self.log_path, "stage A: killing the running solvers; their tasks rerun next start")
                for p in list(self.procs.values()):
                    try:
                        os.killpg(p.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            self.cond.notify_all()

    def go(self):
        signal.signal(signal.SIGTERM, self.signal)
        signal.signal(signal.SIGINT, self.signal)
        log_line(self.log_path, f"stage A: {len(self.cpus)} slots on CPUs {self.cpus}, instance "
                                f"{self.cfg['instance']}, -x {self.cfg['xspec']}, -u {self.cfg['ub']}, sym "
                                f"{self.cfg['sym']}, overrides {self.cfg['env_overrides']}")
        threads = [threading.Thread(target=self._slot, args=(c,), daemon=True) for c in self.cpus]
        for th in threads:
            th.start()
        for th in threads:
            while th.is_alive():
                th.join(timeout=1.0)
        finished = not self.ready and not self.running and not self.failed and len(self.done) == len(self.tasks)
        if finished:
            self.assemble()
            log_line(self.log_path, f"stage A: finished, {len(self.done)} tasks")
            return 0
        log_line(self.log_path, f"stage A: stopped with {len(self.done)} of {len(self.tasks)} known tasks done, "
                                f"{len(self.ready)} ready, {len(self.failed)} failed; run again to resume")
        return 3

    def assemble(self):
        """frontier.state (every frontier node, slim warm starts, the stage's
        summed statistics) and frontier.tsv"""
        D = self.cfg["depth"]
        root = State(os.path.join(self.tdir("r"), "out.state")) if D > 0 else None
        totals = None
        nodes, rows = [], []
        for tid in sorted(self.done, key=lambda s: [int(x) if x.isdigit() else -1 for x in s.split(".")]):
            doc = self.done[tid]
            if doc["depth"] < D:
                continue
            st = State(os.path.join(self.tdir(tid), "out.slim.state"))
            nodes += st.nodes
            for k, f in enumerate(doc["frontier"]):
                rows.append(f"{tid}\t{k}\t{f['spec']}\t{f['bound']!r}\t{f['pred']!r}\t{f['hash']:.12f}")
        # the statistics: the sum over every task (each task's in.state had them zeroed)
        for tid, doc in self.done.items():
            p = os.path.join(self.tdir(tid), "out.slim.state" if doc["depth"] >= D else "out.state")
            st = State(p)
            if totals is None:
                totals = st
            else:
                totals.add_stats(st)
        hdr = root or totals
        totals.obj, totals.found, totals.perm, totals.root_bound = hdr.obj, hdr.found, hdr.perm, hdr.root_bound
        totals.write(os.path.join(self.run, "frontier.state"), nodes)
        with open(os.path.join(self.run, "frontier.tsv"), "w") as f:
            f.write("parent\tindex\tspec\tbound\tpredicted\thash\n" + "\n".join(rows) + "\n")
        log_line(self.log_path, f"stage A: frontier.state with {len(nodes)} nodes")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=os.path.join(D, "runs", "stageA"))
    ap.add_argument("--instance", help="default: fullrun/settings.sh's FULLRUN_INSTANCE")
    ap.add_argument("--sym", help="default: FULLRUN_SYM")
    ap.add_argument("--ub", type=float, help="default: FULLRUN_UB")
    ap.add_argument("--params", default=PARAMS)
    ap.add_argument("--xspec", default=XSPEC)
    ap.add_argument("--depth", type=int, default=2, help="the last depth bounded; its children are the frontier")
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--cpus", default="0,1,2,3", help="one CPU per slot (taskset)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--tau", type=float, default=TAU)
    ap.add_argument("--max-new", type=int, help="start at most this many tasks in this session, then drain")
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--env", action="append", default=[], help="K=V on top of the settings (validation runs only)")
    ap.add_argument("--binary", default=os.path.join(R, "latest", "bnbdrv"))
    a = ap.parse_args()
    for p in ("instance", "sym", "params"):
        v = getattr(a, p)
        if v:
            setattr(a, p, os.path.abspath(v))
    sys.exit(StageA(a).go())


if __name__ == "__main__":
    main()
