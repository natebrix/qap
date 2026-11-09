#!/usr/bin/env python3
"""Stage B of the wil50 sampling plan (README.md; roadmap/tree_size.md
section 3): take a seeded random sample of the stage-A frontier to the
give-up point the way the full run treats a depth-3 node's first solve, then
converge the farthest, and report the shortfall distribution and the class
fractions in the plan's strata.

  stage_b.py [--stage-a runs/stageA] [--run runs/stageB] [--n 400] [--converge 80]
             [--limit K] [--check-F K] [--calibrate SPEC] [--slots 4] [--partial]

The sample: the --n frontier nodes with the smallest unit_hash(seed, node)
(common.py), the seed being stage A's (recorded in its config.json), so the
sample is fixed before stage A runs. --limit K bounds only the first K of it
(the smoke test); the same run directory later reuses them.

Give-up step, one task per node: the node alone in a state file with its
parent's warm start exactly as the frontier carries it (the parent's QPB point,
its cut pool if any, and the parent's DNN state, kept by stage A for the
sampled nodes' parents in out.full.state), then
  bnbdrv -r node.state -d 3 -o out.state   (the settings, -x 1:4:1:300:1)
so the node gets its first DNN solve with the production give-up rule and
nothing else: the cut phase is off at depth 3 (QAP_DNN_CUT_MAXDEPTH=2; with
--phase it runs as in production). A parent stage A did not keep is rebuilt
by rerunning its stage-A task from its depth-1 parent (deterministic; checked
against stage A's record of its children).

Converge step: the --converge nodes farthest from the cutoff at the give-up
point, from the same state file with no give-up (-x 1:4:1:0:1) and no cut
phase: the DNN's own limit (4,000 iterations or its convergence test).

--check-F K: the first K nodes also run with no warm start at all and with
-F (which has none): the two first bounds must be equal.
--calibrate SPEC: the -F node SPEC (e.g. node 5 of tuning/giveup/all32.txt,
transposed: 49:23,15:37,5:6) through both steps, cold, against its recorded
values (329 short at 300 iterations, 315.7 converged after 1,500 iterations).

Resumable (a task with result.json is not rerun); SIGTERM drains.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import signal
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (D, State, disk_free_gb, fullrun_env, log_checks, log_line, node_fields,  # noqa: E402
                    parse_trace, run_bnbdrv, slim_copy, solver_env, spec, stats_of, write_json, write_one)

STRATA = [("fathomed", None, None), ("near (< 150)", 0, 150), ("150-210", 150, 210), ("210-300", 210, 300),
          ("tail (>= 300)", 300, math.inf)]


def wilson(k, n, z=1.0):
    """the 1-sigma Wilson interval of a fraction k/n"""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


class StageB:
    def __init__(self, a):
        self.a = a
        self.A = os.path.abspath(a.stage_a)
        self.run = os.path.abspath(a.run)
        os.makedirs(os.path.join(self.run, "tasks"), exist_ok=True)
        self.log_path = os.path.join(self.run, "launcher.log")
        self.acfg = json.load(open(os.path.join(self.A, "config.json")))
        self.binary = os.path.join(self.A, "bin", "bnbdrv")       # stage A's private copy: the same binary
        settings, _ = fullrun_env(self.binary)
        if settings != self.acfg["settings_env"]:
            raise SystemExit("fullrun/settings.sh changed since stage A ran; stage B must use stage A's settings")
        base = dict(self.acfg["env_overrides"])
        self.env_giveup = solver_env(settings, dict(base, **({} if a.phase else {"QAP_DNN_CUT_MAXDEPTH": "2"})))
        self.env_conv = solver_env(settings, dict(base, QAP_DNN_CUT_ZONE="0"))
        x = self.acfg["xspec"].split(":")
        self.x_giveup = self.acfg["xspec"]
        self.x_conv = ":".join(x[:3] + ["0"] + x[4:])
        self.cpus = [int(c) for c in a.cpus.split(",")][:a.slots]
        self.cpu_lock = threading.Lock()
        self.free_cpus = list(self.cpus)
        self.draining = False
        self.procs = {}
        self.depth = self.acfg["depth"] + 1
        self.sha = hashlib.sha256(open(self.binary, "rb").read()).hexdigest()
        self.t0 = time.time()
        cfg = {"stage_a": self.A, "n": a.n, "seed": self.acfg["seed"], "converge": a.converge, "phase": a.phase,
               "binary": self.binary, "binary_sha256": self.sha,
               "xspec_giveup": self.x_giveup, "xspec_converge": self.x_conv,
               "giveup_overrides": {} if a.phase else {"QAP_DNN_CUT_MAXDEPTH": "2"},
               "converge_overrides": {"QAP_DNN_CUT_ZONE": "0"}, "commit": self.acfg["commit"]}
        p = os.path.join(self.run, "config.json")
        if os.path.exists(p):
            old = json.load(open(p))
            for k in ("stage_a", "n", "seed", "phase", "xspec_giveup", "xspec_converge"):
                if old[k] != cfg[k]:
                    raise SystemExit(f"{p}: {k} was {old[k]!r}, now {cfg[k]!r}; use another --run")
            if old.get("binary_sha256") not in (None, self.sha):
                raise SystemExit(f"{p}: the binary's sha256 was {old['binary_sha256']}, now {self.sha}; use another --run")
            if old.get("binary_sha256") is None:
                # a run from before the sha was recorded (the smoke test): its tasks ran this same file
                cfg["binary_sha256_recorded_later"] = time.strftime("%F %T")
        write_json(p, cfg)

    # ------------------------------------------------------------ the sample
    def frontier(self):
        """every frontier node stage A recorded: (hash, parent task, index, spec, bound, pred)"""
        rows, complete = [], os.path.exists(os.path.join(self.A, "frontier.tsv"))
        tdir = os.path.join(self.A, "tasks")
        for t in os.listdir(tdir):
            p = os.path.join(tdir, t, "done.json")
            if not os.path.exists(p):
                continue
            doc = json.load(open(p))
            for k, f in enumerate(doc.get("frontier", [])):
                rows.append((f["hash"], t, k, f["spec"], f["bound"], f["pred"]))
        if not complete and not self.a.partial:
            raise SystemExit("stage A has not finished (no frontier.tsv); --partial samples what it has so far")
        rows.sort()
        return rows, complete

    def sample(self):
        rows, complete = self.frontier()
        s = rows[:self.a.n]
        old = os.path.join(self.run, "sample.tsv")
        if os.path.exists(old):
            prev = [l.rstrip("\n").split("\t") for l in open(old) if l[0].isdigit()]
            new = [[str(r), f"{x[0]:.12f}", x[1], str(x[2]), x[3], repr(x[4]), repr(x[5])] for r, x in enumerate(s)]
            if prev != new[:len(prev)]:
                raise SystemExit(f"{old}: the sample differs from this run's earlier one; use another --run")
        with open(os.path.join(self.run, "sample.tsv"), "w") as f:
            f.write(f"# {len(rows)} frontier nodes ({'complete' if complete else 'PARTIAL'}), seed "
                    f"{self.acfg['seed']}, the {len(s)} smallest hashes\n")
            f.write("rank\thash\tparent\tindex\tspec\tbound\tpredicted\n")
            for r, x in enumerate(s):
                f.write(f"{r}\t{x[0]:.12f}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]!r}\t{x[5]!r}\n")
        return s, len(rows), complete

    # ------------------------------------------------------------ parents
    def parent_state(self, parent):
        """the parent's full out.state (with its DNN state): stage A's, or rebuilt"""
        p = os.path.join(self.A, "tasks", parent, "out.full.state")
        if os.path.exists(p):
            return p
        q = os.path.join(self.run, "regen", parent, "out.state")
        if os.path.exists(q):
            return q
        raise RuntimeError(f"parent {parent} is not kept; rebuild it first")

    def regen(self, parent, cpu):
        """rerun stage A's task `parent` from its own parent's out.state (stage A keeps
        the depth-1 tasks' out.state); its children must match stage A's record"""
        doc = json.load(open(os.path.join(self.A, "tasks", parent, "done.json")))
        wd = os.path.join(self.run, "regen", parent)
        os.makedirs(wd, exist_ok=True)
        inp, out = os.path.join(wd, "in.state"), os.path.join(wd, "out.state.tmp")
        write_one(os.path.join(self.A, "tasks", doc["parent"], "out.state"), doc["index"], inp, warm="keep")
        c = self.acfg
        env = solver_env(c["settings_env"], c["env_overrides"])
        rc, wall = run_bnbdrv(self.binary, ["-p", c["params"], "-x", c["xspec"], "-s", c["sym"], "-u", repr(c["ub"]),
                                            "-d", str(doc["depth"]), "-r", inp, "-o", out, c["instance"]],
                              env, os.path.join(wd, "solver.log"), os.path.join(wd, "trace.txt"), cpu=cpu,
                              procs=self.procs, tag="regen " + parent)
        if rc != 0:
            raise RuntimeError(f"regen {parent}: exit {rc}")
        st = State(out)
        got = [(spec(node_fields(t, st.n)[3]), node_fields(t, st.n)[1]) for t, _ in st.nodes]
        want = [(f["spec"], f["bound"]) for f in doc["frontier"]]
        if got != want:
            raise RuntimeError(f"regen {parent}: children differ from stage A's record")
        os.replace(out, os.path.join(wd, "out.state"))
        os.remove(inp)
        log_line(self.log_path, f"rebuilt parent {parent} ({wall:.0f} s): children match stage A")

    # ------------------------------------------------------------ one node
    def node_task(self, key, src, index, mode, cpu, F=None):
        """mode: giveup | converge | cold (no warm start) | F (the -F node, cold);
        result.json with the trace row, the statistics and the children count"""
        wd = os.path.join(self.run, "tasks", key)
        res_p = os.path.join(wd, "result.json")
        if os.path.exists(res_p):
            return json.load(open(res_p))
        os.makedirs(wd, exist_ok=True)
        c = self.acfg
        conv = mode.endswith("converge")
        env, xs = (self.env_conv, self.x_conv) if conv else (self.env_giveup, self.x_giveup)
        out, inp = os.path.join(wd, "out.state"), os.path.join(wd, "in.state")
        args = ["-p", c["params"], "-x", xs, "-s", c["sym"], "-u", repr(c["ub"]), "-d", str(self.depth), "-o", out]
        if F:
            args += ["-F", F]
        else:
            write_one(src, index, inp, warm="none" if mode == "cold" else "keep")
            args += ["-r", inp]
        args.append(c["instance"])
        rc, wall = run_bnbdrv(self.binary, args, env, os.path.join(wd, "solver.log"), os.path.join(wd, "trace.txt"),
                              cpu=cpu, procs=self.procs, tag=key)
        if rc != 0 or not os.path.exists(out):
            raise RuntimeError(f"{key}: exit {rc}")
        st = State(out)
        tr = parse_trace(os.path.join(wd, "trace.txt"))
        res = {"key": key, "mode": mode, "cpu_slot": cpu, "wall": wall, "stats": stats_of(st),
               "trace": tr[0] if tr else None, "children": len(st.nodes), "checks": log_checks(os.path.join(wd, "solver.log")),
               "finished": time.strftime("%F %T")}
        slim_copy(out, os.path.join(wd, "out.slim.state"))
        os.remove(out)
        if os.path.exists(inp):
            os.remove(inp)
        write_json(res_p, res)
        return res

    # ------------------------------------------------------------ the pool
    def pool(self, jobs):
        """run jobs [(label, fn)] on the slots, each pinned to a free CPU; stops
        submitting on a drain"""
        results = {}

        def wrap(label, fn):
            with self.cpu_lock:
                cpu = self.free_cpus.pop(0)
            try:
                if not self.draining and disk_free_gb(self.run) < self.a.min_free_gb:
                    self.draining = True
                    log_line(self.log_path, f"stage B: only {disk_free_gb(self.run):.1f} GB free "
                                            f"(< {self.a.min_free_gb}); draining")
                if self.draining:
                    return label, None
                return label, fn(cpu)
            finally:
                with self.cpu_lock:
                    self.free_cpus.append(cpu)

        with cf.ThreadPoolExecutor(max_workers=len(self.cpus)) as ex:
            futs = [ex.submit(wrap, label, fn) for label, fn in jobs]
            for f in cf.as_completed(futs):
                try:
                    label, r = f.result()
                except Exception as e:      # noqa: BLE001
                    log_line(self.log_path, f"stage B: a task failed: {e}")
                    continue
                if r is not None:
                    results[label] = r
                    tr = r.get("trace") if isinstance(r, dict) else None
                    if tr:
                        log_line(self.log_path, f"{label}: {tr['status']} lb {tr['lb']:.2f} ({tr['short']:.1f} short) "
                                                f"{tr['iters']} it cpu {r['stats']['times'][0]:.0f} s warm {tr['warm']} "
                                                f"pool {tr['pool']}")
                    self.guard()
        return results

    def guard(self):
        """the stop rule: the projected wall of the whole stage B (this session's
        elapsed wall plus the remaining CPU on the slots) above --max-wall-h drains
        the run; logged after every task, in status.json too"""
        tdir = os.path.join(self.run, "tasks")
        g, v = [], []
        for t in os.listdir(tdir):
            r = self.load(t)
            if r and t.startswith("g_"):
                g.append(r["stats"]["times"][0])
            elif r and t.startswith("v_"):
                v.append(r["stats"]["times"][0])
        if not g:
            return
        n_conv = min(self.a.converge, self.a.n)
        mg = sum(g) / len(g)
        mv = sum(v) / len(v) if v else 1267.0          # the smoke test's converge until one is measured here
        rem = max(0, self.a.n - len(g)) * mg + max(0, n_conv - len(v)) * mv
        el = (time.time() - self.t0) / 3600
        wall = el + rem / len(self.cpus) / 0.99 / 3600
        doc = {"updated": time.strftime("%F %T"), "giveups_done": len(g), "giveup_mean_cpu": mg,
               "converges_done": len(v), "converge_mean_cpu": mv, "remaining_cpu": rem,
               "session_elapsed_h": el, "projected_wall_h": wall, "disk_free_gb": disk_free_gb(self.run)}
        write_json(os.path.join(self.run, "status.json"), doc)
        log_line(self.log_path, f"projection: {len(g)} give-ups done ({mg:.0f} s each), {len(v)} converges "
                                f"({mv:.0f} s each{'' if v else ', assumed'}); remaining {rem:,.0f} CPU s; wall "
                                f"{el:.2f} h so far + {wall - el:.2f} h = {wall:.2f} h; disk free "
                                f"{doc['disk_free_gb']:.1f} GB")
        if (self.a.max_wall_h and not self.draining and len(g) >= 40 and wall > self.a.max_wall_h):
            self.draining = True
            log_line(self.log_path, f"stage B: the projected wall {wall:.1f} h exceeds --max-wall-h "
                                    f"{self.a.max_wall_h}; draining")

    def signal(self, *_):
        if not self.draining:
            self.draining = True
            log_line(self.log_path, "stage B: draining (running tasks finish)")
        else:
            log_line(self.log_path, "stage B: killing the running solvers")
            for p in list(self.procs.values()):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    # ------------------------------------------------------------ run
    def go(self):
        signal.signal(signal.SIGTERM, self.signal)
        signal.signal(signal.SIGINT, self.signal)
        s, N, complete = self.sample()
        todo = s[:self.a.limit] if self.a.limit else s
        log_line(self.log_path, f"stage B: frontier {N} nodes ({'complete' if complete else 'partial'}), sample "
                                f"{len(s)}, bounding {len(todo)} on {len(self.cpus)} slots; disk free "
                                f"{disk_free_gb(self.run):.1f} GB")
        # parents to rebuild
        need = sorted({x[1] for x in todo if not os.path.exists(os.path.join(self.A, "tasks", x[1], "out.full.state"))
                       and not os.path.exists(os.path.join(self.run, "regen", x[1], "out.state"))})
        if need:
            log_line(self.log_path, f"stage B: rebuilding {len(need)} parents stage A did not keep")
            self.pool([(f"regen {p}", (lambda p: lambda cpu: self.regen(p, cpu) or {})(p)) for p in need])
        jobs = []
        for h, parent, idx, sp, b, pb in todo:
            key = "g_" + sp.replace(":", "_").replace(",", "-")
            jobs.append((sp, (lambda key, parent, idx: lambda cpu: self.node_task(
                key, self.parent_state(parent), idx, "giveup", cpu))(key, parent, idx)))
        for h, parent, idx, sp, b, pb in todo[:self.a.check_F]:
            k = sp.replace(":", "_").replace(",", "-")
            jobs.append((sp + " cold", (lambda k, parent, idx: lambda cpu: self.node_task(
                "c_" + k, self.parent_state(parent), idx, "cold", cpu))(k, parent, idx)))
            jobs.append((sp + " -F", (lambda k, sp: lambda cpu: self.node_task("F_" + k, None, None, "F", cpu, F=sp))(k, sp)))
        if self.a.calibrate:
            k = self.a.calibrate.replace(":", "_").replace(",", "-")
            for mode in ("Fgiveup", "Fconverge"):
                jobs.append((f"calibrate {self.a.calibrate} {mode}", (lambda k, mode: lambda cpu: self.node_task(
                    f"cal{mode}_{k}", None, None, mode, cpu, F=self.a.calibrate))(k, mode)))
        res = self.pool(jobs)
        if self.draining:
            return self.finish(todo, N, complete, partial_run=True)
        # converge the farthest
        give = [(res[x[3]]["trace"]["short"], x) for x in todo if x[3] in res and res[x[3]]["trace"]
                and not res[x[3]]["trace"]["fathomed"]]
        give.sort(key=lambda t: -t[0])
        far = [x for _, x in give[:self.a.converge]]
        log_line(self.log_path, f"stage B: converging the {len(far)} farthest (give-up shortfall "
                                f"{give[0][0] if give else 0:.1f} to {give[len(far) - 1][0] if far else 0:.1f})")
        jobs = [(x[3] + " converge", (lambda key, parent, idx: lambda cpu: self.node_task(
            key, self.parent_state(parent), idx, "converge", cpu))("v_" + x[3].replace(":", "_").replace(",", "-"),
                                                                   x[1], x[2])) for x in far]
        self.pool(jobs)
        return self.finish(todo, N, complete, partial_run=self.draining)

    # ------------------------------------------------------------ the projection
    def projection(self, rows):
        """the full stage B's CPU and wall [e]: n x the mean give-up CPU of the nodes
        bounded so far, plus the converge step at the mean CPU per converged node
        so far (or 4,000 iterations at the give-up runs' CPU per iteration), plus
        the rebuilt parents at their stage-A CPU (exact, from done.json)"""
        try:
            s, _, _ = self.sample()
        except SystemExit:
            return "Projection: no complete frontier yet."
        give = [r["giveup"]["stats"]["times"][0] for r in rows]
        if not give:
            return "Projection: no node bounded yet."
        g = sum(give) / len(give)
        conv = [r["converge"]["stats"]["times"][0] for r in rows if r["converge"]]
        per_it = [r["giveup"]["trace"]["cpu"] / r["giveup"]["trace"]["iters"] for r in rows
                  if r["giveup"]["trace"] and r["giveup"]["trace"]["iters"] >= 100]
        if conv:
            v, vsrc = sum(conv) / len(conv), f"measured on {len(conv)}"
        elif per_it:
            v, vsrc = 4000 * sum(per_it) / len(per_it), "4,000 iterations at the give-up runs' CPU per iteration"
        else:
            v, vsrc = 0.0, "unknown"
        parents = {x[1] for x in s if not os.path.exists(os.path.join(self.A, "tasks", x[1], "out.full.state"))}
        rb = sum(json.load(open(os.path.join(self.A, "tasks", p, "done.json")))["stats"]["times"][0] for p in parents)
        # the full design (--n, and --project-converge converged nodes), whatever this run's --limit/--converge
        nc = self.a.project_converge
        hop = [r["giveup"]["stats"]["times"][0] for r in rows if r["giveup"]["trace"]
               and r["giveup"]["trace"]["status"] == "hopeless"]
        gmax = max(g, sum(hop) / len(hop)) if hop else g
        total = self.a.n * g + nc * v + rb
        upper = self.a.n * gmax + nc * v + rb
        slots = len(self.cpus)
        return (f"Projection for the full stage B [e]: {self.a.n} give-ups x {g:.0f} s (mean of {len(give)}) = "
                f"{self.a.n * g:,.0f} s; {nc} converges x {v:.0f} s ({vsrc}) = {nc * v:,.0f} s; "
                f"{len(parents)} parents to rebuild, {rb:,.0f} s at their stage-A CPU; total {total:,.0f} CPU s, "
                f"{total / slots / 0.99 / 3600:.1f} h of wall on {slots} slots; if every give-up cost what the "
                f"hopeless ones did ({gmax:.0f} s): {upper:,.0f} CPU s, {upper / slots / 0.99 / 3600:.1f} h")

    # ------------------------------------------------------------ the report
    def load(self, key):
        p = os.path.join(self.run, "tasks", key, "result.json")
        return json.load(open(p)) if os.path.exists(p) else None

    def finish(self, todo, N, complete, partial_run=False):
        k = lambda sp: sp.replace(":", "_").replace(",", "-")
        rows = []
        for h, parent, idx, sp, b, pb in todo:
            g = self.load("g_" + k(sp))
            if g is None:
                continue
            rows.append({"spec": sp, "parent": parent, "inherited": b, "giveup": g, "converge": None})
        # the converged set is the --converge misses farthest at the give-up point (as go() chose them);
        # a converge result of any other node (an earlier run's) is listed apart, not used for classes
        miss = sorted((r for r in rows if r["giveup"]["trace"] and not r["giveup"]["trace"]["fathomed"]),
                      key=lambda r: -r["giveup"]["trace"]["short"])
        far = {r["spec"] for r in miss[:self.a.converge]}
        extra_conv = []
        for r in rows:
            v = self.load("v_" + k(r["spec"]))
            if v and r["spec"] in far:
                r["converge"] = v
            elif v:
                extra_conv.append((r["spec"], v))
        cutoff = next((r["giveup"]["trace"]["cutoff"] for r in rows if r["giveup"]["trace"]), None)
        lines = [f"# Stage B: {len(rows)} of {len(todo)} sampled nodes bounded (frontier {N} nodes, "
                 f"{'complete' if complete else 'PARTIAL'}; seed {self.acfg['seed']}; sample of {self.a.n})",
                 f"# binary {self.acfg['commit']['commit'][:10]}, give-up -x {self.x_giveup}, converge -x {self.x_conv}"
                 f", phase {'on' if self.a.phase else 'off'} at depth {self.depth}", ""]
        # per node
        lines.append("| node | parent | inherited | first lb | short | status | iters | CPU s | warm | pool | "
                     "children | converged lb | converged short | conv iters | conv CPU s |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        cls = {name: 0 for name, *_ in STRATA}
        cls_by_conv = 0
        cpu_g = cpu_v = 0.0
        shorts = []
        for r in rows:
            t = r["giveup"]["trace"]
            v = r["converge"]["trace"] if r["converge"] else None
            cpu_g += r["giveup"]["stats"]["times"][0]
            if r["converge"]:
                cpu_v += r["converge"]["stats"]["times"][0]
            if t is None:
                lines.append(f"| {r['spec']} | {r['parent']} | {r['inherited']:.1f} | (QPB fathomed) | | | | | | | 0 | | | | |")
                cls["fathomed"] += 1
                continue
            lines.append(f"| {r['spec']} | {r['parent']} | {r['inherited']:.1f} | {t['lb']:.2f} | {t['short']:.1f} | "
                         f"{t['status']} | {t['iters']} | {r['giveup']['stats']['times'][0]:.0f} | {t['warm']} | "
                         f"{t['pool']} | {r['giveup']['children']} | "
                         + (f"{v['lb']:.2f} | {v['short']:.1f} | {v['iters']} | {r['converge']['stats']['times'][0]:.0f} |"
                            if v else " | | | |"))
            if t["fathomed"]:
                cls["fathomed"] += 1
                continue
            s = v["short"] if v else t["short"]
            if v and v["fathomed"]:
                s = -1.0
            cls_by_conv += bool(v)
            shorts.append(t["short"])
            for name, lo, hi in STRATA[1:]:
                if (lo == 0 and s < hi) or (lo is not None and lo <= s < hi):
                    cls[name] += 1
                    break
        n = len(rows)
        lines += ["", f"Cutoff {cutoff}. Classes (converged shortfall for the {cls_by_conv} converged nodes, the "
                      f"give-up shortfall, an upper bound on it, for the others):", "",
                  "| class | nodes | fraction | 1-sigma interval |", "|---|---|---|---|"]
        for name, *_ in STRATA:
            lo, hi = wilson(cls[name], n)
            lines.append(f"| {name} | {cls[name]} | {cls[name] / n if n else 0:.3f} | {lo:.3f}-{hi:.3f} |")
        if shorts:
            q = sorted(shorts)
            pct = lambda p: q[min(len(q) - 1, int(p * len(q)))]
            lines += ["", f"Give-up shortfall of the {len(q)} misses: min {q[0]:.1f}, 10% {pct(.1):.1f}, median "
                          f"{pct(.5):.1f}, 90% {pct(.9):.1f}, max {q[-1]:.1f}"]
        # the give-up shortfall histogram (every sampled node)
        edges = [0, 50, 100, 150, 210, 300, 350, 400, 450, 500, 1e9]
        names = ["0-50", "50-100", "100-150", "150-210", "210-300", "300-350", "350-400", "400-450", "450-500", ">= 500"]
        fat = sum(1 for r in rows if r["giveup"]["trace"] is None or r["giveup"]["trace"]["fathomed"])
        cnt = [sum(1 for x in shorts if edges[i] <= x < edges[i + 1]) for i in range(len(names))]
        lines += ["", "Give-up shortfall histogram (the first solve, warm, at most 300 iterations when hopeless):", "",
                  "| fathomed | " + " | ".join(names) + " |", "|---|" + "---|" * len(names),
                  f"| {fat} | " + " | ".join(str(x) for x in cnt) + " |"]
        # converged against give-up
        cv = [(r["giveup"]["trace"]["short"], r["converge"]["trace"]["short"], r) for r in rows
              if r["converge"] and r["converge"]["trace"] and r["giveup"]["trace"]]
        if cv:
            gains = sorted(a - b for a, b, _ in cv)
            pg = lambda p: gains[min(len(gains) - 1, int(p * len(gains)))]
            def klass(x):
                for name, lo, hi in STRATA[1:]:
                    if (lo == 0 and x < hi) or (lo is not None and lo <= x < hi):
                        return name
            moved = sum(1 for a, b, _ in cv if klass(a) != klass(b))
            lines += ["", f"Converged against give-up on the {len(cv)} converged nodes (give-up shortfall "
                          f"{min(a for a, _, _ in cv):.1f}-{max(a for a, _, _ in cv):.1f}): the converge gains "
                          f"min {gains[0]:.1f}, 10% {pg(.1):.1f}, median {pg(.5):.1f}, 90% {pg(.9):.1f}, max "
                          f"{gains[-1]:.1f} units; {moved} change class; "
                          f"{sum(1 for _, _, r in cv if r['converge']['trace']['fathomed'])} fathom"]
        # CPU per step
        gc = sorted(r["giveup"]["stats"]["times"][0] for r in rows)
        vc = sorted(r["converge"]["stats"]["times"][0] for r in rows if r["converge"])
        med = lambda xs: xs[len(xs) // 2] if xs else 0.0
        hop = [r["giveup"]["stats"]["times"][0] for r in rows if r["giveup"]["trace"]
               and r["giveup"]["trace"]["status"] == "hopeless"]
        fth = [r["giveup"]["stats"]["times"][0] for r in rows if r["giveup"]["trace"]
               and r["giveup"]["trace"]["fathomed"]]
        lines += ["", f"CPU: give-up {cpu_g:,.0f} s on {n} nodes (mean {cpu_g / n if n else 0:.0f}, median "
                      f"{med(gc):.0f}; fathomed {sum(fth) / len(fth) if fth else 0:.0f} s, hopeless "
                      f"{sum(hop) / len(hop) if hop else 0:.0f} s); converge {cpu_v:,.0f} s on {len(vc)} nodes (mean "
                      f"{cpu_v / len(vc) if vc else 0:.0f}, median {med(vc):.0f})"]
        for sp, v in extra_conv:
            if v["trace"]:
                lines.append(f"Converge outside the farthest {self.a.converge} (an earlier run's; not used for the "
                             f"classes): {sp} {v['trace']['short']:.1f} short after {v['trace']['iters']} iterations")
        extra = [("cold / -F", self.load("c_" + k(sp)), self.load("F_" + k(sp)), sp) for _, _, _, sp, _, _ in
                 todo[:self.a.check_F]]
        for name, c, F, sp in extra:
            if c and F:
                ct, Ft = c["trace"], F["trace"]
                same = (ct is None and Ft is None) or (ct and Ft and ct["lb"] == Ft["lb"] and ct["iters"] == Ft["iters"])
                lines.append(f"Restart check {sp}: no warm start {ct['lb'] if ct else 'QPB'} ({ct['iters'] if ct else 0} "
                             f"it), -F {Ft['lb'] if Ft else 'QPB'} ({Ft['iters'] if Ft else 0} it): "
                             f"{'same' if same else 'DIFFERENT'}")
        if self.a.calibrate:
            kk = k(self.a.calibrate)
            for mode in ("Fgiveup", "Fconverge"):
                r = self.load(f"cal{mode}_{kk}")
                if r and r["trace"]:
                    t = r["trace"]
                    lines.append(f"Calibration {self.a.calibrate} {mode}: lb {t['lb']:.2f}, {t['short']:.1f} short, "
                                 f"{t['iters']} iterations, {t['status']}, CPU {r['stats']['times'][0]:.0f} s")
        lines += ["", self.projection(rows)]
        if partial_run:
            lines.append("(the run was stopped before the end: run the same command again to resume)")
        txt = "\n".join(lines) + "\n"
        with open(os.path.join(self.run, "summary.md"), "w") as f:
            f.write(txt)
        print(txt)
        return 3 if partial_run else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage-a", default=os.path.join(D, "runs", "stageA"))
    ap.add_argument("--run", default=os.path.join(D, "runs", "stageB"))
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--converge", type=int, default=80)
    ap.add_argument("--project-converge", type=int, default=80, help="converged nodes in the projected full stage B")
    ap.add_argument("--limit", type=int, help="bound only the first K sampled nodes (the smoke test)")
    ap.add_argument("--check-F", type=int, default=0)
    ap.add_argument("--calibrate", help="a -F spec run cold through both steps (node 5: 49:23,15:37,5:6)")
    ap.add_argument("--phase", action="store_true", help="let the cut phase run at depth 3 as in production")
    ap.add_argument("--max-wall-h", type=float, help="drain when the projected wall of the stage exceeds this")
    ap.add_argument("--min-free-gb", type=float, default=2.0, help="drain when the disk has less free")
    ap.add_argument("--partial", action="store_true", help="sample an unfinished stage A (smoke tests only)")
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--cpus", default="0,1,2,3")
    a = ap.parse_args()
    sys.exit(StageB(a).go())


if __name__ == "__main__":
    main()
