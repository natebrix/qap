#!/usr/bin/env python3
"""Branching lines at depth-1 nodes (results/depth1_lines.md): at a depth-1
node of the production tree, does a line other than rule 9's give clearly
better depth-2 children?

  depth1_lines.py prep [--run DIR] [--nodes 1,5] [--binary B] [--diag-binary B2]
  depth1_lines.py kids [--run DIR] [--nlines 4] [--slots 4] [--cpus 0,1,2,3] [--max-wall-h 5] [--stop-after-min M]
  depth1_lines.py status [--run DIR]

prep, on this checkout's binary (a private copy) and the full run's settings
(fullrun/settings.sh, fullrun.par; fullrun_env_head below):
  the root, `bnbdrv -d 0 -o root/out.state` (QPB only, its line a row in the 15
  J1 cells); then each depth-1 node r.K (the K-th root child, as stage A numbers
  them) the way stage A runs a task: write_one(root/out.state, K) and
  `bnbdrv -r in.state -d 1 -o out.state`, so its children carry its warm start
  (the root's QPB point and its own fp32 DNN state). Each depth-1 node runs
  twice, on the binary and on the diagnostic binary (the same sources plus a
  stderr print under QAP_XBAR_DEBUG=2 of every free row's and column's rule-8
  score, the sum of squared DNN marginals); the two traces and child sets must
  be equal. Writes prep.json: the root bound and line, each node's trace, rule
  9's line, and the ranking of rows and columns together by rule 8's score.

kids: the candidate lines at each node are rule 9's line and the next --nlines-1
lines of that ranking (rule 9's own line left out). Every line's 49 depth-2
children are built from the node's out.state by the same rewrite (rewrite_line):
each child keeps the node's warm start (its QPB point and DNN state, from which
a child's start is cut at its own assignment, warm_child_offsets), depth 2, the
node's bound inherited; its assignment becomes the line's, and its predicted
bound the inherited bound (the state does not hold U for other lines). Rule 9's
line is rewritten too ("rule9"), so the comparison is paired; its children
straight from the solver's output, with production's predicted bounds, run as
"prod", the check of the rewrite. Each child runs to its first DNN solve only:
  bnbdrv --allow-settings-change -p child.par -r in.state -d 2 -o out.state
with child.par = fullrun.par + cut.maxdepth = 1 (no cut phase at depth 2; the
only setting that differs from the node's state record, hence the flag), the
production give-up rule (dnn.giveup = 300 when hopeless), and -d 2, so a miss
branches (rule 9 at depth 2, recorded) but nothing deeper is bounded. One task
per child, pinned to one CPU, in a seeded random order; resumable (a task with
result.json is not rerun); SIGTERM drains. The big states are deleted after
each task; result.json, trace.txt and solver.log stay.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (D, PARAMS, R, SETTINGS, SOLVER_PREFIXES, THREAD_VARS, State, disk_free_gb, log_checks,  # noqa: E402
                    log_line, node_fields, params_with, parse_trace, run_bnbdrv, solver_env, spec, stats_of,
                    write_json)

SEED = 20261109
CUTOFF_DELTA = 1.95          # the solver's fathoming tolerance on even objectives (cutoff = ub - 1.95)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fullrun_env_head(binary):
    """common.fullrun_env without its fullrun_check, whose golden dump is stale at
    fe31484: 1d40615 added the run setting QAP_DNN_WARM_DIR (not in the settings
    hash), so the dump has 65 entries against the golden's 64. Accepted only when
    that is the whole difference; the diff is returned for the record."""
    base = {k: v for k, v in os.environ.items() if not k.startswith(SOLVER_PREFIXES)}
    script = ('. "$1" || exit 1; diff <(fullrun_dump "$2") "$FULLRUN_ROOT/experiments/wil50_feasibility/fullrun/settings.golden"'
              ' > "$3"; printf "FULLRUN_INSTANCE=%s\\0FULLRUN_SYM=%s\\0FULLRUN_UB=%s\\0FULLRUN_PAR=%s\\0--\\0" '
              '"$FULLRUN_INSTANCE" "$FULLRUN_SYM" "$FULLRUN_UB" "$FULLRUN_PAR"; env -0')
    dpath = binary + ".golden_diff"
    r = subprocess.run(["bash", "-c", script, "bash", SETTINGS, binary, dpath], env=base, capture_output=True)
    if r.returncode != 0:
        raise SystemExit("settings.sh failed:\n" + r.stderr.decode())
    diff = open(dpath).read()
    body = [l for l in diff.splitlines() if l and l[0] in "<>"]
    ok = all(l.startswith(("< settings 65 hash ", "> settings 64 hash ", '< setting QAP_DNN_WARM_DIR "" default run'))
             for l in body)
    hashes = {l.split()[-1] for l in body if " hash " in l}
    if not ok or len(hashes) > 1:
        raise SystemExit(f"{binary}: the settings under fullrun.par differ from settings.golden beyond the "
                         f"QAP_DNN_WARM_DIR entry:\n{diff}")
    head, _, rest = r.stdout.decode().partition("--\0")
    info = dict(x.split("=", 1) for x in head.split("\0") if "=" in x)
    env = {}
    for x in rest.split("\0"):
        if "=" not in x:
            continue
        k, v = x.split("=", 1)
        if k.startswith(SOLVER_PREFIXES) or k in THREAD_VARS:
            env[k] = v
    return env, info, diff


def child_line(own, kids):
    """'row i' / 'col j' from the children's assignments (stage_a.child_line)"""
    new = [sorted(set(k) - set(own)) for k in kids]
    new = [x[0] for x in new if len(x) == 1]
    if len(new) < 2:
        return "?" if new else "-"
    if len({i for i, _ in new}) == 1:
        return f"row {new[0][0]}"
    if len({j for _, j in new}) == 1:
        return f"col {new[0][1]}"
    return "mixed"


def read_state_any(path):
    """common.write_one accepts state versions 6-8; HEAD writes 9 (the @path DNN
    references, unused here: no QAP_DNN_WARM_DIR, so every token is inline)"""
    st = State(path)
    if st.magic[1] not in ("6", "7", "8", "9"):
        raise ValueError(f"{path}: state version {st.magic[1]}")
    return st


def write_node(st, node, out_path):
    """a state file holding `node` (text, warm) alone, with st's header and zeroed statistics"""
    st.zero_stats()
    st.write(out_path, [node])


def rewrite_line(st, line):
    """the depth-2 children of a depth-1 node on `line` ('row f' or 'col c'), from
    the node's out.state `st` (the solver's children on its own line): every
    field of the solver's first child is kept (depth, inherited bound, the warm
    start) except the assignment, its history (history[1] = the facility fixed
    at depth 2) and the predicted bound, which becomes the inherited bound.
    Children in increasing order of the free location (row) or facility (col)."""
    n = st.n
    text0, w = st.nodes[0]
    tok = text0.split()
    if int(tok[1]) != 2 or int(tok[4]) != 2:
        raise RuntimeError(f"the node's children are not depth-2 nodes with two assignments: {text0[:80]}")
    if any(x[1] is not w for x in st.nodes):
        raise RuntimeError("the node's children do not share one warm start")
    ps = [[int(x) for x in t.split()[5:5 + n]] for t, _ in st.nodes]
    hs = [t.split()[5 + n:5 + 2 * n] for t, _ in st.nodes]
    # the node's own assignment: what every child shares, and history[0]
    par = [col[0] if len(set(col)) == 1 else -1 for col in zip(*ps)]
    fixed = [(i, j) for i, j in enumerate(par) if j >= 0]
    if len(fixed) != 1 or any(h[0] != str(fixed[0][0]) for h in hs):
        raise RuntimeError(f"the node's assignment is not one fix: {fixed}")
    f0, c0 = fixed[0]
    kind, idx = line.split()
    idx = int(idx)
    pairs = []
    if kind == "row":
        if par[idx] >= 0:
            raise RuntimeError(f"row {idx} is fixed at the node")
        pairs = [(idx, j) for j in range(n) if j not in par]
    else:
        if idx in par:
            raise RuntimeError(f"col {idx} is fixed at the node")
        pairs = [(i, idx) for i in range(n) if par[i] < 0]
    nodes = []
    for i, j in pairs:
        p = list(par)
        p[i] = j
        hist = [str(f0), str(i)] + ["-1"] * (n - 2)
        t = list(tok)
        t[3] = t[2]
        nodes.append((" ".join(t[:5] + [str(x) for x in p] + hist), w))
    return nodes, (f0, c0)


class Run:
    def __init__(self, a):
        self.a = a
        self.run = os.path.abspath(a.run)
        os.makedirs(self.run, exist_ok=True)
        self.log_path = os.path.join(self.run, "launcher.log")
        cfg_p = os.path.join(self.run, "config.json")
        if os.path.exists(cfg_p):
            self.cfg = json.load(open(cfg_p))
        else:
            self.cfg = self.make_config()
            write_json(cfg_p, self.cfg)
        c = self.cfg
        self.env = solver_env(c["settings_env"])
        self.cpus = [int(x) for x in a.cpus.split(",")][:a.slots]
        self.procs = {}
        self.draining = False

    def make_config(self):
        a = self.a
        bins = {}
        for tag, src in (("bin", a.binary), ("bin_diag", a.diag_binary)):
            d = os.path.join(self.run, tag)
            os.makedirs(d, exist_ok=True)
            dst = os.path.join(d, "bnbdrv")
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
            bins[tag] = {"path": dst, "sha256": sha256(dst), "source": os.path.abspath(src)}
        settings, info, gdiff = fullrun_env_head(bins["bin"]["path"])
        settings_d, _, _ = fullrun_env_head(bins["bin_diag"]["path"])
        if settings != settings_d:
            raise SystemExit("the two binaries resolve different settings")
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=R, capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "latest"], cwd=R).returncode != 0
        params = PARAMS            # the root and the depth-1 nodes: fullrun.par itself
        child_par = os.path.join(self.run, "child.par")
        params_with(PARAMS, {"cut.maxdepth": "1"}, child_par, bins["bin"]["path"])
        return {"instance": info["FULLRUN_INSTANCE"], "sym": info["FULLRUN_SYM"], "ub": float(info["FULLRUN_UB"]),
                "params": params, "child_params": child_par, "params_base": PARAMS, "settings_env": settings,
                "binaries": bins, "golden_diff": gdiff, "commit": commit, "latest_modified": dirty, "nodes": a.nodes, "seed": SEED,
                "created": time.strftime("%F %T")}

    # ------------------------------------------------------------ one solver task
    def task(self, wd, args, cpu, binary="bin", extra_env=None, tag=None):
        os.makedirs(wd, exist_ok=True)
        c = self.cfg
        env = dict(self.env)
        env.update(extra_env or {})
        full = args + ["-s", c["sym"], "-u", repr(c["ub"]), c["instance"]]
        rc, wall = run_bnbdrv(c["binaries"][binary]["path"], full, env, os.path.join(wd, "solver.log"),
                              os.path.join(wd, "trace.txt"), cpu=cpu, procs=self.procs, tag=tag or wd)
        return rc, wall

    # ------------------------------------------------------------ prep
    def prep(self):
        c = self.cfg
        rootd = os.path.join(self.run, "root")
        root_out = os.path.join(rootd, "out.state")
        if not os.path.exists(os.path.join(rootd, "done.json")):
            rc, wall = self.task(rootd, ["-p", c["params"], "-d", "0", "-o", root_out], self.cpus[0])
            if rc != 0:
                raise SystemExit(f"root: exit {rc}")
            st = read_state_any(root_out)
            kids = [node_fields(t, st.n) for t, _ in st.nodes]
            doc = {"root_bound": float(st.root_bound), "children": [{"spec": spec(k[3]), "bound": k[1], "pred": k[2]}
                                                                     for k in kids],
                   "line": child_line([], [k[3] for k in kids]), "stats": stats_of(st), "wall": wall,
                   "checks": log_checks(os.path.join(rootd, "solver.log"))}
            write_json(os.path.join(rootd, "done.json"), doc)
        root = json.load(open(os.path.join(rootd, "done.json")))
        log_line(self.log_path, f"root: bound {root['root_bound']:.3f}, line {root['line']}, "
                                f"{len(root['children'])} children")
        st = read_state_any(root_out)
        jobs = []
        for K in c["nodes"]:
            for b in ("bin", "bin_diag"):
                wd = os.path.join(self.run, "d1", f"r.{K}", b)
                if os.path.exists(os.path.join(wd, "done.json")):
                    continue
                jobs.append((K, b, wd))
        lock = threading.Lock()
        free = list(self.cpus)

        def one(K, b, wd):
            with lock:
                cpu = free.pop(0)
            try:
                os.makedirs(wd, exist_ok=True)
                inp, out = os.path.join(wd, "in.state"), os.path.join(wd, "out.state")
                s2 = read_state_any(root_out)
                write_node(s2, s2.nodes[K], inp)
                rc, wall = self.task(wd, ["-p", c["params"], "-d", "1", "-r", inp, "-o", out], cpu, binary=b,
                                     extra_env={"QAP_XBAR_DEBUG": "2"} if b == "bin_diag" else None)
                if rc != 0:
                    raise RuntimeError(f"r.{K} {b}: exit {rc}")
                so = read_state_any(out)
                kids = [node_fields(t, so.n) for t, _ in so.nodes]
                fixes = node_fields(st.nodes[K][0], st.n)[3]
                logtxt = open(os.path.join(wd, "solver.log")).read()
                rank = [(m.group(1), int(m.group(2)), float(m.group(3))) for m in
                        re.finditer(r"^xbar9rank: depth 1 (row|col) (\d+) sumsq (\S+)$", logtxt, re.M)]
                x9 = re.search(r"^xbar9: depth 1 (row|col) (\d+) \(row (\S+) col (\S+)\)$", logtxt, re.M)
                doc = {"K": K, "binary": b, "spec": spec(fixes), "wall": wall, "stats": stats_of(so),
                       "trace": parse_trace(os.path.join(wd, "trace.txt")),
                       "trace_raw": open(os.path.join(wd, "trace.txt")).read(),
                       "line": child_line(fixes, [k[3] for k in kids]),
                       "children": [{"spec": spec(k[3]), "bound": k[1], "pred": k[2]} for k in kids],
                       "warm_sha256": hashlib.sha256(so.nodes[0][1].text.encode()).hexdigest() if so.nodes else None,
                       "rank": rank, "xbar9": x9.groups() if x9 else None,
                       "checks": log_checks(os.path.join(wd, "solver.log"))}
                os.remove(inp)
                if b == "bin_diag":
                    os.remove(out)          # the children are built from the binary's own out.state
                write_json(os.path.join(wd, "done.json"), doc)
                log_line(self.log_path, f"r.{K} [{doc['spec']}] {b}: {doc['trace'][0]['status']} lb "
                                        f"{doc['trace'][0]['lb']:.4f} {doc['trace'][0]['iters']} it, line {doc['line']}, "
                                        f"cpu {doc['stats']['times'][0]:.1f} s")
            finally:
                with lock:
                    free.append(cpu)

        if jobs:
            with cf.ThreadPoolExecutor(max_workers=len(self.cpus)) as ex:
                for f in [ex.submit(one, *j) for j in jobs]:
                    f.result()
        # checks and the candidate lines
        out = {"root": {k: root[k] for k in ("root_bound", "line")}, "nodes": {}}
        for K in c["nodes"]:
            d = {b: json.load(open(os.path.join(self.run, "d1", f"r.{K}", b, "done.json"))) for b in ("bin", "bin_diag")}
            strip = lambda s: re.sub(r"pid \d+ ", "", re.sub(r" cpu \S+", "", s))
            same_trace = strip(d["bin"]["trace_raw"]) == strip(d["bin_diag"]["trace_raw"])
            same_kids = d["bin"]["children"] == d["bin_diag"]["children"]
            same_warm = d["bin"]["warm_sha256"] == d["bin_diag"]["warm_sha256"]
            rank = sorted(d["bin_diag"]["rank"], key=lambda x: -x[2])
            rule = d["bin"]["line"]
            others = [f"{k} {i}" for k, i, _ in rank if f"{k} {i}" != rule]
            out["nodes"][f"r.{K}"] = {
                "spec": d["bin"]["spec"], "trace": d["bin"]["trace"][0], "rule9_line": rule,
                "rule9_lookahead": d["bin_diag"]["xbar9"], "diag_same_trace": same_trace,
                "diag_same_children": same_kids, "diag_same_warm": same_warm,
                "rank": [{"line": f"{k} {i}", "sumsq": s} for k, i, s in rank],
                "rule9_rank": next((r for r, (k, i, _) in enumerate(rank) if f"{k} {i}" == rule), None),
                "candidates": [rule] + others[:self.a.nlines - 1],
                "cpu": {b: d[b]["stats"]["times"][0] for b in d}}
            log_line(self.log_path, f"r.{K} [{d['bin']['spec']}]: rule 9 {rule} (rank {out['nodes'][f'r.{K}']['rule9_rank']} "
                                    f"by rule 8's score), diagnostic build same trace {same_trace}, children {same_kids}, "
                                    f"warm start {same_warm}; candidates {out['nodes'][f'r.{K}']['candidates']}")
        write_json(os.path.join(self.run, "prep.json"), out)
        return out

    # ------------------------------------------------------------ the children
    def kid_jobs(self):
        prep = json.load(open(os.path.join(self.run, "prep.json")))
        jobs = []
        for node, nd in prep["nodes"].items():
            src = os.path.join(self.run, "d1", node, "bin", "out.state")
            st = read_state_any(src)
            sets = {"prod": [(t, w) for t, w in st.nodes]}
            rule_nodes, _ = rewrite_line(st, nd["rule9_line"])
            got = sorted(spec(node_fields(t, st.n)[3]) for t, _ in rule_nodes)
            want = sorted(spec(node_fields(t, st.n)[3]) for t, _ in st.nodes)
            if got != want:
                raise SystemExit(f"{node}: the rewrite of rule 9's line does not give the solver's children")
            for line in nd["candidates"]:
                tag = "rule9" if line == nd["rule9_line"] else line.replace(" ", "")
                sets[tag] = rewrite_line(st, line)[0]
            for tag, nodes in sets.items():
                for k, nodeline in enumerate(nodes):
                    key = f"{node}/{tag}/{k:02d}"
                    jobs.append((key, node, tag, k, nodeline))
        jobs.sort(key=lambda j: hashlib.sha256(f"{SEED}:{j[0]}".encode()).hexdigest())
        return jobs, prep

    def kid_task(self, key, node, tag, k, nodeline, cpu):
        wd = os.path.join(self.run, "kids", key)
        res_p = os.path.join(wd, "result.json")
        if os.path.exists(res_p):
            return json.load(open(res_p))
        os.makedirs(wd, exist_ok=True)
        c = self.cfg
        src = os.path.join(self.run, "d1", node, "bin", "out.state")
        st = read_state_any(src)
        text, w = nodeline
        # the warm start object must be this file's (rewrite_line's came from another parse)
        w = st.nodes[0][1]
        inp, out = os.path.join(wd, "in.state"), os.path.join(wd, "out.state")
        write_node(st, (text, w), inp)
        rc, wall = self.task(wd, ["--allow-settings-change", "-p", c["child_params"], "-d", "2", "-r", inp, "-o", out],
                             cpu, tag=key)
        if rc != 0 or not os.path.exists(out):
            raise RuntimeError(f"{key}: exit {rc}")
        so = read_state_any(out)
        fixes = node_fields(text, st.n)[3]
        kids = [node_fields(t, so.n)[3] for t, _ in so.nodes]
        tr = parse_trace(os.path.join(wd, "trace.txt"))
        logtxt = open(os.path.join(wd, "solver.log"), errors="replace").read()
        warn = next((l for l in logtxt.splitlines() if "warning:" in l and "settings" in l), None)
        _, b, pb, _ = node_fields(text, st.n)
        res = {"key": key, "node": node, "line": tag, "k": k, "spec": spec(fixes), "inherited": b, "predicted": pb,
               "cpu_slot": cpu, "wall": wall, "stats": stats_of(so), "trace": tr[0] if tr else None,
               "children": len(kids), "child_line": child_line(fixes, kids) if kids else None,
               "checks": log_checks(os.path.join(wd, "solver.log")), "settings_warning": warn,
               "finished": time.strftime("%F %T")}
        for p in (inp, out):
            if os.path.exists(p):
                os.remove(p)
        write_json(res_p, res)
        return res

    def kids(self):
        jobs, prep = self.kid_jobs()
        todo = [j for j in jobs if not os.path.exists(os.path.join(self.run, "kids", j[0], "result.json"))]
        log_line(self.log_path, f"kids: {len(jobs)} tasks, {len(todo)} to run on {len(self.cpus)} slots; lines "
                                + "; ".join(f"{n}: {', '.join(v['candidates'])} + prod" for n, v in prep["nodes"].items()))
        signal.signal(signal.SIGTERM, self.signal)
        signal.signal(signal.SIGINT, self.signal)
        lock = threading.Lock()
        free = list(self.cpus)
        t0 = time.time()
        done_cpu = []
        n_total = len(jobs)
        n_done0 = len(jobs) - len(todo)

        def one(j):
            with lock:
                cpu = free.pop(0)
            try:
                if self.draining:
                    return None
                if self.a.stop_after_min and time.time() - t0 > 60 * self.a.stop_after_min:
                    return None             # a session slice (foreground calls): resume with the same command
                if disk_free_gb(self.run) < self.a.min_free_gb:
                    self.draining = True
                    log_line(self.log_path, f"kids: only {disk_free_gb(self.run):.1f} GB free; draining")
                    return None
                return self.kid_task(*j, cpu)
            finally:
                with lock:
                    free.append(cpu)

        with cf.ThreadPoolExecutor(max_workers=len(self.cpus)) as ex:
            futs = {ex.submit(one, j): j for j in todo}
            for f in cf.as_completed(futs):
                j = futs[f]
                try:
                    r = f.result()
                except Exception as e:      # noqa: BLE001
                    log_line(self.log_path, f"kids: {j[0]} failed: {e}")
                    continue
                if r is None:
                    continue
                done_cpu.append(r["stats"]["times"][0])
                tr = r["trace"]
                what = (f"{tr['status']} lb {tr['lb']:.2f} ({tr['short']:.1f} short) {tr['iters']} it"
                        if tr else "no DNN (QPB fathomed)")
                el = (time.time() - t0) / 3600
                rem = (len(todo) - len(done_cpu)) * (sum(done_cpu) / len(done_cpu))
                proj = el + rem / len(self.cpus) / 3600
                log_line(self.log_path, f"{r['key']} [{r['spec']}] cpu {r['stats']['times'][0]:.0f} s: {what}; "
                                        f"{r['children']} children; done {n_done0 + len(done_cpu)}/{n_total}, wall "
                                        f"{el:.2f} h + {proj - el:.2f} h; disk {disk_free_gb(self.run):.1f} GB")
                if self.a.max_wall_h and len(done_cpu) >= 20 and proj > self.a.max_wall_h and not self.draining:
                    self.draining = True
                    log_line(self.log_path, f"kids: projected {proj:.1f} h exceeds --max-wall-h {self.a.max_wall_h}; "
                                            "draining")
        left = sum(1 for j in jobs if not os.path.exists(os.path.join(self.run, "kids", j[0], "result.json")))
        log_line(self.log_path, f"kids: stopped, {len(jobs) - left} of {len(jobs)} done")
        return 0 if left == 0 else 3

    def signal(self, *_):
        if not self.draining:
            self.draining = True
            log_line(self.log_path, "kids: draining (running tasks finish; signal again to kill them)")
        else:
            log_line(self.log_path, "kids: killing the running solvers")
            for p in list(self.procs.values()):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["prep", "kids", "status"])
    ap.add_argument("--run", default=os.path.join(D, "runs", "depth1_lines"))
    ap.add_argument("--nodes", type=lambda v: [int(x) for x in v.split(",")], default=[1, 5],
                    help="the depth-1 nodes, by root-child index as stage A numbers them (r.1 = 21:1, r.5 = 21:10)")
    ap.add_argument("--nlines", type=int, default=4, help="rule 9's line plus this many minus one by rule 8's score")
    ap.add_argument("--binary", default=os.path.join(R, "latest", "bnbdrv"))
    ap.add_argument("--diag-binary", help="the same sources plus the QAP_XBAR_DEBUG=2 print (prep, first time)")
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--cpus", default="0,1,2,3")
    ap.add_argument("--max-wall-h", type=float, default=5.0)
    ap.add_argument("--stop-after-min", type=float, default=0,
                    help="start no task after this many minutes (running ones finish); 0 = no limit")
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    a = ap.parse_args()
    if a.what == "prep" and not a.diag_binary and not os.path.exists(os.path.join(a.run, "config.json")):
        ap.error("prep needs --diag-binary the first time")
    r = Run(a)
    if a.what == "prep":
        r.prep()
    elif a.what == "kids":
        sys.exit(r.kids())
    else:
        jobs, prep = r.kid_jobs()
        n = sum(1 for j in jobs if os.path.exists(os.path.join(r.run, "kids", j[0], "result.json")))
        print(f"{n} of {len(jobs)} child tasks done")


if __name__ == "__main__":
    main()
