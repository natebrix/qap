#!/usr/bin/env python3
"""Solve one QAPLIB instance with several solver processes.

  python3 bench/parallel.py qaplib/nug20.dat --jobs 4 --ub 2570
  python3 bench/parallel.py qaplib/nug20.dat --jobs 4            # no upper bound

Uses the solver's state files (latest/state.h):

 1. A frontier pass runs the search to depth D (-d D -o F); the open nodes
    left in F are split into one-node tasks. Its statistics are the base.
 2. Up to --jobs workers run at a time, each resuming a task file holding a
    chunk of nodes (-r T -o T). A worker that finishes contributes its
    statistics. When a worker is idle and the pool is empty, the busiest
    worker is sent SIGUSR1: it saves its remaining nodes and stops at the
    next node (exit 3), and those nodes go back to the pool. So a long
    subtree does not hold up the end.
 3. New tasks start from the best incumbent found so far.
 4. Statistics from all pieces are summed and printed in the solver's own
    RESULTS format (bench.py parses it).
 5. With --checkpoint F, every --checkpoint-every seconds all workers are
    stopped and their nodes merged into F, an ordinary state file holding
    every open node and the totals so far. --resume F continues from it (or
    finish it serially: bnbdrv -r F).

With a fixed upper bound (e.g. the known optimum) the work done on a node
does not depend on which process does it, so node, fathom and FW totals are
identical to a serial run. Without one, workers only learn of a better
incumbent when they start a new task, so totals can differ.
"""

import argparse
import collections
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "latest", "bnbdrv")


# ---------------------------------------------------------------- state files

class Warm:
    """A parent's Frank-Wolfe point (a "warm" line), shared by the nodes that
    start from it; written once per file, wherever its nodes go."""
    __slots__ = ("text",)

    def __init__(self, text):
        self.text = text


class State:
    """A parsed state file (format in latest/state.h). nodes is a list of
    (node line without its warm start index, Warm or None)."""

    ARRAYS = ["n_level", "f_level", "l_level", "e_level", "t_level", "n_strat",
              "t_strat", "rg_sum", "rg_sumsq"]

    def __init__(self, path):
        with open(path) as f:
            lines = [l.split() for l in f.read().splitlines() if l.strip()]
        it = iter(lines)
        self.magic = next(it)
        hdr = next(it)
        self.n, self.levels, self.strategies = int(hdr[1]), int(hdr[3]), int(hdr[5])
        inc = next(it)
        self.obj, self.found, self.perm = float(inc[1]), int(inc[2]), inc[3:]
        self.root_bound = next(it)[1]
        c = next(it)
        self.counters = [int(x) for x in c[1:5]]
        t = next(it)
        self.times = [float(x) for x in t[1:3]]
        self.arrays = collections.OrderedDict()
        line = next(it)
        while line[0] not in ("warm_starts", "nodes"):
            conv = float if line[0] in ("t_level", "t_strat", "rg_sum", "rg_sumsq") else int
            self.arrays[line[0]] = [conv(x) for x in line[1:]]
            line = next(it)
        warms = []
        if line[0] == "warm_starts":            # version 2
            warms = [Warm(" ".join(next(it)[1:])) for _ in range(int(line[1]))]
            line = next(it)
        has_warm = int(self.magic[1]) >= 2
        self.nodes = []
        for _ in range(int(line[1])):
            f = next(it)
            w = int(f.pop()) if has_warm else -1
            self.nodes.append((" ".join(f), warms[w] if w >= 0 else None))
        assert next(it)[0] == "end"

    def zero_stats(self):
        self.counters = [0, 0, 0, 0]
        self.times = [0.0, 0.0]
        for k, v in self.arrays.items():
            self.arrays[k] = [type(x)(0) for x in v]

    def add_stats(self, other):
        self.counters = [a + b for a, b in zip(self.counters, other.counters)]
        self.times = [a + b for a, b in zip(self.times, other.times)]
        for k in self.arrays:
            self.arrays[k] = [a + b for a, b in zip(self.arrays[k], other.arrays[k])]

    def write(self, path, nodes=None):
        nodes = self.nodes if nodes is None else nodes
        out = ["QAP-STATE 2",
               f"n {self.n} levels {self.levels} strategies {self.strategies}",
               " ".join(["incumbent", repr(self.obj), str(self.found)] + list(self.perm)),
               f"root_bound {self.root_bound}",
               "counters " + " ".join(map(str, self.counters)),
               "times " + " ".join(repr(x) for x in self.times)]
        for k, v in self.arrays.items():
            out.append(" ".join([k] + [repr(x) if isinstance(x, float) else str(x) for x in v]))
        ids = {}
        for _, w in nodes:
            if w is not None and id(w) not in ids:
                ids[id(w)] = (len(ids), w)
        out.append(f"warm_starts {len(ids)}")
        out += ["warm " + w.text for _, w in sorted(ids.values(), key=lambda x: x[0])]
        out.append(f"nodes {len(nodes)}")
        out += [f"{t} {ids[id(w)][0] if w is not None else -1}" for t, w in nodes]
        out.append("end")
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            f.write("\n".join(out) + "\n")
        os.replace(tmp, path)


# ---------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("instance", help="QAPLIB .dat file")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    p.add_argument("--ub", help="initial upper bound (default: none)")
    p.add_argument("--params", default="PARAM/default.par", help="relative to latest/")
    p.add_argument("--depth", type=int, default=0,
                   help="frontier depth for the first split (default 0: only the root "
                        "is solved serially; deeper frontiers serialize the expensive "
                        "lookahead levels, e.g. nug17 x3.74 CPU/wall at 0, x2.0 at 2)")
    p.add_argument("--slice", type=float, default=0,
                   help="optional CPU-seconds limit per task (default: none; idle "
                        "workers take work from busy ones instead)")
    p.add_argument("--steal-after", type=float, default=1.0,
                   help="only take work from a worker that has run this many "
                        "seconds (default 1)")
    p.add_argument("--solver", default=SOLVER)
    p.add_argument("--extra", default="", help="extra solver flags, e.g. '-w 5'")
    p.add_argument("--checkpoint", help="write a merged state file here every "
                   "--checkpoint-every seconds (all open nodes plus totals)")
    p.add_argument("--checkpoint-every", type=float, default=600,
                   help="wall seconds between checkpoints (default 600)")
    p.add_argument("--resume", help="continue from a checkpoint written by --checkpoint "
                   "(skips the frontier pass)")
    p.add_argument("--workdir", help="keep task files here (default: a temp dir, removed)")
    args = p.parse_args()

    solver = os.path.abspath(args.solver)
    sdir = os.path.dirname(solver)
    inst = os.path.abspath(args.instance)
    work = args.workdir or tempfile.mkdtemp(prefix="qap-par-")
    os.makedirs(work, exist_ok=True)
    base_cmd = [solver, "-p", args.params] + args.extra.split()
    if args.ub is not None:
        base_cmd += ["-u", args.ub]
    wall0 = time.perf_counter()

    # 1. frontier (or a parallel checkpoint, which is an ordinary state file)
    if args.resume:
        frontier = os.path.abspath(args.resume)
        print(f"parallel: resuming from {frontier}", flush=True)
    else:
        frontier = os.path.join(work, "frontier.state")
        r = subprocess.run(base_cmd + ["-d", str(args.depth), "-o", frontier, inst],
                           cwd=sdir, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"frontier pass failed (rc={r.returncode}):\n{r.stdout[-2000:]}{r.stderr}")
        # the frontier pass's output first: bench.py takes the last "New objective
        # value" line as the answer, so workers' improvements must come after it
        print(r.stdout.split("  ---- RESULTS ----")[0].rstrip(), flush=True)
    total = State(frontier)
    best_obj, best_perm, best_found = total.obj, total.perm, total.found
    template = State(frontier)
    template.zero_stats()

    pool = collections.deque(total.nodes)
    cpu_before = total.times[0]          # CPU already spent before a --resume
    print(f"parallel: frontier at depth {args.depth} has {len(pool)} open nodes; "
          f"{args.jobs} workers", flush=True)

    # 2-3. workers
    running = {}                    # pid -> [task path, start time, stop requested?, Popen]
    task_id = 0
    tasks_run = steals = 0

    def launch(free_slots):
        nonlocal task_id
        # chunks shrink as the pool does; at least one node per task
        size = max(1, len(pool) // (4 * args.jobs))
        nodes = [pool.popleft() for _ in range(min(size, len(pool)))]
        task_id += 1
        path = os.path.join(work, f"task{task_id}.state")
        template.obj, template.perm, template.found = best_obj, best_perm, best_found
        template.write(path, nodes)
        cmd = base_cmd + ["-r", path, "-o", path]
        if args.slice > 0:
            cmd += ["-t", str(args.slice)]
        log = open(path + ".log", "w")
        proc = subprocess.Popen(cmd + [inst], cwd=sdir, stdout=log, stderr=subprocess.STDOUT)
        log.close()
        # keep the Popen object: if it is garbage-collected, subprocess reaps
        # the child itself and os.waitpid below never sees it
        running[proc.pid] = [path, time.perf_counter(), False, proc]

    def steal():
        """Ask the longest-running worker that has run --steal-after seconds to
        save its remaining nodes and stop (SIGUSR1)."""
        nonlocal steals
        now = time.perf_counter()
        cand = [(t0, pid) for pid, (path, t0, asked, proc) in running.items()
                if not asked and now - t0 >= args.steal_after]
        if not cand:
            return False
        pid = min(cand)[1]
        os.kill(pid, signal.SIGUSR1)
        running[pid][2] = True
        steals += 1
        return True

    def stop_workers(signum, frame):
        for v in running.values():
            v[3].kill()
        for v in running.values():
            v[3].wait()
        sys.exit(f"interrupted (signal {signum}); workers stopped"
                 + (f"; resume with --resume {args.checkpoint}" if args.checkpoint else ""))
    signal.signal(signal.SIGTERM, stop_workers)
    signal.signal(signal.SIGINT, stop_workers)

    next_checkpoint = time.perf_counter() + args.checkpoint_every
    checkpointing = False

    while pool or running:
        # Parallel checkpoint: stop every worker (each saves its nodes), merge,
        # write one state file holding the totals and every open node, go on.
        if args.checkpoint and not checkpointing and time.perf_counter() >= next_checkpoint:
            checkpointing = True
            for pid, v in running.items():
                if not v[2]:
                    os.kill(pid, signal.SIGUSR1)
                    v[2] = True
        if checkpointing and not running:
            total.obj, total.perm, total.found = best_obj, best_perm, best_found
            total.write(os.path.abspath(args.checkpoint), list(pool))
            print(f"parallel: checkpoint with {len(pool)} open nodes, {total.counters[0]} "
                  f"nodes so far -> {args.checkpoint}", flush=True)
            checkpointing = False
            next_checkpoint = time.perf_counter() + args.checkpoint_every
        while pool and len(running) < args.jobs and not checkpointing:
            launch(args.jobs - len(running))
        if not running:
            continue
        # idle slots and nothing queued: take work back from a busy worker
        idle = args.jobs - len(running)
        asked = sum(1 for v in running.values() if v[2])
        while not checkpointing and idle > asked and steal():
            asked += 1
        # wait for a worker to finish (poll only while a steal is still owed)
        if idle > asked or (args.checkpoint and not checkpointing):
            pid, status = os.waitpid(-1, os.WNOHANG)
            if pid == 0:
                time.sleep(0.02)
                continue
        else:
            pid, status = os.waitpid(-1, 0)
        if pid not in running:
            continue
        rc = os.waitstatus_to_exitcode(status)
        path = running.pop(pid)[0]
        tasks_run += 1
        if rc not in (0, 3):
            sys.exit(f"worker failed (rc={rc}) on {path}; see {path}.log")
        st = State(path)
        total.add_stats(st)
        if st.obj < best_obj:
            best_obj, best_perm, best_found = st.obj, st.perm, st.found
            print(f"New objective value: {st.obj:.17g}", flush=True)
        if rc == 3:
            pool.extend(st.nodes)          # unfinished nodes go back to the pool
        os.remove(path)
        os.remove(path + ".log")

    wall = time.perf_counter() - wall0

    # 4. report in the solver's format
    total.obj, total.perm, total.found = best_obj, best_perm, best_found
    node, fathom, enumerate_, fw = total.counters
    print("  ---- RESULTS ----\n")
    print(f"    Total Time: {total.times[0]:.6f} seconds ")
    print(f"         Nodes: {node} / {fathom} / {enumerate_}")
    print(f"    root bound: {total.root_bound}")
    print(f" fw iterations: {fw}")
    session_cpu = total.times[0] - cpu_before
    print(f"\nparallel: {args.jobs} workers, {tasks_run} tasks, {steals} steals, wall {wall:.2f} s, "
          f"CPU {session_cpu:.2f} s (x{session_cpu / wall:.2f})"
          + (f", {total.times[0]:.2f} s including before the resume" if cpu_before else "") + ", "
          f"best objective {best_obj:.17g}")
    total.write(os.path.join(work, "final.state"), [])
    if not args.workdir:
        shutil.rmtree(work)


if __name__ == "__main__":
    main()
