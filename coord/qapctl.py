#!/usr/bin/env python3
"""qapctl: the coordinator, the worker supervisor and the run controls
(docs/distributed/plan.md). Backends: local:DIR (one machine; SQLite and
files) or aws:BUCKET[/PREFIX][?region=R&name=N] (S3, SQS, DynamoDB).

  python3 -m coord.qapctl init RUN --backend B --instance qaplib/nug20.dat \\
        --params PARAM/default.par --ub 2570 --slots 4 [--extra "-x 0:4:1:300:1"] \\
        [--env DNN_R32=1 ...] [--budget 1800] [--trace]
      init RUN ... --sample-from SRC --sample 3:150,4:150 [--sample-by SRC-dives] [--sample-seed 1]
        (a sampled run, coord/sample.py: no frontier pass; its first units are a sample of the paused
        run SRC's open nodes at those depths, drawn in proportion to the dives run's estimates of
        their subtrees; the same data file and flags as SRC, its own params. Arms with the same
        sample options solve the same nodes: experiments/estimates_nov19/paired.py compares them)
  python3 -m coord.qapctl coordinator RUN --backend B
  python3 -m coord.qapctl worker RUN --backend B --slots 4
  python3 -m coord.qapctl status RUN --backend B
  python3 -m coord.qapctl publish RUN --backend B     (the public status now; publish.py)
  python3 -m coord.qapctl pause|resume RUN --backend B
  python3 -m coord.qapctl local RUN --backend local:DIR ...   (init, then the
        coordinator here and one worker process with --slots slots)
  python3 -m coord.qapctl aws-setup --backend aws:...        (bucket, queues, tables)
  python3 -m coord.qapctl requeue-leased RUN --backend B     (after every worker died, e.g. a
        machine restart: their units go back to the queue at once)
  python3 -m coord.qapctl nodes RUN --backend B --out FILE[.gz]   (a run initialized with --trace:
        its whole node trace, the frontier's and each merged unit's, one line per node)
  python3 -m coord.qapctl export RUN --backend B --out DIR   (a --trace run's whole node and DNN
        traces, nodes.txt.gz and dnn.txt.gz, one line per node and per DNN solve: coord/export.py)
  python3 -m coord.qapctl stop RUN --backend B        (give the run up: status stopped, a final
        public status; its coordinator and workers exit, other runs' workers delete its queued
        units' messages. Terminate its instances apart: infra/phase1/cleanup.sh)
  python3 -m coord.qapctl purge-states RUN --backend B   (a stopped or finished run: delete its
        unit state files, the bulk of its store; ledger, totals, results, logs, traces stay)
  python3 -m coord.qapctl results RUN --backend B   (a paused or stopped solve: its RESULTS block as
        far as it went, results.txt's format, to runs/RUN/results_partial.txt and stdout; the estimate
        flow and qapctl stop write it themselves since Nov 21, this is for runs from before)
  python3 -m coord.qapctl dives SRC --as NEW --backend B [--dive-b 0.4] [--budget S]
        (SRC paused, its units returned and merged: run NEW, one unit of survivor dives per queued
        unit of SRC, S wall seconds each; then the coordinator and workers of NEW: coord/dives.py)
  python3 -m coord.qapctl estimate NEW --backend B [--out DIR]   (SRC's work so far plus the
        dives' estimate of the rest, by depth: estimate.txt and estimate.json)
  python3 -m coord.qapctl estimate-flow RUN --backend B --solve-wall S --dive-budget T [--dive-b 0.4]
        [--decide 0.15]
        (in place of `coordinator RUN`: the solve for S s of its wall time, then pause, dives of T s a
        unit as RUN-dives, coordinated to the end; workers started with `worker RUN --follow` go
        from the one to the other. --decide F: the decision point, once a tenth of the dives units
        have merged: when the estimator's CPU to the end of the dives would be over F of the
        estimated rest of the tree, the dives stop and the solve goes on to the end, its workers
        coming back to it; 0: off)

Unit sizing and timing (init): --budget (CPU s a unit may run before it
returns its open nodes, default 1800), --target-unit-cpu (estimated CPU per
new unit when splitting returned nodes, default 1800), --reserve (keep this
many x slots units queued, default 2), --steal-after (s a unit must have run
before an idle slot may ask it to return, default 30), --checkpoint-every
(the solver's -k, default 60), --lease (s, default 1800), --heartbeat (s,
default 60). The public status (coord/publish.py): --publish-every (s, default
300; 0: only at the end), --visibility public|private (private: never).
"""
import argparse
import json
import os
import shlex
import subprocess
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "coord"

from coord.backend import open_backend  # noqa: E402
from coord.common import ROOT, SOLVER, run_key, unit_key  # noqa: E402
from coord.coordinator import Coordinator, drop_queues, ensure_queues, init_run  # noqa: E402
from coord.worker import Worker  # noqa: E402


def log(msg):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def cfg_from(a):
    return {"slots": a.slots, "budget": a.budget, "target_unit_cpu": a.target_unit_cpu, "reserve": a.reserve,
            "steal_after": a.steal_after, "checkpoint_every": a.checkpoint_every, "lease_s": a.lease,
            "heartbeat_s": a.heartbeat, "max_attempts": a.max_attempts, "trace": a.trace,
            "large_depth": a.large_depth, "status_every": a.status_every, "poll_s": a.poll,
            "ckpt_upload_s": a.ckpt_upload, "publish_every": a.publish_every}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["init", "coordinator", "worker", "status", "pause", "resume", "local",
                                       "aws-setup", "requeue-leased", "publish", "nodes", "export", "stop",
                                       "purge-states", "dives", "estimate", "estimate-flow", "readout", "dives-export", "results"])
    p.add_argument("run", nargs="?")
    p.add_argument("--backend", required=True)
    p.add_argument("--instance")
    p.add_argument("--params", default="PARAM/default.par", help="relative to latest/")
    p.add_argument("--extra", default="", help="solver flags, e.g. '-F 39:27 -x 0:4:1:300:1'")
    p.add_argument("--env", action="append", default=[],
                   help="K=V for the solver's environment (repeat): experiment switches and diagnostics; "
                        "policy settings are keys in the params file (RUNNING.md, Settings)")
    p.add_argument("--ub", type=float)
    p.add_argument("--frontier-depth", type=int, default=0)
    p.add_argument("--sample-from", metavar="SRC", help="init: a sampled run of the paused run SRC's open nodes")
    p.add_argument("--sample", metavar="D:N,...", help="init --sample-from: how many open nodes at each depth")
    p.add_argument("--sample-by", metavar="DIVES_RUN", help="init --sample-from: the dives run whose estimates size "
                                                             "the nodes (default SRC's own)")
    p.add_argument("--sample-seed", type=int, default=1)
    p.add_argument("--sample-mix", type=float, default=0.2,
                   help="init --sample-from: the uniform share mixed into the size-proportional draw")
    p.add_argument("--sample-power", type=float, default=0.5,
                   help="init --sample-from: draw in proportion to the estimate to this power (0.5: sqrt, the "
                        "least variance of a ratio between arms at a fixed cost; 0: uniform)")
    p.add_argument("--solver", default=SOLVER)
    p.add_argument("--slots", type=int, default=os.cpu_count() or 1)
    p.add_argument("--budget", type=float, default=1800)
    p.add_argument("--target-unit-cpu", type=float, default=1800)
    p.add_argument("--reserve", type=float, default=2.0)
    p.add_argument("--steal-after", type=float, default=30)
    p.add_argument("--checkpoint-every", type=float, default=60)
    p.add_argument("--lease", type=float, default=1800)
    p.add_argument("--heartbeat", type=float, default=60)
    p.add_argument("--max-attempts", type=int, default=3)
    p.add_argument("--ckpt-upload", type=float, default=300, help="s between uploads of a running unit's checkpoint")
    p.add_argument("--large-depth", type=int, default=4)
    p.add_argument("--status-every", type=float, default=30)
    p.add_argument("--publish-every", type=float, default=300,
                   help="s between public status bundles (coord/publish.py); 0: only at the end")
    p.add_argument("--visibility", choices=["public", "private"], default="public",
                   help="private: publish no public status for this run")
    p.add_argument("--poll", type=float, default=2.0)
    p.add_argument("--trace", action="store_true",
                   help="QAP_DNN_TRACE and QAP_NODE_TRACE per unit, uploaded with the result as trace.txt and "
                        "nodes.txt (the frontier pass's under runs/RUN/frontier/)")
    p.add_argument("--workdir")
    p.add_argument("--out", help="nodes: the file to write (.gz: gzipped); export: the directory to write "
                                 "nodes.txt.gz, dnn.txt.gz and manifest.json in (coord/export.py)")
    p.add_argument("--as", dest="as_run", help="dives: the new run's name")
    p.add_argument("--dive-b", default="0.4",
                   help="dives: b of the survivor weights exp(b x), x = shortfall / the depth's median (0: uniform); "
                        "kQ: Knuth's dives, one child a level by g^Q of its QPB bound (bushy trees)")
    p.add_argument("--dive-sample", type=float, default=1.0,
                   help="dives: the share of the open nodes sampled (each weighted 1/pi; 1: all)")
    p.add_argument("--dive-weights", choices=["unif", "qpb"], default="unif",
                   help="dives, --dive-sample: the open nodes' sampling weights (qpb: g^2 of their own QPB bound)")
    p.add_argument("--follow", action="store_true",
                   help="worker: when the run pauses for an estimate (estimate-flow), go on with its dives run")
    p.add_argument("--save", action="store_true", help="readout: also write it to the store as runs/RUN/readout.md")
    p.add_argument("--dive-env", action="append", default=[], metavar="QAP_EST_X=V",
                   help="dives, estimate-flow: another estimator setting for every dives unit "
                        "(QAP_EST_MIN_DIVES, QAP_EST_SAMPLE_DEPTHS, ...; repeatable)")
    p.add_argument("--solve-wall", type=float, help="estimate-flow: s of the run's wall time the solve gets")
    p.add_argument("--dive-budget", type=float, help="estimate-flow: wall s of survivor dives per dives unit")
    p.add_argument("--decide-every", type=float, default=120.0,
                   help="estimate-flow: seconds between decision checks (each reads the merged units' dives)")
    p.add_argument("--decide", type=float, default=None,
                   help="estimate-flow: stop the dives and solve when they would cost over this share of the "
                        "estimated rest (default 0.15, coord/dives.py DECIDE; 0: off)")
    p.add_argument("--exit-when-idle", type=float)
    p.add_argument("--max-wall", type=float)
    a = p.parse_args(argv)
    b = open_backend(a.backend)

    if a.command == "aws-setup":
        b.setup()
        log("aws: bucket, queues and tables ready")
        return 0
    if not a.run:
        p.error("RUN is required")

    if a.command in ("init", "local"):
        env = dict(kv.split("=", 1) for kv in a.env)
        sample = None
        if a.sample_from:
            from coord.sample import parse_plan
            if not a.sample:
                p.error("--sample-from needs --sample D:N,...")
            sample = {"src": a.sample_from, "by": a.sample_by, "plan": parse_plan(a.sample),
                      "seed": a.sample_seed, "mix": a.sample_mix, "power": a.sample_power}
        init_run(b, a.run, a.instance, a.params, shlex.split(a.extra), env, a.ub, cfg_from(a),
                 a.solver, a.frontier_depth, log=log, visibility=a.visibility, sample=sample)
        if a.command == "init":
            return 0
        # one worker process with all the slots, the coordinator here
        wcmd = [sys.executable, "-m", "coord.qapctl", "worker", a.run, "--backend", a.backend,
                "--slots", str(a.slots), "--solver", a.solver]
        if a.workdir:
            wcmd += ["--workdir", a.workdir]
        wproc = subprocess.Popen(wcmd, cwd=ROOT)
        out = None
        try:
            out = Coordinator(b, a.run, log=log).loop(poll=a.poll, max_wall=a.max_wall)
        finally:
            if out is None:
                wproc.terminate()           # drain: open nodes go back as returned units
            try:
                wproc.wait(timeout=180)
            except subprocess.TimeoutExpired:
                wproc.kill()
        return 0 if out is not None else 3

    if a.command == "coordinator":
        out = Coordinator(b, a.run, log=log).loop(poll=a.poll, max_wall=a.max_wall)
        return 0 if out is not None else 3
    if a.command == "worker":
        run = a.run
        while True:
            w = Worker(b, run, a.slots, solver=a.solver, workdir=a.workdir, log=log,
                       exit_when_idle=a.exit_when_idle, follow=a.follow)
            w.run_forever()
            r = b.get_run(run) or {}
            nxt = r.get("dives_run")
            dec = r.get("decision") if r.get("config", {}).get("kind") == "dives" else None
            if a.follow and not w.draining.is_set() and r.get("status") == "stopped" and dec \
                    and dec.get("verdict") == "solve":
                # the estimate's decision point stopped the dives for the solve (coord/dives.py): back
                # to the source run, which the estimate-flow resumed before stopping this one
                src = r["config"]["source_run"]
                for _ in range(720):
                    if (b.get_run(src) or {}).get("status") == "running":
                        break
                    time.sleep(5)
                else:
                    log(f"worker: {run} stopped for the solve, but {src} did not resume in an hour")
                    return 1
                log(f"worker: {run}'s estimate decided to solve: back to {src}")
                run = src
                continue
            if not a.follow or w.draining.is_set() or r.get("status") != "paused" or not nxt:
                return 0
            # --follow: the run paused for an estimate (qapctl estimate-flow): on to its dives run,
            # which exists once the run's units have returned and merged
            for _ in range(720):
                if b.get_run(nxt) is not None:
                    break
                time.sleep(5)
            else:
                log(f"worker: {run} paused for {nxt}, which did not appear in an hour")
                return 1
            log(f"worker: {run} paused for its estimate: on to {nxt}")
            run = nxt
    if a.command == "dives-export":
        from coord.dives import export_dives
        if not a.out:
            p.error("dives-export needs --out DIR")
        return export_dives(b, a.run, a.out, log=log)
    if a.command == "readout":
        # where the run stands and what it found (coord/readout.py); --save writes runs/RUN/readout.md
        from coord.readout import readout, write_readout
        text = write_readout(b, a.run, log=log) if a.save else readout(b, a.run)
        if a.out:
            with open(a.out, "w") as f:
                f.write(text or "")
        print(text)
        return 0
    if a.command == "estimate-flow":
        from coord.dives import estimate_flow
        if a.solve_wall is None or a.dive_budget is None:
            p.error("estimate-flow: --solve-wall and --dive-budget are required")
        from coord.dives import DECIDE
        out = estimate_flow(b, a.run, a.solve_wall, a.dive_budget, dive_b=a.dive_b, poll=a.poll, log=log,
                            sample=a.dive_sample, weights=a.dive_weights, env=dive_env(a, p),
                            decide=DECIDE if a.decide is None else a.decide, check_every=a.decide_every)
        return 0 if out is not None else 3
    if a.command == "status":
        s = b.get_bytes(run_key(a.run, "status.json"))
        print(s.decode() if s else json.dumps(b.get_run(a.run), indent=1))
        print("units:", b.count_by_status(a.run))
        return 0
    if a.command == "publish":
        c = Coordinator(b, a.run, log=log)
        if c.publisher.publish(b.count_by_status(a.run), final=c.run.get("status") in ("finished", "stopped")):
            log(f"run {a.run}: public status written under public/runs/{a.run}/")
        else:
            log(f"run {a.run} is private: nothing published")
        return 0
    if a.command == "nodes":
        if not a.out:
            p.error("nodes needs --out FILE")
        return write_nodes(b, a.run, a.out, log=log)
    if a.command == "dives":
        from coord.dives import init_dives
        if not a.as_run:
            p.error("dives: --as NEW is required")
        init_dives(b, a.run, a.as_run, dive_b=a.dive_b, budget=a.budget, log=log, sample=a.dive_sample,
                   weights=a.dive_weights, env=dive_env(a, p))
        return 0
    if a.command == "estimate":
        from coord.dives import estimate_run
        estimate_run(b, a.run, out=a.out, log=log)
        return 0
    if a.command == "export":
        if not a.out:
            p.error("export needs --out DIR")
        from .export import export_run
        return export_run(b, a.run, a.out, log=log)
    if a.command == "stop":
        return stop_run(b, a.run, log=log)
    if a.command == "purge-states":
        return purge_states(b, a.run, log=log)
    if a.command == "results":
        return partial_results(b, a.run, log=log)
    if a.command == "requeue-leased":
        Coordinator(b, a.run, log=log).requeue_leased()
        return 0
    if a.command in ("pause", "resume"):
        b.update_run(a.run, {"status": "paused" if a.command == "pause" else "running"})
        if a.command == "resume":
            # each queued unit sent again: the run's queues may be new (a run from before per-run
            # queues) or its messages past SQS's 14-day retention; a second copy is harmless
            n = ensure_queues(b, a.run, log=log, resend=True)
            log(f"run {a.run}: resumed ({n} queued units sent again)")
        else:
            log(f"run {a.run}: paused")
        return 0


def dive_env(a, p):
    env = {}
    for kv in a.dive_env:
        k, sep, v = kv.partition("=")
        if not sep:
            p.error(f"--dive-env {kv}: K=V")
        env[k] = v
    return env


def stop_run(b, run_id, log=print):
    """give up a run for good: status stopped, its wall time the last the
    coordinator recorded, and a final public status (the catalog and the site
    show it stopped, no longer live). Its coordinator and workers exit at their
    next look at the run, and its queues are deleted with the messages in them.
    Its instances are terminated separately (infra/phase1/cleanup.sh)."""
    run = b.get_run(run_id)
    if run is None:
        log(f"no run {run_id}")
        return 1
    if run.get("status") == "finished":
        log(f"run {run_id} is finished: nothing to stop")
        return 1
    last = json.loads(b.get_bytes(run_key(run_id, "status.json")) or b"{}")
    fields = {"status": "stopped", "stopped_at": time.time()}
    if run.get("status") != "stopped":   # a second stop keeps the first one's times
        fields.update({"wall_s": last.get("wall_s"), "finished_at": last.get("updated")})
    b.update_run(run_id, fields)
    drop_queues(b, run_id)   # its coordinator and workers exit at their next look; nothing runs it again
    from coord.readout import write_readout
    write_readout(b, run_id, log=log)
    # the saved totals, without replaying the merges after them: a replay fetches
    # those units' output states, gigabytes each in a run like tho40-a
    c = Coordinator(b, run_id, log=log, replay=False)
    if c.publisher.publish(b.count_by_status(run_id), final=True):
        log(f"run {run_id}: stopped; final public status written under public/runs/{run_id}/")
    else:
        log(f"run {run_id}: stopped (private: nothing published)")
    # the stopped solve's RESULTS block, as far as it went (collect.sh puts it in git)
    try:
        c.write_partial("the solve stopped (qapctl stop)")
    except Exception as e:
        log(f"run {run_id}: results_partial.txt failed: {e!r}")
    return 0


def partial_results(b, run_id, log=print):
    """results_partial.txt for a solve its operation left unfinished (a run paused or stopped before
    Nov 21, when the estimate flow and qapctl stop began to write it); a finished run has results.txt"""
    run = b.get_run(run_id)
    if run is None:
        log(f"no run {run_id}")
        return 1
    st = run.get("status")
    if st == "finished":
        log(f"run {run_id} is finished: its RESULTS block is results.txt")
        return 1
    if st not in ("paused", "stopped"):
        log(f"run {run_id} is {st}: results_partial.txt is for a paused or stopped solve")
        return 1
    why = (f"the solve paused for the estimate {run['dives_run']}" if st == "paused" and run.get("dives_run")
           else "the solve stopped (qapctl stop)" if st == "stopped" else "the solve paused (qapctl pause)")
    print(Coordinator(b, run_id, log=log, replay=False).write_partial(why))
    return 0


def purge_states(b, run_id, log=print):
    """delete a stopped or finished run's unit state files (in.state, each
    attempt's out.state, checkpoints) and parent states (warm/): the bulk of
    its store, which nothing reads once the run is over. The ledger, the
    totals, results, logs and node traces stay."""
    run = b.get_run(run_id)
    if run is None or run.get("status") not in ("stopped", "finished"):
        log(f"run {run_id} is {run and run.get('status')}: purge-states is for a stopped or finished run")
        return 1
    keys = [k for k in b.list(run_key(run_id, "units") + "/")
            if k.endswith(".state") or k.endswith("ckpt.nodes.txt") or k.endswith("ckpt.trace.txt")]
    keys += b.list(run_key(run_id, "warm") + "/")   # the parent states (coord/warm.py)
    b.delete_many(keys)
    log(f"run {run_id}: {len(keys)} unit state files and parent states deleted")
    return 0


def write_nodes(b, run_id, out, log=print):
    """a --trace run's tree, one QAP_NODE_TRACE line per node (latest/bnb.C): the
    frontier pass's lines, then each merged unit's, from the attempt the ledger
    records (coord/README.md, Traces: each node is in exactly one of them).
    Units not merged yet are left out and counted; 1 if the run has no trace."""
    import gzip
    run = b.get_run(run_id)
    if run is None or not run["config"].get("trace"):
        log(f"run {run_id}: no frontier/nodes.txt (not initialized with --trace?)")
        return 1
    # a frontier pass that bounds no node (a -F run at depth 0) leaves no trace
    frontier = b.get_bytes(run_key(run_id, "frontier", "nodes.txt")) or b""
    units = sorted(b.all_units(run_id), key=lambda u: (u.get("merged_seq") or 0, u["uid"]))
    merged = [u for u in units if u.get("status") == "merged"]
    tmp = out + ".tmp"
    lines = frontier.count(b"\n")
    missing = []
    with (gzip.open(tmp, "wb") if out.endswith(".gz") else open(tmp, "wb")) as f:
        f.write(frontier)
        for u in merged:
            data = b.get_bytes(unit_key(run_id, u["uid"], f"a{u['attempt']}", "nodes.txt"))
            if data is None:
                missing.append(u["uid"])
                continue
            f.write(data)
            lines += data.count(b"\n")
    os.replace(tmp, out)
    log(f"run {run_id}: {lines} node lines (frontier and {len(merged) - len(missing)} merged units) in {out}; "
        f"{len(units) - len(merged)} units not merged{', no nodes.txt for ' + ' '.join(missing) if missing else ''}")
    return 2 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
