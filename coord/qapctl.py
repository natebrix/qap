#!/usr/bin/env python3
"""qapctl: the coordinator, the worker supervisor and the run controls
(docs/distributed/plan.md). Backends: local:DIR (one machine; SQLite and
files) or aws:BUCKET[/PREFIX][?region=R&name=N] (S3, SQS, DynamoDB).

  python3 -m coord.qapctl init RUN --backend B --instance qaplib/nug20.dat \\
        --params PARAM/default.par --ub 2570 --slots 4 [--extra "-x 0:4:1:300:1"] \\
        [--env DNN_R32=1 ...] [--budget 1800] [--trace]
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
from coord.common import ROOT, SOLVER, run_key  # noqa: E402
from coord.coordinator import Coordinator, init_run  # noqa: E402
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
                                       "aws-setup", "requeue-leased", "publish"])
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
        init_run(b, a.run, a.instance, a.params, shlex.split(a.extra), env, a.ub, cfg_from(a),
                 a.solver, a.frontier_depth, log=log, visibility=a.visibility)
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
        Worker(b, a.run, a.slots, solver=a.solver, workdir=a.workdir, log=log,
               exit_when_idle=a.exit_when_idle).run_forever()
        return 0
    if a.command == "status":
        s = b.get_bytes(run_key(a.run, "status.json"))
        print(s.decode() if s else json.dumps(b.get_run(a.run), indent=1))
        print("units:", b.count_by_status(a.run))
        return 0
    if a.command == "publish":
        c = Coordinator(b, a.run, log=log)
        if c.publisher.publish(b.count_by_status(a.run), final=c.run.get("status") == "finished"):
            log(f"run {a.run}: public status written under public/runs/{a.run}/")
        else:
            log(f"run {a.run} is private: nothing published")
        return 0
    if a.command == "requeue-leased":
        Coordinator(b, a.run, log=log).requeue_leased()
        return 0
    if a.command in ("pause", "resume"):
        b.update_run(a.run, {"status": "paused" if a.command == "pause" else "running"})
        log(f"run {a.run}: {a.command}d")
        return 0


if __name__ == "__main__":
    sys.exit(main())
