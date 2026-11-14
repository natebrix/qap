#!/usr/bin/env python3
"""The coordinator (coord/, docs/distributed/plan.md) against the serial
solver: with a fixed incumbent every scenario must reproduce the serial
node, fathom, enumerate and Frank-Wolfe totals exactly.

  normal     coordinator + one worker, small budgets (units return, steals)
  kill9      the worker is killed (SIGKILL, solvers included) mid-run; a new
             worker takes the units over after their leases run out
  kill9ckpt  the same with 4-second units, so the killed units resume from
             their uploaded checkpoints
  coordkill  the coordinator is killed mid-run and restarted
  crash_marked, crash_applied, crash_saving
             the coordinator stops (QAP_COORD_CRASH) in its 20th merge right
             after marking a unit merged or after applying it, or in its 3rd
             save between writing the totals and totals.json (saves are once
             per cycle, so a short run has far fewer saves than merges); a new
             one replays and finishes
  drain      the worker gets SIGTERM mid-run (spot notice): running units
             return their open nodes; a new worker finishes
  noub       no incumbent given: the run finds the optimum and publishes it
             (totals then depend on timing; only the objective is checked)
  aws        `normal` on the AWS backend against a local moto server
             (needs boto3 and moto: --python VENV/bin/python)
  tworuns    two runs on one backend at the same time: they share the
             queues, so each run's worker also receives the other's messages;
             both must reproduce the serial totals

Every scenario also checks the public status (coord/publish.py), published
every second: status.json final, with the run's totals; totals.state the same
statistics; report.txt the solver's report; history.json growing in time
(across a coordinator restart too); public/runs/index.json lists the run once
with its final totals. noub runs with --visibility private and must publish
nothing, and is not in the catalog.

  python3 tests/coord_test.py [--instance nug17 --opt 1732] [--only kill9,drain] [--python PY]
"""
import argparse
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
FIXES = ["-F", "0:0"]   # the "fixes" scenario's subtree
FAST = ["--ckpt-upload", "0.5", "--budget", "0.3", "--target-unit-cpu", "0.5", "--steal-after", "0.5", "--heartbeat", "0.5",
        "--checkpoint-every", "0.2", "--poll", "0.3", "--status-every", "5", "--slots", "3", "--publish-every", "1"]


def public_problems(r, private=False):
    """what is wrong with a finished run's public status (coord/publish.py)"""
    import json
    from coord.backend import open_backend
    from coord.common import State
    b = open_backend(r.backend)
    key = lambda name: f"public/runs/{r.run_id}/{name}"
    status = b.get_bytes(key("status.json"))
    catalog = json.loads(b.get_bytes("public/runs/index.json") or b'{"runs": []}')["runs"]
    entry = [c for c in catalog if c["run"] == r.run_id]
    if private:
        published = status is not None or b.list(f"public/runs/{r.run_id}/") or entry
        return ["a private run published"] if published else []
    if status is None:
        return ["no public status.json"]
    s, problems = json.loads(status), []
    want = totals(r.results())
    if not (s["status"] == "finished" and s["final"]):
        problems.append(f"status {s['status']}, final {s['final']}")
    got = (s["nodes"]["evaluated"], s["nodes"]["fathomed"], s["nodes"]["enumerated"])
    if got != want[0] or s["fw_iterations"] != want[1] or s["nodes"]["open"] != 0:
        problems.append(f"status.json nodes {got} fw {s['fw_iterations']} open {s['nodes']['open']}")
    tmp = os.path.join(r.tmp, "public-totals.state")
    if not b.get_file(key("totals.state"), tmp) or tuple(State(tmp).counters[:3]) != want[0]:
        problems.append("totals.state does not hold the run's totals")
    report = (b.get_bytes(key("report.txt")) or b"").decode()
    if totals(report)[:2] != want[:2] or "\nlevel " not in report:
        problems.append("report.txt is not the run's report with the solver's tables")
    if len(entry) != 1 or entry[0]["status"] != "finished" or not entry[0]["final"] or \
            (entry[0]["nodes"], entry[0]["fathomed"]) != want[0][:2]:
        problems.append(f"public/runs/index.json: {len(entry)} entries for the run, not its final totals")
    h = json.loads(b.get_bytes(key("history.json")) or b"[]")
    if len(h) < 2 or h[-1]["nodes"] != want[0][0] or any(x["wall_s"] > y["wall_s"] for x, y in zip(h, h[1:])):
        problems.append(f"history.json: {len(h)} samples, not ending on the totals or not in time order")
    return problems


def totals(text):
    m = re.search(r"Nodes: (\d+) / (\d+) / (\d+)", text)
    f = re.search(r"fw iterations: (\d+)", text)
    b = re.search(r"best objective (\S+)", text)
    return (tuple(map(int, m.groups())) if m else None, int(f.group(1)) if f else None,
            float(b.group(1)) if b else None)


def serial(inst, opt, flags=()):
    r = subprocess.run([os.path.join(ROOT, "latest", "bnbdrv"), "-p", "PARAM/default.par", "-u", str(opt)]
                       + list(flags) + [os.path.join(ROOT, "qaplib", inst + ".dat")], cwd=os.path.join(ROOT, "latest"),
                       capture_output=True, text=True)
    return totals(r.stdout)


class Run:
    def __init__(self, py, backend, run_id, inst, opt, extra=(), lease=4):
        self.py, self.backend, self.run_id = py, backend, run_id
        self.args = ["--backend", backend]
        self.init = ["init", run_id] + self.args + ["--instance", os.path.join(ROOT, "qaplib", inst + ".dat"),
                                                    "--lease", str(lease)] + FAST + list(extra)
        if opt is not None:
            self.init += ["--ub", str(opt)]
        self.procs = []
        self.tmp = tempfile.mkdtemp(prefix="coordtest-")

    def ctl(self, *a, wait=True, log=None, env=None):
        cmd = [self.py, "-m", "coord.qapctl"] + list(a)
        out = open(log, "a") if log else subprocess.DEVNULL
        e = dict(os.environ, **(env or {}))
        p = subprocess.Popen(cmd, cwd=ROOT, stdout=out, stderr=subprocess.STDOUT, start_new_session=True, env=e)
        if wait:
            p.wait()
            return p
        self.procs.append(p)
        return p

    def worker(self, tag):
        return self.ctl("worker", self.run_id, *self.args, "--slots", "3", "--workdir",
                        os.path.join(self.tmp, tag), wait=False, log=os.path.join(self.tmp, tag + ".log"))

    def coordinator(self, tag="coord", env=None):
        return self.ctl("coordinator", self.run_id, *self.args, "--poll", "0.3", wait=False,
                        log=os.path.join(self.tmp, tag + ".log"), env=env)

    def results(self):
        from coord.backend import open_backend
        b = open_backend(self.backend)
        data = b.get_bytes(f"runs/{self.run_id}/results.txt")
        return data.decode() if data else ""

    def cleanup(self):
        for p in self.procs:
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGKILL)


def wait_all(procs, timeout):
    t0 = time.time()
    for p in procs:
        p.wait(timeout=max(1, timeout - (time.time() - t0)))


def two_runs(name, py, base, inst, opt):
    be = f"local:{os.path.join(base, name)}"
    runs = [Run(py, be, f"{name}_{x}", inst, opt) for x in "ab"]
    try:
        for r in runs:
            r.ctl(*r.init, log=os.path.join(r.tmp, "init.log"))
        procs = [p for r in runs for p in (r.coordinator(), r.worker("w1"))]
        try:
            wait_all(procs, 120)
        except subprocess.TimeoutExpired:
            pass                                # a stranded run never finishes: reported below
        got = [totals(r.results()) for r in runs]
        from coord.backend import open_backend
        units = [u for r in runs for u in open_backend(be).all_units(r.run_id)]
        r0 = runs[0]
        r0.stats = {"units": len(units), "retried": sum(u["attempt"] > 1 for u in units),
                    "returned": sum(u.get("outcome") == "returned" for u in units),
                    "resumed": sum(bool(u.get("resumed_from_ckpt")) for u in units)}
        r0.tmp = " and ".join(r.tmp for r in runs)
        return (got[0] if got[0] == got[1] else (None, None, None)), r0
    finally:
        for r in runs:
            r.cleanup()


def scenario(name, py, base, inst, opt, aws_endpoint=None):
    if name == "tworuns":
        return two_runs(name, py, base, inst, opt)
    if name == "aws":
        backend = f"aws:qaptest/{name}?region=eu-north-1&name=qt{int(time.time())}&endpoint={aws_endpoint}"
        r = Run(py, backend, name, inst, opt)
        r.ctl("aws-setup", *r.args)
    else:
        extra = ["--budget", "4", "--target-unit-cpu", "4", "--steal-after", "3"] if name == "kill9ckpt" else []
        extra += ["--visibility", "private"] if name == "noub" else []
        # fixes: a subtree run, -F in the solver flags (the units resume state
        # files, which the solver does not take together with -F)
        extra += ["--extra", " ".join(FIXES)] if name == "fixes" else []
        r = Run(py, f"local:{os.path.join(base, name)}", name, inst, opt if name != "noub" else None, extra)
    try:
        r.ctl(*r.init, log=os.path.join(r.tmp, "init.log"))
        crash = name.startswith("crash_")
        # crash_*: the coordinator exits at a chosen point of its 20th merge
        # (its 3rd save for crash_saving: one save per cycle, after several merges)
        n_crash = 3 if name == "crash_saving" else 20
        c = r.coordinator(env={"QAP_COORD_CRASH": f"{name[6:]}:{n_crash}"} if crash else None)
        w = r.worker("w1")
        if name == "kill9ckpt":
            # kill once a running unit has uploaded a checkpoint
            from coord.backend import open_backend
            b = open_backend(r.backend)
            t_end = time.time() + 30
            while time.time() < t_end and not any(
                    b.exists(f"runs/{r.run_id}/units/{u['uid']}/ckpt.state")
                    for u in b.units_by_status(r.run_id, "leased")):
                time.sleep(0.1)
            time.sleep(0.3)
        if name in ("kill9", "kill9ckpt"):
            if name == "kill9":
                time.sleep(3)
            os.killpg(w.pid, signal.SIGKILL)       # the worker and its solvers
            w.wait()
            time.sleep(1)
            w = r.worker("w2")
        elif name == "coordkill":
            time.sleep(3)
            os.killpg(c.pid, signal.SIGKILL)
            c.wait()
            time.sleep(2)
            c = r.coordinator("coord2")
        elif crash:
            if c.wait(timeout=120) != 9:
                raise RuntimeError(f"{name}: the coordinator did not stop at its crash point")
            time.sleep(1)
            c = r.coordinator("coord2")
        elif name == "drain":
            time.sleep(3)
            os.killpg(w.pid, signal.SIGTERM)       # the supervisor drains; solvers get SIGUSR1
            w.wait(timeout=60)
            w = r.worker("w2")
        wait_all([c, w], 400)
        from coord.backend import open_backend
        units = open_backend(r.backend).all_units(r.run_id)
        r.stats = {"units": len(units), "retried": sum(u["attempt"] > 1 for u in units),
                   "returned": sum(u.get("outcome") == "returned" for u in units),
                   "resumed": sum(bool(u.get("resumed_from_ckpt")) for u in units)}
        r.public = public_problems(r, private=name == "noub")
        return totals(r.results()), r
    finally:
        r.cleanup()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instance", default="nug17")
    ap.add_argument("--opt", type=float, default=1732)
    ap.add_argument("--only", default="normal,kill9,kill9ckpt,coordkill,crash_marked,crash_applied,crash_saving,"
                                      "drain,noub,fixes,aws,tworuns")
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    ref = serial(a.instance, int(a.opt))
    print(f"serial: nodes {ref[0]}, fw {ref[1]}")
    base = tempfile.mkdtemp(prefix="coordtest-backends-")
    server = None
    endpoint = None
    bad = 0
    for name in a.only.split(","):
        if name == "aws" and server is None:
            try:
                import importlib
                subprocess.run([a.python, "-c", "import boto3, moto"], check=True)
                s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
                env = dict(os.environ)
                server = subprocess.Popen([a.python, "-c", f"from moto.server import ThreadedMotoServer; "
                                          f"s = ThreadedMotoServer(ip_address='127.0.0.1', port={port}); s.start(); "
                                          f"import time; time.sleep(3600)"], env=env,
                                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                endpoint = f"http://127.0.0.1:{port}"
                time.sleep(2)
                for k, v in (("AWS_ACCESS_KEY_ID", "test"), ("AWS_SECRET_ACCESS_KEY", "test"),
                             ("NO_PROXY", "127.0.0.1,localhost"), ("no_proxy", "127.0.0.1,localhost")):
                    os.environ[k] = v
                os.environ.pop("AWS_CA_BUNDLE", None)
            except Exception as e:
                print(f"aws: skipped ({e})")
                continue
        t0 = time.time()
        got, r = scenario(name, a.python, base, a.instance, a.opt, endpoint)
        if name == "fixes":
            # the subtree under FIXES against the serial run with the same -F
            # (the bound is the incumbent passed with -u: the subtree may not hold the optimum)
            sub = serial(a.instance, int(a.opt), FIXES)
            ok = got[0] == sub[0] and got[1] == sub[1]
        else:
            ok = got[2] == a.opt if name == "noub" else (got[0] == ref[0] and got[1] == ref[1] and got[2] == a.opt)
        public = getattr(r, "public", [])
        ok = ok and not public
        bad += not ok
        print(f"{name:10s} {'ok ' if ok else 'BAD'} nodes {got[0]} fw {got[1]} best {got[2]}  ({time.time() - t0:.0f} s)"
              + f"  units {r.stats['units']}, retried {r.stats['retried']} ({r.stats['resumed']} from a checkpoint), "
              f"returned {r.stats['returned']}"
              + (f"  public status: {'; '.join(public)}" if public else "")
              + ("" if ok else f"  logs in {r.tmp}"))
    if server:
        server.kill()
    shutil.rmtree(base, ignore_errors=True)
    print("all coordinator scenarios passed" if not bad else f"{bad} scenario(s) failed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
