"""The worker supervisor (docs/distributed/plan.md 4.1-4.3, 8.4): one per
machine, `slots` solver processes, each slot a thread that

  receives a unit from the queue (the message's visibility is the lease),
  leases it in the ledger (conditional: queued, or a lease that ran out),
  runs  bnbdrv -p P <flags> -u INC -U inc.txt -r in.state -o out.state -t BUDGET -k CKPT  instance,
  heartbeats every heartbeat_s while it runs: extends the message and the
    lease, rewrites inc.txt and sends SIGUSR2 when a better incumbent is
    published, sends SIGUSR1 (save the open nodes, stop at the next node)
    when the coordinator asks the unit to return or the supervisor drains,
  uploads out.state, the solver log, the node trace and result.json under the
    attempt's prefix, records the outcome (conditional on still holding the
    lease) and deletes the message.

SIGTERM or SIGINT drains: every running solver is asked to stop at its next
node; one still inside a node after `grace` seconds is killed and its last
checkpoint (-k) is returned instead (it holds every node not yet finished);
a unit with no checkpoint is released for another worker.
"""
import json
import os
import shutil
import signal
import socket
import subprocess
import threading
import time

from .common import LANES, ROOT, SOLVER, State, run_key, unit_key

# every run on a backend shares its queues: a worker that receives another
# run's message hides it this long (and waits a poll) and leaves it for that
# run's workers
FOREIGN_HIDE_S = 1


class Worker:
    def __init__(self, backend, run_id, slots, solver=SOLVER, workdir=None, log=print,
                 exit_when_idle=None, grace=90.0):
        self.b = backend
        self.run_id = run_id
        self.slots = slots
        self.solver = os.path.abspath(solver)
        self.log = log
        self.exit_when_idle = exit_when_idle
        self.grace = grace
        self.workdir = os.path.abspath(workdir or os.path.join("/tmp", f"qapworker-{run_id}-{os.getpid()}"))
        os.makedirs(self.workdir, exist_ok=True)
        self.name = f"{socket.gethostname()}/{os.getpid()}"
        self.draining = threading.Event()
        self.lock = threading.Lock()
        self.run = self.b.get_run(run_id)
        if self.run is None:
            raise SystemExit(f"no run {run_id}")
        self.cfg = self.run["config"]
        self.sol = self.run["solver"]
        # the instance and the parameters, as the run stored them
        self.inst = os.path.join(self.workdir, self.sol["instance"])
        self.b.get_file(run_key(run_id, "instance.dat"), self.inst)
        self.params = os.path.join(self.workdir, self.sol["params"])
        self.b.get_file(run_key(run_id, "params.par"), self.params)
        self.idle_since = time.time()
        self.busy = 0
        self.units_run = 0
        # the host, for runs that land on different machines
        model = "?"
        try:
            with open("/proc/cpuinfo") as f:
                info = f.read()
            model = next((l.split(":", 1)[1].strip() for l in info.splitlines() if l.startswith("model name")), "?")
            model += " amx" if "amx_bf16" in info else ""
        except OSError:
            pass
        self.host = model
        self.b.put_bytes(run_key(run_id, "workers", self.name.replace("/", "_") + ".json"),
                         json.dumps({"name": self.name, "host": model, "slots": slots, "started": time.time()}).encode())

    def _refresh_run(self):
        self.run = self.b.get_run(self.run_id)
        return self.run

    # ------------------------------------------------------------ one unit
    def _lease(self, msg):
        uid = msg.body["uid"]
        now = time.time()
        unit = self.b.get_unit(self.run_id, uid)
        if unit is None or unit["status"] not in ("queued", "leased"):
            self.b.delete_message(msg)          # finished or gone: a duplicate
            return None
        lease = self.cfg["lease_s"]
        doc = self.b.update_unit(
            self.run_id, uid,
            {"status": "leased", "attempt": unit["attempt"] + 1, "worker": self.name,
             "leased_at": now, "lease_until": now + lease},
            expect={"status": ("queued", "leased"), "attempt": unit["attempt"]},
            check=lambda d: d["status"] == "queued" or d.get("lease_until", 0) < now)
        if doc is None:
            self.b.delete_message(msg)          # someone else holds it: a duplicate
            return None
        return doc

    def _run_unit(self, slot, msg, unit):
        uid, attempt = unit["uid"], unit["attempt"]
        wd = os.path.join(self.workdir, f"slot{slot}", f"{uid}.a{attempt}")
        os.makedirs(wd, exist_ok=True)
        inp, out = os.path.join(wd, "in.state"), os.path.join(wd, "out.state")
        incf, logf, tracef = os.path.join(wd, "inc.txt"), os.path.join(wd, "solver.log"), os.path.join(wd, "trace.txt")
        # resume from the unit's last uploaded checkpoint when an earlier
        # attempt left one (every checkpoint is a valid continuation of
        # in.state, statistics included), else from in.state
        ckpt_key = unit_key(self.run_id, uid, "ckpt.state")
        resumed = self.b.get_file(ckpt_key, inp)
        if not resumed and not self.b.get_file(unit["in_key"], inp):
            raise RuntimeError(f"unit {uid}: no {unit['in_key']}")
        inc_data, inc_tag = self.b.get_if_changed(run_key(self.run_id, "incumbent.json"), None)
        inc = json.loads(inc_data) if inc_data else {}
        self._write_inc(incf, inc)
        # a unit resumes a state file, which already fixes its nodes'
        # assignments: the run's -F belongs to init_run's frontier only (the
        # solver rejects -F together with -r)
        flags = list(self.sol["flags"])
        while "-F" in flags:
            k = flags.index("-F")
            del flags[k:k + 2]
        cmd = [self.solver, "-p", self.params] + flags
        ub = min(x for x in (self.sol.get("ub"), inc.get("obj")) if x is not None) \
            if (self.sol.get("ub") is not None or inc.get("obj") is not None) else None
        if ub is not None and ub < 1e99:
            cmd += ["-u", repr(ub)]
        cmd += ["-U", incf, "-r", inp, "-o", out, "-t", str(unit["budget"]),
                "-k", str(self.cfg["checkpoint_every"]), self.inst]
        env = dict(os.environ)
        env.update(self.sol.get("env", {}))
        env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"})
        if self.cfg.get("trace"):
            env["QAP_DNN_TRACE"] = tracef
        t0 = time.time()
        with open(logf, "w") as lf:
            # its own session: a signal to the supervisor's process group (a
            # terminal's ^C, `kill -TERM -pgid`) reaches the supervisor only
            proc = subprocess.Popen(cmd, cwd=os.path.join(ROOT, "latest"), stdout=lf, stderr=subprocess.STDOUT,
                                    env=env, start_new_session=True)
        asked = False
        killed = False
        drain_at = None
        hb = self.cfg["heartbeat_s"]
        last_hb = time.time()
        last_ckpt, ckpt_mtime = time.time(), None
        while True:
            try:
                rc = proc.wait(timeout=min(1.0, hb))
                break
            except subprocess.TimeoutExpired:
                pass
            now = time.time()
            if self.draining.is_set() and drain_at is None:
                if not asked:
                    proc.send_signal(signal.SIGUSR1)
                    asked = True
                drain_at = now
            if drain_at is not None and now - drain_at > self.grace and not killed:
                proc.kill()
                killed = True
            if now - last_hb < hb:
                continue
            last_hb = now
            # heartbeat: the message and the lease, the incumbent, the control flag
            try:
                self.b.extend(msg, self.cfg["lease_s"])
                if self.b.update_unit(self.run_id, uid, {"lease_until": now + self.cfg["lease_s"]},
                                      expect={"status": "leased", "attempt": attempt}) is None:
                    proc.kill()                 # the lease is gone: someone else runs it
                    killed = True
                elif now - last_ckpt >= self.cfg.get("ckpt_upload_s", 300) and os.path.exists(out):
                    # the solver's -k checkpoint, for a later attempt to resume from
                    m = os.stat(out).st_mtime_ns
                    if m != ckpt_mtime:
                        snap = out + ".up"
                        shutil.copyfile(out, snap)
                        self.b.put_file(ckpt_key, snap)
                        ckpt_mtime = m
                    last_ckpt = now
                data, inc_tag = self.b.get_if_changed(run_key(self.run_id, "incumbent.json"), inc_tag)
                if data:
                    new = json.loads(data)
                    if new.get("obj", float("inf")) < inc.get("obj", float("inf")):
                        inc = new
                        self._write_inc(incf, inc)
                        proc.send_signal(signal.SIGUSR2)
                if not asked:
                    ctl = self.b.get_bytes(unit_key(self.run_id, uid, "control.json"))
                    if ctl and json.loads(ctl).get("return_now") and json.loads(ctl).get("attempt") == attempt:
                        proc.send_signal(signal.SIGUSR1)
                        asked = True
            except Exception as e:      # a backend hiccup must not kill the unit
                self.log(f"worker {self.name} slot {slot}: heartbeat error {e!r}")
        wall = time.time() - t0
        have_out = os.path.exists(out)
        if self.draining.is_set() and rc not in (0, 3):
            killed = True                # died while draining (a SIGTERM to the whole tree): use the checkpoint
        if rc == 0:
            outcome = "done"
        elif rc == 3 or (killed and self.draining.is_set() and have_out):
            outcome = "returned"         # exit 3, or the last checkpoint of a killed solver
        elif killed and not self.draining.is_set():
            outcome = None               # lost the lease: discard
        else:
            outcome = "failed"
        if outcome in ("done", "returned") and not have_out:
            outcome = "failed"
        if killed and self.draining.is_set() and not have_out:
            # nothing saved: give the unit back untouched
            self.b.update_unit(self.run_id, uid, {"status": "queued", "worker": None},
                               expect={"status": "leased", "attempt": attempt})
            self.b.release(msg)
            return "released"
        prefix = unit_key(self.run_id, uid, f"a{attempt}")
        cpu = None
        if have_out:
            try:
                cpu = State(out).times[0]
            except Exception:
                pass
        result = {"uid": uid, "attempt": attempt, "worker": self.name, "slot": slot, "rc": rc,
                  "outcome": outcome, "wall": wall, "cpu": cpu, "killed": killed, "resumed_from_ckpt": resumed,
                  "host": self.host, "finished_at": time.time()}
        # results before the ledger: the store is the truth
        if have_out and outcome in ("done", "returned"):
            self.b.put_file(prefix + "/out.state", out)
        self.b.put_file(prefix + "/solver.log", logf)
        if os.path.exists(tracef):
            self.b.put_file(prefix + "/trace.txt", tracef)
        self.b.put_bytes(prefix + "/result.json", json.dumps(result).encode())
        if outcome is not None:
            ok = self.b.update_unit(self.run_id, uid, {"status": outcome, "rc": rc, "wall": wall, "cpu": cpu,
                                                       "finished_at": result["finished_at"], "host": self.host,
                                                       "resumed_from_ckpt": resumed},
                                    expect={"status": "leased", "attempt": attempt})
            if ok is None:
                self.log(f"worker {self.name}: unit {uid} attempt {attempt} lost its lease; result discarded")
        self.b.delete_message(msg)
        return outcome

    @staticmethod
    def _write_inc(path, inc):
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            if inc.get("obj") is not None and inc.get("obj") < 1e99:
                f.write(repr(float(inc["obj"])) + (" " + " ".join(map(str, inc["perm"])) if inc.get("perm") else "") + "\n")
        os.replace(tmp, path)

    # ------------------------------------------------------------ slots
    def _slot(self, slot):
        poll = self.cfg.get("poll_s", 2.0)
        last_run_check = 0
        while not self.draining.is_set():
            now = time.time()
            if now - last_run_check > 10:
                last_run_check = now
                if self._refresh_run().get("status") == "finished":
                    return
            if self.run.get("status") == "paused":
                time.sleep(poll)
                continue
            msg = self.b.receive(list(LANES), self.cfg["lease_s"])
            if msg is None:
                with self.lock:
                    if self.busy == 0 and self.exit_when_idle and time.time() - self.idle_since > self.exit_when_idle:
                        return
                time.sleep(poll)
                continue
            if msg.body.get("run", self.run_id) != self.run_id:
                # another run's unit: deleting the message would strand it (a
                # coordinator sends a queued unit once)
                self.b.extend(msg, FOREIGN_HIDE_S)
                time.sleep(poll)
                continue
            unit = self._lease(msg)
            if unit is None:
                continue
            with self.lock:
                self.busy += 1
            try:
                outcome = self._run_unit(slot, msg, unit)
                self.log(f"worker {self.name} slot {slot}: unit {unit['uid']} attempt {unit['attempt']} -> {outcome}")
            except Exception as e:
                self.log(f"worker {self.name} slot {slot}: unit {unit['uid']} error {e!r}")
                self.b.update_unit(self.run_id, unit["uid"], {"status": "failed", "rc": -1, "error": repr(e)},
                                   expect={"status": "leased", "attempt": unit["attempt"]})
                try:
                    self.b.delete_message(msg)
                except Exception:
                    pass
            finally:
                with self.lock:
                    self.busy -= 1
                    self.units_run += 1
                    self.idle_since = time.time()

    def drain(self, *_):
        self.draining.set()

    def run_forever(self, install_signals=True):
        if install_signals:
            signal.signal(signal.SIGTERM, self.drain)
            signal.signal(signal.SIGINT, self.drain)
        threads = [threading.Thread(target=self._slot, args=(k,), daemon=True) for k in range(self.slots)]
        for t in threads:
            t.start()
        for t in threads:
            while t.is_alive():
                t.join(timeout=1.0)
        self.log(f"worker {self.name}: {self.units_run} units, exiting")
