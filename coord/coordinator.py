"""The coordinator (docs/distributed/plan.md 2.3-2.5, 4.1, 4.4, 5.3): one
process per run, holding no state that is not in the ledger or the store, so
it can be killed and restarted at any time.

Each tick:
  1. merge finished units exactly once, in merge-sequence order: a unit is
     first marked merged with the next sequence number (a conditional update),
     then its result is added to the run totals and the open nodes it
     returned are split into new units (ids derived from the parent's, so
     re-creating them after a crash is idempotent); the totals are saved with
     the last sequence they include, and a restart replays every unit merged
     after it;
  2. re-queue failed units (up to max_attempts), re-send units whose message
     may be lost;
  3. verify and publish a better incumbent found by a unit (the objective
     recomputed from the original data, the permutation checked);
  4. when units are waiting for nothing to do (idle slots, empty queue), ask
     the longest-running unit to return its open nodes (control.json, read by
     the worker's heartbeat, which sends the solver SIGUSR1), the
     equivalent of parallel.py's steal, but without a global barrier;
  5. write status.json; every publish_every seconds the public status
     (publish.py); when no unit is pending, the final report.
"""
import json
import os
import tempfile
import time

from .common import (LANES, SOLVER, State, fmt_secs, new_unit, plan_split, read_instance, run_key,
                     solver_tables, unit_key, verify_incumbent, write_unit_state)
from . import warm
from .backend import check_queue_name
from .publish import Publisher

PENDING = ("queued", "leased", "done", "returned", "failed")


class Coordinator:
    def __init__(self, backend, run_id, log=print, replay=True):
        self.b = backend
        self.run_id = run_id
        self.log = log
        self.replay = replay   # False: the saved totals as they are (qapctl stop: no unit files fetched)
        self.run = self.b.get_run(run_id)
        if self.run is None:
            raise SystemExit(f"no run {run_id}")
        self.cfg = self.run["config"]
        self.tmp = tempfile.mkdtemp(prefix=f"qapcoord-{run_id}-")
        inst_text = self.b.get_bytes(run_key(run_id, "instance.dat")).decode()
        self.inst = read_instance(inst_text)
        self.spent = []   # units merged since the totals were last saved: their state files go after the save
        self._load_totals()
        self.last_status = 0
        self.last_warm_gc = time.time()
        self.t0 = self.run.get("started_at") or time.time()
        self.publisher = Publisher(self)

    def _path(self, name):
        return os.path.join(self.tmp, name)

    # ------------------------------------------------------------ totals
    def _load_totals(self):
        # totals.json names the totals file that matches its last_seq (written
        # last, in one put, so the pair is always consistent; runs created
        # before Nov 5 have no totals_key and use totals.state)
        meta = json.loads(self.b.get_bytes(run_key(self.run_id, "totals.json")) or b"{}")
        self.totals_key = meta.get("totals_key", run_key(self.run_id, "totals.state"))
        p = self._path("totals.state")
        if not self.b.get_file(self.totals_key, p):
            raise SystemExit(f"run has no {self.totals_key}")
        self.totals = State(p)
        # the header new units are written with; before the replay, which may
        # create units (a restart once crashed here and left the workers idle)
        self.template = State(p)
        self.template.zero_stats()
        self.last_seq = meta.get("last_seq", 0)
        self.busy_wall = meta.get("busy_wall", 0.0)
        self.units_merged = meta.get("units_merged", 0)
        self.dirty = False
        # replay units merged after the saved totals (a crash between marking a
        # unit merged and saving the totals): their stats are not in the file
        replay = self.b.units_merged_after(self.run_id, self.last_seq) if self.replay else []
        for u in sorted(replay, key=lambda d: d["merged_seq"]):
            self.log(f"coordinator: replaying unit {u['uid']} (merge {u['merged_seq']})")
            self._apply(u)

    def _save_totals(self):
        """the totals go to a new file named by the last merge they include, then
        totals.json points at it: one put switches both, so a coordinator killed
        in between restarts from the old pair and replays the merges after it
        (two separate puts once let a kill double-count three units)"""
        if not self.dirty:
            return
        p = self._path("totals.state.out")
        self.totals.write(p, [])
        key = run_key(self.run_id, f"totals-{self.last_seq}.state")
        self.b.put_file(key, p)
        self._crash_point("saving")
        self.b.put_bytes(run_key(self.run_id, "totals.json"), json.dumps(
            {"last_seq": self.last_seq, "totals_key": key, "busy_wall": self.busy_wall,
             "units_merged": self.units_merged}).encode())
        if self.totals_key != key:
            try:
                self.b.delete(self.totals_key)
            except Exception:
                pass                      # an orphan file is harmless
            self.totals_key = key
        self.dirty = False
        self._drop_spent()

    def _drop_spent(self):
        """a merged unit's state files (its input, checkpoint and output) are
        in the saved totals and its open nodes in the units made from them:
        delete them (cfg keep_unit_states keeps them). Results, solver logs and
        node traces stay."""
        spent, self.spent = self.spent, []
        if self.cfg.get("keep_unit_states"):
            return
        for u in spent:
            for key in (u.get("in_key"), unit_key(self.run_id, u["uid"], "ckpt.state"),
                        unit_key(self.run_id, u["uid"], "ckpt.nodes.txt"),
                        unit_key(self.run_id, u["uid"], f"a{u['attempt']}", "out.state")):
                if key:
                    try:
                        self.b.delete(key)
                    except Exception:
                        pass              # left for purge-states

    def _warm_gc(self, force=False):
        """every warm_gc_every seconds (default 600): the parent states nothing can need go"""
        every = self.cfg.get("warm_gc_every", 600)
        if not force and time.time() - self.last_warm_gc < every:
            return
        self.last_warm_gc = time.time()
        kept, deleted = warm.gc(self.b, self.run_id, self.b.all_units(self.run_id))
        if deleted:
            self.log(f"coordinator: parent states: {deleted} deleted, {kept} kept")

    def _crash_point(self, where):
        """tests only (tests/coord_test.py, crash_*): QAP_COORD_CRASH=WHERE:N exits
        the coordinator at once the Nth time it reaches WHERE: marked (a unit
        marked merged, nothing applied), applied (applied, totals not saved),
        saving (the totals file written, totals.json not)"""
        spec = os.environ.get("QAP_COORD_CRASH")
        if not spec or spec.split(":")[0] != where:
            return
        self._crash_hits = getattr(self, "_crash_hits", 0) + 1
        if self._crash_hits >= int(spec.split(":")[1]):
            self.log(f"coordinator: QAP_COORD_CRASH at {where}")
            os._exit(9)

    # ------------------------------------------------------------ merging
    def merge_ready(self):
        n = 0
        for status in ("done", "returned"):
            for u in self.b.units_by_status(self.run_id, status):
                seq = self.b.next_seq(self.run_id)
                doc = self.b.update_unit(self.run_id, u["uid"], {"status": "merged", "outcome": status,
                                                                  "merged_seq": seq},
                                         expect={"status": status, "attempt": u["attempt"]})
                if doc is None:
                    continue           # changed under us; next tick
                self._crash_point("marked")
                self._apply(doc)
                self._crash_point("applied")
                n += 1
        if n:
            self._save_totals()
        return n

    def _apply(self, u):
        """add a merged unit's result to the totals and queue its open nodes"""
        if self.cfg.get("kind") == "dives":
            # a dives unit (coord/dives.py) adds no nodes and returns none; its state file is the
            # source run's, which must keep it (never in self.spent); `qapctl estimate` reads its dives
            self.busy_wall += u.get("wall", 0.0)
            self.units_merged += 1
            self.last_seq = max(self.last_seq, u["merged_seq"])
            self.dirty = True
            return
        p = self._path("out.state")
        key = unit_key(self.run_id, u["uid"], f"a{u['attempt']}", "out.state")
        if not self.b.get_file(key, p):
            raise RuntimeError(f"unit {u['uid']} merged without {key}")
        st = warm.rebase(State(p), self.tmp)
        self.totals.add_stats(st)
        self.busy_wall += u.get("wall", 0.0)
        self.units_merged += 1
        self._candidate(st, u["uid"])
        if u.get("outcome") == "returned" and st.nodes:
            self._make_units(st.nodes, parent=u["uid"])
        self.last_seq = max(self.last_seq, u["merged_seq"])
        self.spent.append(u)
        self.dirty = True

    def _make_units(self, nodes, parent):
        """split nodes into units parent.0, parent.1, ... The split is recorded
        first (split.json) and reused on a replay, and an existing unit is
        never rewritten, so a restart re-creates exactly the same units."""
        plan_key = unit_key(self.run_id, parent or "root", "split.json")
        plan = self.b.get_bytes(plan_key)
        if plan is not None:
            sizes = json.loads(plan)["sizes"]
        else:
            queued = self.b.count_by_status(self.run_id).get("queued", 0)
            sizes = [len(c) for c in plan_split(nodes, self.totals, self.cfg, queued)]
            self.b.put_bytes(plan_key, json.dumps({"sizes": sizes}).encode())
        assert sum(sizes) == len(nodes), (parent, sizes, len(nodes))
        pos = 0
        for k, size in enumerate(sizes):
            chunk = nodes[pos:pos + size]
            pos += size
            uid = f"{parent}.{k}" if parent else f"u{k}"
            have = self.b.get_unit(self.run_id, uid)
            if have is None:
                doc = new_unit(self.run_id, uid, parent, chunk, self.cfg)
                doc["warm"] = warm.node_refs(chunk)   # the parent states it starts from (coord/warm.py)
                path = self._path("in.state")
                self.template.obj, self.template.perm, self.template.found = \
                    self.totals.obj, self.totals.perm, self.totals.found
                write_unit_state(self.template, chunk, path)
                self.b.put_file(doc["in_key"], path)
                self.b.put_unit(self.run_id, doc)
                have = doc
            if not have.get("sent_at"):
                self.b.send(self.run_id, have["lane"], {"run": self.run_id, "uid": uid})
                self.b.update_unit(self.run_id, uid, {"sent_at": time.time()})
        return len(sizes)

    # ------------------------------------------------------------ incumbent
    def _candidate(self, st, uid):
        inc = self.run.get("incumbent", {})
        if st.obj >= inc.get("obj", float("inf")):
            return
        ok, msg = verify_incumbent(self.inst, st.obj, st.perm)
        if not ok:
            self.log(f"coordinator: unit {uid} reports {st.obj} but {msg}; not published")
            return
        new = {"version": inc.get("version", 0) + 1, "obj": st.obj, "perm": [int(x) for x in st.perm],
               "found_by": uid, "found_at": time.time()}
        self.b.put_bytes(run_key(self.run_id, "incumbent.json"), json.dumps(new).encode())
        self.run = self.b.update_run(self.run_id, {"incumbent": new})
        self.totals.obj, self.totals.perm, self.totals.found = st.obj, st.perm, st.found
        self.log(f"New objective value: {st.obj:.17g} (unit {uid}, verified, version {new['version']})")

    # ------------------------------------------------------------ repair
    def repair(self):
        now = time.time()
        for u in self.b.units_by_status(self.run_id, "failed"):
            if u["attempt"] < self.cfg["max_attempts"]:
                if self.b.update_unit(self.run_id, u["uid"], {"status": "queued", "sent_at": now},
                                      expect={"status": "failed", "attempt": u["attempt"]}):
                    self.b.send(self.run_id, u["lane"], {"run": self.run_id, "uid": u["uid"]})
                    self.log(f"coordinator: unit {u['uid']} failed (rc {u.get('rc')}), re-queued")
            else:
                self.b.update_unit(self.run_id, u["uid"], {"status": "dead"}, expect={"status": "failed"})
                self.log(f"coordinator: unit {u['uid']} failed {u['attempt']} times; marked dead")
        lease = self.cfg["lease_s"]
        for u in self.b.units_by_status(self.run_id, "queued"):
            if not u.get("sent_at"):        # created, crash before the send
                self.b.send(self.run_id, u["lane"], {"run": self.run_id, "uid": u["uid"]})
                self.b.update_unit(self.run_id, u["uid"], {"sent_at": now}, expect={"status": "queued"})
        for u in self.b.units_by_status(self.run_id, "leased"):
            # a lease that ran out long ago: its message should have come back;
            # send another in case it was lost (duplicates are harmless)
            if u.get("lease_until", now) < now - 2 * lease and u.get("resent_for") != u["attempt"]:
                self.b.send(self.run_id, u["lane"], {"run": self.run_id, "uid": u["uid"]})
                self.b.update_unit(self.run_id, u["uid"], {"resent_for": u["attempt"]},
                                   expect={"status": "leased", "attempt": u["attempt"]})
                self.log(f"coordinator: unit {u['uid']} lease expired at attempt {u['attempt']}; re-sent")

    # ------------------------------------------------------------ load balancing
    def balance(self, counts):
        """steal: idle slots and an empty queue -> ask the longest-running
        unit(s) to return their open nodes"""
        if self.run.get("status") == "paused":
            for u in self.b.units_by_status(self.run_id, "leased"):
                self._ask_return(u)
            return
        leased = self.b.units_by_status(self.run_id, "leased")
        idle = self.cfg["slots"] - len(leased)
        queued = sum(self.b.depth(self.run_id, l) for l in LANES)
        if idle <= 0 or queued > 0:
            return
        asked = sum(1 for u in leased if u.get("return_requested"))
        now = time.time()
        cand = sorted((u for u in leased if not u.get("return_requested")
                       and now - u.get("leased_at", now) >= self.cfg["steal_after"]),
                      key=lambda u: u.get("leased_at", now))
        for u in cand[:max(0, idle - asked)]:
            self._ask_return(u)

    def requeue_leased(self):
        """after a crash of every worker (a machine restart): give their units
        back at once instead of waiting for the leases to run out"""
        n = 0
        for u in self.b.units_by_status(self.run_id, "leased"):
            if self.b.update_unit(self.run_id, u["uid"], {"status": "queued", "worker": None, "sent_at": time.time()},
                                  expect={"status": "leased", "attempt": u["attempt"]}):
                self.b.send(self.run_id, u["lane"], {"run": self.run_id, "uid": u["uid"]})
                n += 1
        self.log(f"coordinator: {n} leased units re-queued")
        return n

    def _ask_return(self, u):
        if u.get("return_requested") == u["attempt"]:
            return
        self.b.put_bytes(unit_key(self.run_id, u["uid"], "control.json"),
                         json.dumps({"return_now": True, "attempt": u["attempt"]}).encode())
        self.b.update_unit(self.run_id, u["uid"], {"return_requested": u["attempt"]},
                           expect={"status": "leased", "attempt": u["attempt"]})

    # ------------------------------------------------------------ status
    def status(self, counts, force=False):
        now = time.time()
        if not force and now - self.last_status < self.cfg.get("status_every", 30):
            return
        self.last_status = now
        t = self.totals
        dc = t.arrays.get("dnn_counts", [0, 0])
        wall = now - self.t0
        doc = {"run": self.run_id, "status": self.run.get("status"), "updated": now, "wall_s": wall,
               "units": counts, "queue": {l: self.b.depth(self.run_id, l) for l in LANES},
               "nodes": t.counters[0], "fathomed": t.counters[1], "dnn_nodes": dc[0], "dnn_fathomed": dc[1],
               "cpu_s": t.times[0], "busy_wall_s": self.busy_wall, "units_merged": self.units_merged,
               "busy_fraction_merged": self.busy_wall / (self.cfg["slots"] * wall) if wall > 0 else 0,
               "incumbent": self.run.get("incumbent", {}).get("obj")}
        self.b.put_bytes(run_key(self.run_id, "status.json"), json.dumps(doc, indent=1).encode())
        self.log(f"coordinator: {fmt_secs(wall)} units {counts} nodes {t.counters[0]} "
                 f"(DNN {dc[0]}, fathomed {dc[1]}) CPU {t.times[0]:.0f} s")

    def _solver_tables(self, state_path):
        """the solver's per-level tables for the merged statistics (the
        solver reprints a state file with no open nodes); "" without a solver"""
        sol = self.run.get("solver", {})
        params = self._path(sol.get("params") or "params.par")
        inst = self._path(sol.get("instance") or "instance.dat")
        for key, path in ((run_key(self.run_id, "params.par"), params), (run_key(self.run_id, "instance.dat"), inst)):
            if not os.path.exists(path) and not self.b.get_file(key, path):
                return ""
        return solver_tables(SOLVER, params, state_path, inst)

    def report(self, state_path, wall):
        """the solver's RESULTS block for the merged statistics written at
        state_path: the head, the solver's own tables, the coordinator's line"""
        t = self.totals
        inc = self.run.get("incumbent", {})
        node, fathom, enum, fw = t.counters
        lines = ["  ---- RESULTS ----", "",
                 f"    Total Time: {t.times[0]:.6f} seconds ",
                 f"         Nodes: {node} / {fathom} / {enum}",
                 f"    root bound: {t.root_bound}",
                 f" fw iterations: {fw}"]
        dc = t.arrays.get("dnn_counts", [0, 0])
        if dc[0]:
            lines.append(f"     DNN bound: {dc[0]} nodes, {dc[1]} fathomed, {t.arrays['dnn_time'][0]:.2f} s")
        if sum(t.arrays.get("cut_n_level", [])):   # state version 7
            a = t.arrays
            lines.append(f"     cut phase: {sum(a['cut_n_level'])} nodes, {sum(a['cut_f_level'])} fathomed, "
                         f"{sum(a['cut_t_level']):.2f} s")
        tables = self._solver_tables(state_path)
        if tables:
            lines += [""] + tables.rstrip("\n").split("\n")
        busy = self.busy_wall / (self.cfg["slots"] * wall) if wall > 0 else 0
        lines += ["", f"coordinator: {self.cfg['slots']} slots, {self.units_merged} units, wall {wall:.2f} s, "
                      f"CPU {t.times[0]:.2f} s (x{t.times[0] / wall if wall else 0:.2f}), slots busy {100 * busy:.1f}%, "
                      f"best objective {inc.get('obj', t.obj):.17g}"]
        return "\n".join(lines)

    def finish(self):
        now = time.time()
        wall = now - self.t0
        t = self.totals
        inc = self.run.get("incumbent", {})
        if inc.get("perm"):
            t.obj, t.perm = inc["obj"], [str(x) for x in inc["perm"]]
        p = self._path("final.state")
        t.write(p, [])
        self.b.put_file(run_key(self.run_id, "final.state"), p)
        if self.cfg.get("kind") == "dives":
            from .dives import estimate_run, report as dives_report
            report = dives_report(estimate_run(self.b, self.run_id, log=lambda *_: None))
        else:
            report = self.report(p, wall)
        self.b.put_bytes(run_key(self.run_id, "results.txt"), (report + "\n").encode())
        self.run = self.b.update_run(self.run_id, {"status": "finished", "finished_at": now, "wall_s": wall})
        try:
            drop_queues(self.b, self.run_id)   # its workers see it finished; nothing is sent to it again
        except Exception as e:
            self.log(f"coordinator: deleting the run's queues failed: {e!r}")
        # the run is over: no node will start from a parent state again
        for key in self.b.list(run_key(self.run_id, "warm") + "/"):
            try:
                self.b.delete(key)
            except Exception:
                pass
        # the last status.json says finished, with the final totals (a run
        # shorter than status_every would otherwise leave its first one)
        counts = self.b.count_by_status(self.run_id)
        self.status(counts, force=True)
        try:
            self.publisher.publish(counts, final=True)
        except Exception as e:
            self.log(f"coordinator: publishing the final status failed: {e!r}")
        self.log(report)
        return report

    # ------------------------------------------------------------ loop
    def tick(self):
        self.run = self.b.get_run(self.run_id)
        if self.run.get("status") == "stopped":
            return self.b.count_by_status(self.run_id)   # loop returns at its next look
        merged = self.merge_ready()
        self.repair()
        counts = self.b.count_by_status(self.run_id)
        self.balance(counts)
        self.status(counts, force=merged > 0 and self.cfg.get("status_on_merge", False))
        try:
            self._warm_gc()
        except Exception as e:          # the store must never stop a run
            self.log(f"coordinator: the parent-state cleanup failed: {e!r}")
        if self.publisher.due():
            try:
                self.publisher.publish(counts)
            except Exception as e:      # the public status must never stop a run
                self.log(f"coordinator: publishing the status failed: {e!r}")
        return counts

    def loop(self, poll=2.0, max_wall=None):
        if not self.run.get("started_at"):
            self.run = self.b.update_run(self.run_id, {"started_at": self.t0})
        if self.run.get("status") in ("running", "paused"):
            ensure_queues(self.b, self.run_id, log=self.log)
            self.run = self.b.get_run(self.run_id)
        while True:
            if self.run.get("status") == "stopped":
                self.log(f"coordinator: run {self.run_id} was stopped (qapctl stop)")
                return None
            if self.run.get("status") == "finished":
                return self.b.get_bytes(run_key(self.run_id, "results.txt")).decode()
            counts = self.tick()
            if not any(counts.get(s, 0) for s in PENDING):
                self._save_totals()
                return self.finish()
            if counts.get("dead"):
                self._save_totals()
                raise SystemExit(f"coordinator: {counts['dead']} unit(s) dead; see the ledger")
            if max_wall and time.time() - self.t0 > max_wall:
                self._save_totals()
                return None
            time.sleep(poll)


def ensure_queues(b, run_id, log=print, resend=False):
    """the run's queues, created if missing (idempotent). A run from before per-run queues (its
    record has no "queues": its messages are in the old shared lanes, which no worker reads now),
    a run whose queues were deleted, or resend: each queued unit is sent again. A unit sent twice
    is harmless: the second copy finds it taken and is deleted (coord/worker.py, _lease)."""
    b.create_queues(run_id)
    run = b.get_run(run_id)
    if run is None or (run.get("queues") == "run" and not resend):
        return 0
    n = 0
    for u in b.units_by_status(run_id, "queued"):
        b.send(run_id, u["lane"], {"run": run_id, "uid": u["uid"]})
        n += 1
    b.update_run(run_id, {"queues": "run"})
    if n:
        log(f"coordinator: run {run_id}: {n} queued units sent to its queues")
    return n


def drop_queues(b, run_id):
    """delete a run's queues (finished or stopped: nothing will run again); the record says so,
    so that a later ensure_queues sends its queued units again"""
    b.delete_queues(run_id)
    if b.get_run(run_id) is not None:
        b.update_run(run_id, {"queues": None})


def cache_key(key, name=""):
    """the root cache's prefix: under runs/, which every role may read and write, as .cache (a run
    name has no dot)"""
    return f"runs/.cache/frontier/{key}" + (f"/{name}" if name else "")


def frontier_key(solver, sdir, instance, params, flags, env, ub, depth, trace):
    """what a frontier pass's output depends on: the solver binary, the CPU (OpenBLAS picks its kernels
    by it, and the rounding with them), the data, the files the flags name (a symmetry file), the
    flags, the incumbent, the depth, the traces, and the solver's resolved settings under the params
    file and the environment (bnbdrv --print-settings: every setting but the run class, with where
    it came from)"""
    import hashlib
    import subprocess
    h = hashlib.sha256()

    def file(path):
        with open(path, "rb") as f:
            h.update(hashlib.sha256(f.read()).digest())
    file(solver)
    file(os.path.abspath(instance))
    for tok in flags:
        p = tok if os.path.isabs(tok) else os.path.join(sdir, tok)
        if os.path.isfile(p):
            file(p)
    cpu = ""
    try:
        with open("/proc/cpuinfo") as f:
            cpu = next((l.split(":", 1)[1].strip() for l in f if l.startswith("model name")), "")
    except OSError:
        pass
    r = subprocess.run([solver, "--print-settings", "-p", params] + flags, cwd=sdir, capture_output=True, text=True, env=env)
    # the run class (paths, traces, the estimator's) does not change the tree, as for the settings hash
    settings = [l for l in r.stdout.splitlines() if l.startswith("settings ") or
                (l.startswith("setting ") and not l.endswith(" run") and not l.startswith("setting par.file "))]
    h.update(json.dumps([cpu, flags, repr(ub), depth, trace, settings]).encode())
    return h.hexdigest()[:32]


def frontier_to_cache(b, key, tmp, traces, stdout):
    """the frontier pass's state file, the parent states it names, its traces and its report"""
    frontier = os.path.join(tmp, "frontier.state")
    for name in warm.refs_in(frontier):
        b.put_file(cache_key(key, f"warm/{name}"), os.path.join(tmp, "warm", name))
    for name in traces:
        if os.path.exists(os.path.join(tmp, name)):
            b.put_file(cache_key(key, name), os.path.join(tmp, name))
    b.put_bytes(cache_key(key, "stdout.txt"), stdout.encode())
    b.put_file(cache_key(key, "frontier.state"), frontier)   # last: its presence says the rest is there


def frontier_from_cache(b, key, tmp, traces):
    """True when the cache had it: tmp then holds what the frontier pass would have written"""
    frontier = os.path.join(tmp, "frontier.state")
    if not b.get_file(cache_key(key, "frontier.state"), frontier):
        return False
    os.makedirs(os.path.join(tmp, "warm"), exist_ok=True)
    for name in warm.refs_in(frontier):
        if not b.get_file(cache_key(key, f"warm/{name}"), os.path.join(tmp, "warm", name)):
            os.remove(frontier)
            return False
    for name in traces:
        b.get_file(cache_key(key, name), os.path.join(tmp, name))
    return True


def init_run(backend, run_id, instance, params, flags, env, ub, cfg, solver, frontier_depth=0, log=print,
             visibility="public"):
    """the frontier pass, run here (as bench/parallel.py does): the root, or the
    -F node, expanded to frontier_depth; its open nodes become the first units"""
    import shutil
    import subprocess
    if backend.get_run(run_id):
        raise SystemExit(f"run {run_id} exists")
    check_queue_name(run_id)   # before the frontier pass, which can take hours
    tmp = tempfile.mkdtemp(prefix=f"qapinit-{run_id}-")
    sdir = os.path.dirname(os.path.abspath(solver))
    frontier = os.path.join(tmp, "frontier.state")
    cmd = [solver, "-p", params] + flags + (["-u", repr(ub)] if ub is not None else []) + \
          ["-d", str(frontier_depth), "-o", frontier, os.path.abspath(instance)]
    e = dict(os.environ)
    e.update(env)
    traces = {"QAP_NODE_TRACE": "nodes.txt", "QAP_DNN_TRACE": "trace.txt"} if cfg.get("trace") else {}
    for var, name in traces.items():
        e[var] = os.path.join(tmp, name)   # the frontier's nodes are in no unit's trace
    e["QAP_DNN_WARM_DIR"] = os.path.join(tmp, "warm")   # the frontier's parents, once each (coord/warm.py)
    # the root cache: a frontier pass already made on this solver, CPU, data, flags, settings and
    # incumbent is read back instead of being run again (tai50b's root: 72 minutes, tai60b's: hours;
    # a relaunch after an interruption during init, or another estimate of the same problem)
    key = None if os.environ.get("QAP_ROOT_CACHE", "on") == "off" else \
        frontier_key(solver, sdir, instance, params, flags, e, ub, frontier_depth, bool(traces))
    if key and frontier_from_cache(backend, key, tmp, list(traces.values())):
        stdout = (backend.get_bytes(cache_key(key, "stdout.txt")) or b"").decode()
        log(f"coordinator: the frontier pass from the root cache ({cache_key(key)}), not run again")
    else:
        r = subprocess.run(cmd, cwd=sdir, capture_output=True, text=True, env=e)
        if r.returncode != 0:
            raise SystemExit(f"frontier pass failed (rc={r.returncode}):\n{r.stdout[-2000:]}{r.stderr}")
        stdout = r.stdout
        if key:
            try:
                frontier_to_cache(backend, key, tmp, list(traces.values()), stdout)
            except Exception as ex:      # the cache must never stop a run
                log(f"coordinator: writing the root cache failed: {ex!r}")
    log(stdout.split("  ---- RESULTS ----")[0].rstrip())
    for name in traces.values():
        if os.path.exists(os.path.join(tmp, name)):
            backend.put_file(run_key(run_id, "frontier", name), os.path.join(tmp, name))
    warm.push(backend, run_id, frontier, set(), log=log)
    st = State(frontier)
    backend.put_file(run_key(run_id, "instance.dat"), os.path.abspath(instance))
    params_path = params if os.path.isabs(params) else os.path.join(sdir, params)
    backend.put_file(run_key(run_id, "params.par"), params_path)
    st.write(os.path.join(tmp, "totals.state"), [])
    backend.put_file(run_key(run_id, "totals.state"), os.path.join(tmp, "totals.state"))
    backend.put_bytes(run_key(run_id, "totals.json"), json.dumps({"last_seq": 0}).encode())
    inc = {"version": 0, "obj": st.obj}
    if st.perm and all(x != "-1" for x in st.perm):
        inc["perm"] = [int(x) for x in st.perm]
    backend.put_bytes(run_key(run_id, "incumbent.json"), json.dumps(inc).encode())
    backend.create_queues(run_id)
    backend.put_run(run_id, {"run_id": run_id, "status": "running", "created": time.time(), "incumbent": inc,
                             "visibility": visibility, "queues": "run",
                             "config": cfg, "solver": {"flags": flags, "env": env, "ub": ub,
                                                       "params": os.path.basename(params),
                                                       "instance": os.path.basename(instance)}})
    c = Coordinator(backend, run_id, log=log)
    n = c._make_units(warm.rebase(st, c.tmp).nodes, parent=None)
    log(f"coordinator: run {run_id}: frontier at depth {frontier_depth} has {len(st.nodes)} open nodes, "
        f"{n} units, {cfg['slots']} slots")
    shutil.rmtree(tmp)
    return c
