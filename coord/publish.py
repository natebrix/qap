"""The public status of a run (docs/distributed/plan.md section 7.1), which
brixius-web reads. Every publish_every seconds (init --publish-every, default
300) and once at the end, the coordinator writes under public/runs/<run>/ in
the store:

  status.json    what the run is doing: status, units, incumbent, CPU, and
                 the open nodes by depth (as handed out in queued and running
                 units); the names of the other files
  totals.state   the merged statistics (QAP-STATE 7, no open nodes): every
                 per-level table of the solver's report, the DNN ones too
  report.txt     the solver's report for them: the RESULTS head, the
                 solver's own tables, the coordinator's line
  history.json   one sample per publish, for charts over time
  params.par     the params file the run was given (its levels and its
                 [settings]), for the run page's branching rules
  tree/          a --trace run's tree, every merged node (coord/tree.py):
                 meta.json and the shards that changed since the last publish

status.json goes last, so a reader that sees a new one finds the others at
least as new. Then the catalog of public runs, public/runs/index.json
(brixius-web's runs page): one entry per run, upserted at every publish from
a private master copy, runs/index.json, which a run's role may read. Two runs
publishing at once can drop each other's entry; the next publish puts it back,
and the final one checks that its entry is there. Nothing in the bundle names the bucket, an instance or a
credential. A run created with --visibility private publishes nothing. The
store keeps the history privately as well (runs/<run>/history.json), so a
restarted coordinator goes on with it. With the AWS backend the objects carry
a content type and Cache-Control max-age=60 for CloudFront
(infra/phase1/publish_setup.sh serves public/ and nothing else).
"""
import json
import os
import time
from datetime import datetime, timezone

from .common import run_key
from .tree import TreePublisher, grid_of

PUBLIC = "public/runs"
CATALOG = "runs/index.json"     # the master copy; public/runs/index.json is what readers fetch
MAX_HISTORY = 2000          # thinned to every other sample beyond this
CACHE = "max-age=60"
OPEN_STATUSES = ("queued", "leased", "failed")


def public_key(run_id, name):
    return f"{PUBLIC}/{run_id}/{name}"


def iso(t):
    return datetime.fromtimestamp(t, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") if t else None


def open_by_depth(backend, run_id, levels):
    """open nodes by depth: the nodes handed out in units not yet finished
    (a running unit's are counted as handed out, not as it has left them)"""
    out = [0] * (levels + 1)
    for status in OPEN_STATUSES:
        for u in backend.units_by_status(run_id, status):
            by_depth = u.get("by_depth")
            if by_depth:
                for d, count in by_depth.items():
                    out[min(int(d), levels)] += count
            else:                              # units made before by_depth existed
                out[min(u.get("depth_min", 0), levels)] += u.get("nodes", 0)
    return out


class Publisher:
    def __init__(self, coord):
        self.c = coord
        self.last = 0.0
        data = coord.b.get_bytes(run_key(coord.run_id, "history.json"))
        self.history = json.loads(data) if data else []
        # a --trace run's tree (coord/tree.py), published with the status
        self.tree = None
        if coord.cfg.get("trace"):
            n, a, b, _ = coord.inst
            sol = coord.run.get("solver", {})
            commit = os.environ.get("QAP_COMMIT")
            self.tree = TreePublisher(coord.b, coord.run_id, os.path.splitext(sol.get("instance") or "")[0] or None,
                                      n, grid_of(n, b, a),
                                      f"coordinator run {coord.run_id}, QAP_NODE_TRACE" +
                                      (f" (qap {commit[:7]})" if commit else ""), log=coord.log)

    def due(self, now=None):
        every = self.c.cfg.get("publish_every", 300)
        return every > 0 and (now or time.time()) - self.last >= every

    def private(self):
        return self.c.run.get("visibility", "public") == "private"

    def publish(self, counts, final=False):
        """write the bundle; False (and nothing written) for a private run"""
        c, b, run_id = self.c, self.c.b, self.c.run_id
        if self.private():
            return False
        now = time.time()
        self.last = now
        t = c.totals
        inc = c.run.get("incumbent", {})
        if inc.get("perm"):
            t.obj, t.perm = inc["obj"], [str(x) for x in inc["perm"]]
        totals_path = c._path("public-totals.state")
        t.write(totals_path, [])
        # a finished run's wall time is the one finish() recorded; a later
        # publish (qapctl publish, a backfill) must not count the time since
        wall = c.run["wall_s"] if final and c.run.get("wall_s") is not None else now - c.t0
        report = c.report(totals_path, wall)
        node, fathom, enum, fw = t.counters
        open_d = open_by_depth(b, run_id, t.levels)
        dc = t.arrays.get("dnn_counts", [0, 0])
        sol = c.run.get("solver", {})
        obj = inc.get("obj")
        status = {
            "schema": 1,
            "run": run_id,
            "status": c.run.get("status"),
            "final": final,
            "generated_at": iso(now),
            "started_at": iso(c.t0),
            "finished_at": iso(c.run.get("finished_at")) if final else None,
            "wall_s": round(wall, 1),
            "instance": os.path.splitext(sol.get("instance") or "")[0] or None,
            "n": t.n,
            "solver": {"params": sol.get("params"), "flags": sol.get("flags", []), "env": sol.get("env", {}),
                       "ub": sol.get("ub"), "commit": os.environ.get("QAP_COMMIT") or None},
            "incumbent": {"objective": obj if obj is not None and obj < 1e99 else None,
                          "version": inc.get("version", 0), "found_at": iso(inc.get("found_at")),
                          "found_by_search": obj is not None and obj < 1e99 and
                          (sol.get("ub") is None or obj < sol["ub"])},
            "root_bound": float(t.root_bound),
            "nodes": {"evaluated": node, "fathomed": fathom, "enumerated": enum, "open": sum(open_d),
                      "open_by_depth": open_d},
            "fw_iterations": fw,
            "dnn": {"nodes": dc[0], "fathomed": dc[1], "cpu_s": t.arrays.get("dnn_time", [0.0])[0]},
            "cpu_s": t.times[0],
            "units": counts,
            "slots": c.cfg["slots"],
            "busy_fraction": round(c.busy_wall / (c.cfg["slots"] * wall), 4) if wall > 0 else 0,
            "files": {"totals": "totals.state", "report": "report.txt", "history": "history.json"},
        }
        # the params file the solver ran with (its levels and [settings]): the
        # branching rules a run page lists; the file is in the public repo anyway
        params = b.get_bytes(run_key(run_id, "params.par"))
        if params is not None:
            status["files"]["params"] = "params.par"
        if self.tree is not None:
            try:
                if self.tree.publish(final, now):
                    status["files"]["tree"] = "tree/meta.json"
            except Exception as e:      # the tree must never stop the status, nor the run
                c.log(f"coordinator: publishing the tree failed: {e!r}")
        sample = {"t": iso(now), "wall_s": round(wall, 1), "nodes": node, "fathomed": fathom,
                  "open": sum(open_d), "cpu_s": round(t.times[0], 1),
                  "dnn_nodes": dc[0], "dnn_fathomed": dc[1],
                  "incumbent": status["incumbent"]["objective"], "units_merged": c.units_merged}
        if final:
            # a re-publish of a finished run: nothing after the finish (an
            # earlier re-publish counted the time since), and the final sample once
            self.history = [h for h in self.history if h["wall_s"] <= sample["wall_s"]]
            last = self.history[-1] if self.history else None
            if last is not None and (last["wall_s"], last["nodes"]) == (sample["wall_s"], sample["nodes"]):
                self.history.pop()
        self.history.append(sample)
        if len(self.history) > MAX_HISTORY:   # keep the newest half whole, thin the older half
            half = len(self.history) // 2
            self.history = self.history[:half:2] + self.history[half:]
        history = json.dumps(self.history).encode()
        b.put_bytes(run_key(run_id, "history.json"), history)
        with open(totals_path, "rb") as f:
            b.put_bytes(public_key(run_id, "totals.state"), f.read(), "text/plain; charset=utf-8", CACHE)
        b.put_bytes(public_key(run_id, "report.txt"), (report + "\n").encode(), "text/plain; charset=utf-8", CACHE)
        b.put_bytes(public_key(run_id, "history.json"), history, "application/json", CACHE)
        if params is not None:
            b.put_bytes(public_key(run_id, "params.par"), params, "text/plain; charset=utf-8", CACHE)
        b.put_bytes(public_key(run_id, "status.json"), json.dumps(status, indent=1).encode(),
                    "application/json", CACHE)
        update_catalog(b, catalog_entry(status), verify=final)
        return True


def catalog_entry(status):
    """a run's line in the catalog: what a list of runs shows, from status.json"""
    inc = status["incumbent"]
    return {"run": status["run"], "instance": status["instance"], "n": status["n"], "status": status["status"],
            "final": status["final"], "started_at": status["started_at"], "finished_at": status["finished_at"],
            "generated_at": status["generated_at"], "wall_s": status["wall_s"], "cpu_s": round(status["cpu_s"], 1),
            "nodes": status["nodes"]["evaluated"], "fathomed": status["nodes"]["fathomed"],
            "open": status["nodes"]["open"], "incumbent": inc["objective"], "found_by_search": inc["found_by_search"],
            "slots": status["slots"], "tree": "tree" in status["files"]}


def update_catalog(b, entry, verify=False, tries=3):
    """upsert entry into the catalog (newest start first) and publish it; with
    verify, read it back and try again if another run's write dropped it"""
    for _ in range(tries):
        data = b.get_bytes(CATALOG)
        runs = json.loads(data).get("runs", []) if data else []
        runs = [r for r in runs if r.get("run") != entry["run"]] + [entry]
        runs.sort(key=lambda r: (r.get("started_at") or "", r.get("run")), reverse=True)
        body = json.dumps({"schema": 1, "generated_at": entry["generated_at"], "runs": runs}, indent=1).encode()
        b.put_bytes(CATALOG, body)
        b.put_bytes(f"{PUBLIC}/index.json", body, "application/json", CACHE)
        if not verify:
            return
        back = b.get_bytes(CATALOG)
        if back and any(r.get("run") == entry["run"] and r.get("final") == entry["final"]
                        for r in json.loads(back).get("runs", [])):
            return
