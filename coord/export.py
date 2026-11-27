"""qapctl export: a --trace run's whole traces, out of the store, for analysis
(docs/tho40b_results_plan.md).

  nodes.txt.gz   every node once, in merge order: the QAP_NODE_TRACE lines of
                 the frontier pass and of each merged unit's recorded attempt
                 (what qapctl nodes writes; format in latest/bnb.C)
  dnn.txt.gz     every DNN solve of those nodes once: the QAP_DNN_TRACE line
                 (iterations, give-up or fathom, warm start, the cut phase) of
                 each node in nodes.txt.gz whose DNN bound ran, in the same order
  manifest.json  what was read and what is missing

A unit's node trace is exact (coord/worker.py: a resumed attempt starts from
its checkpoint's lines, a killed one is cut to its checkpoint). Its DNN trace
is carried across checkpoints too (ckpt.trace.txt) but not cut, so it may
hold a node twice: solved before a kill and again after the checkpoint the
next attempt resumed. The DNN lines are matched to the nodes by their fixes
(the DNN trace lists every fixed facility in facility order, the node trace
the fixes in branching order: the same set), each node taking the last line
for its fixes in the recorded attempt's trace.txt, else in the unit's
ckpt.trace.txt or an earlier attempt's (runs from before Nov 15 did not carry
the DNN trace: a SIGKILLed worker's lines before its checkpoint were lost). A
node of the frontier pass takes the frontier's line.

The manifest names no backend: an S3 bucket's name holds the AWS account
number, and the collect commits the manifest's numbers to the repository.
"""
import gzip
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor

from .common import run_key, unit_key


def node_fixes(line):
    """the fixes of a QAP_NODE_TRACE line as a sorted tuple of (facility, location)"""
    f = line.split(None, 3)
    if len(f) < 3 or f[0] != "node" or f[2] == "-":
        return ()
    return tuple(sorted(tuple(int(x) for x in pair.split(":")) for pair in f[2].split(",")))


def node_ran_dnn(line):
    """whether a QAP_NODE_TRACE line records a DNN bound (its dnn field is not '-')"""
    f = line.split(None, 5)
    return len(f) > 4 and f[4] != "-"


def dnn_fixes(line):
    """the fixes of a QAP_DNN_TRACE line ('pid P depth D fixes a:b c:d | ...')"""
    head = line.split("|", 1)[0].split()
    try:
        k = head.index("fixes")
    except ValueError:
        return None
    return tuple(sorted(tuple(int(x) for x in pair.split(":")) for pair in head[k + 1:]))


def _lines(data):
    return [ln for ln in data.decode(errors="replace").splitlines() if ln.strip()] if data else []


def export_run(b, run_id, out_dir, log=print, threads=64):
    """write nodes.txt.gz, dnn.txt.gz and manifest.json under out_dir; 0 when
    both traces are whole, 2 when some unit's trace is missing, 1 when the run
    has no trace.

    The reads wait on the round trip, not on the bytes: tai50b-261116's export (7,151 units,
    8,225 nodes, almost no data) took 396 s on 16 threads, 55 ms a unit. So 64 threads (16 until
    Nov 24); the S3 client keeps as many connections (coord/backend.py). The rest is one thread's
    work, about 28 us a node (2/3 of it gzip), 80 s for tho40-261122's 2.7 million nodes. The
    output does not depend on the count: map returns the units in merge order."""
    run = b.get_run(run_id)
    if run is None or not run["config"].get("trace"):
        log(f"run {run_id}: not initialized with --trace, nothing to export")
        return 1
    os.makedirs(out_dir, exist_ok=True)
    units = sorted(b.all_units(run_id), key=lambda u: (u.get("merged_seq") or 0, u["uid"]))
    merged = [u for u in units if u.get("status") == "merged"]

    def fetch(u):
        rec = u["attempt"]
        nodes = b.get_bytes(unit_key(run_id, u["uid"], f"a{rec}", "nodes.txt"))
        # the recorded attempt's DNN trace first, then the carried one, then the earlier attempts'
        traces = [b.get_bytes(unit_key(run_id, u["uid"], f"a{rec}", "trace.txt"))]
        if rec > 0 or u.get("resumed_from_ckpt"):
            traces.append(b.get_bytes(unit_key(run_id, u["uid"], "ckpt.trace.txt")))
            traces += [b.get_bytes(unit_key(run_id, u["uid"], f"a{a}", "trace.txt")) for a in range(rec - 1, -1, -1)]
        return nodes, traces

    t0 = time.time()
    counts = {"node_lines": 0, "dnn_lines": 0, "dnn_nodes": 0, "dnn_missing": 0,
              "units_merged": len(merged), "units_without_nodes": 0, "units_not_merged": len(units) - len(merged)}
    missing_units = []

    def emit(nodes_lines, traces, nf, df):
        # within a trace the last line for a node's fixes (its counted solve); a trace
        # earlier in the list wins over the ones after it
        by_fixes = {}
        for t in traces:
            last = {}
            for ln in t:
                k = dnn_fixes(ln)
                if k is not None:
                    last[k] = ln
            for k, ln in last.items():
                by_fixes.setdefault(k, ln)
        for ln in nodes_lines:
            nf.write(ln + "\n")
            counts["node_lines"] += 1
            if node_ran_dnn(ln):
                counts["dnn_nodes"] += 1
                d = by_fixes.pop(node_fixes(ln), None)
                if d is None:
                    counts["dnn_missing"] += 1
                else:
                    df.write(d + "\n")
                    counts["dnn_lines"] += 1

    ntmp, dtmp = os.path.join(out_dir, "nodes.txt.gz.tmp"), os.path.join(out_dir, "dnn.txt.gz.tmp")
    with gzip.open(ntmp, "wt") as nf, gzip.open(dtmp, "wt") as df:
        emit(_lines(b.get_bytes(run_key(run_id, "frontier", "nodes.txt"))),
             [_lines(b.get_bytes(run_key(run_id, "frontier", "trace.txt")))], nf, df)
        with ThreadPoolExecutor(threads) as pool:
            # map keeps the merge order; fetching runs ahead on the pool
            for u, (nodes, traces) in zip(merged, pool.map(fetch, merged)):
                if nodes is None:
                    counts["units_without_nodes"] += 1
                    missing_units.append(u["uid"])
                    continue
                emit(_lines(nodes), [_lines(t) for t in traces], nf, df)
    os.replace(ntmp, os.path.join(out_dir, "nodes.txt.gz"))
    os.replace(dtmp, os.path.join(out_dir, "dnn.txt.gz"))
    cfg = run["config"]
    manifest = {"run": run_id, "status": run.get("status"), "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "instance": (run.get("solver") or {}).get("instance"), "flags": (run.get("solver") or {}).get("flags"),
                "trace": bool(cfg.get("trace")), **counts, "units_missing_nodes": missing_units[:50],
                "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    log(f"run {run_id}: {counts['node_lines']} node lines and {counts['dnn_lines']} DNN lines "
        f"({counts['dnn_missing']} of {counts['dnn_nodes']} DNN nodes without one) from the frontier and "
        f"{len(merged) - counts['units_without_nodes']} merged units in {out_dir}; "
        f"{counts['units_not_merged']} units not merged")
    return 2 if (counts["units_without_nodes"] or counts["dnn_missing"]) else 0
