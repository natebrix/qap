"""A run's readout: one Markdown page that says where the run stands and what it found, written to
the store as runs/RUN/readout.md (docs/review_nov18_runs.md, items 4 and 5: a readout anyone can read
without a person collecting it by hand).

- The coordinator writes it when the run finishes or is stopped (Coordinator.finish, qapctl stop); a
  dives run's finish writes its source's (an estimate run's readout is its estimate).
- `qapctl readout RUN [--save] [--out FILE]` makes it at any time, part-way included.
- `infra/phase1/collect.sh RUN` fetches it into the run's record, and the GitHub workflow
  phase1-readouts commits the readouts of the runs in the store (infra/phase1/readouts.py).

Its first line after the title is the verdict: COMPLETE, PART-WAY, RUNNING, PAUSED, STOPPED, with
the warnings below it (an estimate's flags, coord/dives.py; a status not written for a while; units
the coordinator gave up on). The rest is the record: what was run, when, the work, the results.
"""
import json
import re
import time

from .common import run_key

STALE_S = 900      # a running run whose coordinator has not written its status for this long


def _when(t):
    return time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime(t)) if t else "-"


def _h(s):
    return f"{(s or 0) / 3600:,.1f}"


def _wall(s):
    s = int(s or 0)
    return f"{s} s" if s < 600 else f"{s // 3600} h {s % 3600 // 60:02d} min"


def _cpu(s):
    s = s or 0.0
    return f"{s / 3600:,.1f} h" if s >= 360 else f"{s:,.1f} s"


def _status(b, run_id):
    try:
        return json.loads(b.get_bytes(run_key(run_id, "status.json")) or b"{}")
    except ValueError:
        return {}


def _results(b, run_id):
    data = b.get_bytes(run_key(run_id, "results.txt"))
    return data.decode() if data else ""


def _pick(text, pattern):
    m = re.search(pattern, text, re.M)
    return m.group(0).strip() if m else None


def _estimate(b, dives_run):
    """the dives run's estimate: its finished results, or computed now part-way"""
    from .dives import estimate_run
    d = b.get_run(dives_run)
    if d is None:
        return None, None
    if d.get("status") == "finished":
        txt = (b.get_bytes(run_key(dives_run, "estimate.txt")) or b"").decode()
        js = b.get_bytes(run_key(dives_run, "estimate.json"))
        if txt and js:
            return json.loads(js), txt
    e = estimate_run(b, dives_run, log=lambda *_: None, save=False)
    return e, e.get("text", "")


def readout(b, run_id, now=None):
    """the readout of RUN (of its source, for a dives run), as Markdown"""
    now = now or time.time()
    run = b.get_run(run_id)
    if run is None:
        return f"# {run_id}\n\n**NO RUN**: the store has no run {run_id}.\n"
    cfg = run.get("config", {})
    if cfg.get("kind") == "dives" and cfg.get("source_run"):
        return readout(b, cfg["source_run"], now)
    sol = run.get("solver", {})
    status = run.get("status")
    st = _status(b, run_id)
    counts = b.count_by_status(run_id)
    results = _results(b, run_id)
    dives_run = run.get("dives_run")
    est, est_txt = _estimate(b, dives_run) if dives_run else (None, None)
    d_status = (b.get_run(dives_run) or {}).get("status") if dives_run else None
    warnings = []

    # the verdict
    total = sum(counts.values())
    merged = counts.get("merged", 0)
    if status == "finished":
        best = _pick(results, r"best objective \S+") if results else None   # (no backslash in an f-string: Python 3.11)
        verdict = "**COMPLETE**: solved" + (f", {best}" if best else "") + \
            (f" (and its estimate, {dives_run}, below)" if dives_run and d_status == "finished" else "")
    elif dives_run and d_status == "finished":
        verdict = f"**COMPLETE**: the estimate ({dives_run} finished)"
    elif dives_run and est is not None:
        verdict = (f"**PART-WAY**: the estimate's dives run {dives_run} is {d_status}, "
                   f"{est['units'] - est['units_pending']} of {est['units']} dives units merged")
    elif status == "stopped":
        verdict = f"**STOPPED** part-way: {merged} of {total} units merged (qapctl stop)"
    elif status == "paused":
        verdict = f"**PAUSED**: {merged} of {total} units merged" + \
            (f"; its estimate's dives run {dives_run} is not made yet" if dives_run else "")
    else:
        verdict = f"**RUNNING**: {merged} of {total} units merged"
    if status in ("running", "paused") and st.get("updated") and now - st["updated"] > STALE_S \
            and not (dives_run and d_status == "finished"):
        warnings.append(f"STALE: the coordinator last wrote its status {(now - st['updated']) / 60:.0f} minutes ago "
                        f"(an interrupted coordinator? collect.sh shows the instances and why they stopped)")
    if counts.get("dead"):
        warnings.append(f"DEAD UNITS: {counts['dead']} units failed every attempt (the ledger has their errors)")
    if est is not None:
        warnings += est.get("flags", [])

    lines = [f"# {run_id}", "", verdict, ""]
    if warnings:
        lines += [f"- {w}" for w in warnings] + [""]

    # what was run
    lines += ["## The run", "",
              "| | |", "|---|---|",
              f"| problem | {sol.get('instance', '?')}{(' ' + ' '.join(sol['flags'])) if sol.get('flags') else ''} |",
              f"| params | {sol.get('params', '?')} |",
              f"| incumbent | {sol.get('ub') if sol.get('ub') is not None else 'none'} |",
              f"| slots | {cfg.get('slots', '?')}; units of {cfg.get('budget', '?')} CPU s"
              f"{', traced' if cfg.get('trace') else ''} |",
              f"| created | {_when(run.get('created'))} |",
              f"| finished | {_when(run.get('finished_at'))} |" if run.get("finished_at") else f"| status | {status} |"]
    if dives_run:
        dv = (b.get_run(dives_run) or {}).get("config", {}).get("dive", {})
        lines.append(f"| estimate | dives run {dives_run}: {dv.get('b', '?')}, sample {dv.get('sample', 1)} "
                     f"({dv.get('weights', 'unif')}){''.join(f', {k}={v}' for k, v in sorted((dv.get('env') or {}).items()))} |")
    lines += [""]

    # the work so far
    lines += ["## The work", "",
              "| | |", "|---|---|",
              f"| wall | {_wall(st.get('wall_s') or run.get('wall_s'))} |",
              f"| nodes | {st.get('nodes', 0):,} ({st.get('dnn_nodes', 0):,} DNN, {st.get('dnn_fathomed', 0):,} DNN-fathomed) |",
              f"| CPU | {_cpu(st.get('cpu_s'))} |",
              f"| units | {', '.join(f'{k} {v}' for k, v in sorted(counts.items()))} |",
              f"| status written | {_when(st.get('updated'))} |", ""]
    if results:
        keep = [l for l in (_pick(results, r"^\s*Nodes:.*$"), _pick(results, r"^\s*DNN bound:.*$"),
                            _pick(results, r"^.*best objective.*$"), _pick(results, r"^\s*Total Time:.*$")) if l]
        lines += ["## Results", "", "```"] + keep + ["```", ""]
    if est_txt:
        lines += ["## The estimate", "", "```"] + [l for l in est_txt.splitlines() if not l.startswith("|")] + ["```", ""]
        table = [l for l in est_txt.splitlines() if l.startswith("|")]
        if table:
            lines += table + [""]
    return "\n".join(lines)


def write_readout(b, run_id, log=print):
    """runs/RUN/readout.md in the store (RUN's source for a dives run); never fatal"""
    try:
        run = b.get_run(run_id) or {}
        target = run.get("config", {}).get("source_run") if run.get("config", {}).get("kind") == "dives" else run_id
        text = readout(b, target)
        b.put_bytes(run_key(target, "readout.md"), text.encode(), content_type="text/markdown")
        return text
    except Exception as e:       # a readout must never stop a run's end
        log(f"readout: {run_id}: {e!r}")
        return None
