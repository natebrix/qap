#!/usr/bin/env python3
"""Every run's readout, for the record: the GitHub workflow phase1-readouts runs this hourly with
a read-only role and commits what it writes to the branch phase1-readouts (no approval: it reads,
and writes only to that branch). Also usable from the laptop.

  python3 infra/phase1/readouts.py --backend aws:BUCKET?region=R --out DIR [--since-days 30]

Writes DIR/board.md (one line per run: status, verdict, when) and DIR/runs/RUN/readout.md for each
run (coord/readout.py): the coordinator's own readout for a finished or stopped run, else one made
now (a running estimate's part-way estimate included). Dives runs are their sources' readouts.
Nothing is written to the store. The output names no bucket (the account number is in its name).
"""
import argparse
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from coord.backend import open_backend  # noqa: E402
from coord.common import run_key  # noqa: E402
from coord.readout import _primary, readout  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--since-days", type=float, default=30, help="runs created or finished within this many days")
    a = ap.parse_args()
    b = open_backend(a.backend)
    now = time.time()
    rows = []
    for run in b.list_runs():
        rid = run.get("run_id", "")
        # a dives run that is its source's estimate is in its source's readout; a second one (a re-dive of
        # a paused frontier, start.sh FROM= DIVES=) has its own (coord/readout.py)
        if not rid or (run.get("config", {}).get("kind") == "dives" and _primary(b, rid, run)):
            continue
        last = max(run.get("finished_at") or 0, run.get("created") or 0, run.get("stopped_at") or 0)
        if run.get("status") in ("finished", "stopped") and now - last > a.since_days * 86400:
            continue
        stored = b.get_bytes(run_key(rid, "readout.md")) if run.get("status") in ("finished", "stopped") else None
        try:
            text = stored.decode() if stored else readout(b, rid, now)
        except Exception as e:      # one run's readout must not stop the others'
            text = f"# {rid}\n\n**UNKNOWN**: the readout could not be made ({e!r})\n"
        d = os.path.join(a.out, "runs", rid)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "readout.md"), "w") as f:
            f.write(text)
        verdict = next((l for l in text.splitlines() if l.startswith("**")), "")
        flags = [l[2:].split(":")[0] for l in text.splitlines() if l.startswith("- ") and l[2:3].isupper()]
        rows.append((run.get("created") or 0, rid, run.get("status", "?"), verdict, flags,
                     run.get("finished_at") or run.get("stopped_at")))
    rows.sort(reverse=True)
    when = lambda t: time.strftime("%Y-%m-%d %H:%M", time.gmtime(t)) if t else "-"
    lines = ["# Phase 1 runs", "",
             f"Every run of the last {a.since_days:g} days, newest first (`infra/phase1/readouts.py`); "
             f"each run's readout is `runs/RUN/readout.md`. Times UTC.", "",
             "| run | status | verdict | warnings | created | ended |", "|---|---|---|---|---|---|"]
    for created, rid, status, verdict, flags, ended in rows:
        lines.append(f"| [{rid}](runs/{rid}/readout.md) | {status} | {verdict} | {', '.join(flags) or '-'} | "
                     f"{when(created)} | {when(ended)} |")
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "board.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"{len(rows)} runs: {os.path.join(a.out, 'board.md')}")


if __name__ == "__main__":
    main()
