#!/usr/bin/env python3
"""Report certified bound, cost and stopping outcomes from replay.py runs."""

import argparse
import csv
import json
from pathlib import Path
import re
import sys


FIRST = re.compile(r"depth (\d+) fixes(.*?) \| first lb (\S+) primal (\S+) status (\w+) iters (\d+) cpu (\S+) warm (\d+) pool (\d+) kept (\d+) dropped (\d+)")
PHASE = re.compile(r"\| cut 1 best (\S+) rounds (\d+) iters (\d+) cpu (\S+) status (\w+) cp (\d+) active (\d+) stalled (\d+) gaveup (\d+)")
CUTOFF = re.compile(r"\| cutoff (\S+) fathomed (\d+) pool_out (\d+)")
CHECKS = (100, 200, 400, 1000)


def dnn(line):
    match = FIRST.search(line)
    if not match:
        raise ValueError(f"unrecognized first DNN solve line {line!r}")
    depth, fixes, lb, primal, status, iterations, cpu, warm, pool, kept, dropped = match.groups()
    cutoff, fathomed, pool_out = CUTOFF.search(line).groups()
    result = dict(depth=int(depth), fixes=",".join(fixes.split()), first_bound=float(lb), first_primal=float(primal),
                  first_status=status, first_iters=int(iterations), first_cpu=float(cpu), warm=bool(int(warm)),
                  pool=bool(int(pool)), kept=int(kept), dropped=int(dropped), cutoff=float(cutoff),
                  fathomed=bool(int(fathomed)), pool_out_bytes=int(pool_out))
    if phase := PHASE.search(line):
        bound, rounds, iterations, cpu, status, cp, active, stalled, gaveup = phase.groups()
        result.update(phase_bound=float(bound), phase_rounds=int(rounds), phase_iters=int(iterations), phase_cpu=float(cpu),
                      phase_status=status, cp=bool(int(cp)), phase_active=int(active), stalled=bool(int(stalled)),
                      gaveup=bool(int(gaveup)), shortfall=float(cutoff) - float(bound))
    if " | hist " in line:
        hist = line.split(" | hist ", 1)[1].split()
        checks = {int(it): float(lb) for it, lb in (entry.split(":") for entry in hist)}
        for query in CHECKS:
            eligible = [iteration for iteration in checks if iteration <= query]
            if eligible and query <= result["phase_iters"]:
                latest = max(eligible)
                result[f"bound_at_{query}"] = checks[latest]
                result[f"checked_at_{query}"] = latest
    return result


def summarize(folder):
    metadata = json.loads((folder / "run.json").read_text())
    result = {"arm": folder.parent.name, "node": folder.name, "commit": metadata["commit"],
              "input_state_sha256": metadata["input_state_sha256"], "wall_seconds": metadata.get("wall_seconds"),
              "process_cpu": metadata.get("user_seconds", 0.) + metadata.get("system_seconds", 0.),
              "max_rss_kib": metadata.get("max_rss_kib"), "returncode": metadata.get("returncode")}
    trace = folder / "nodes.txt"
    if not trace.exists():
        return result
    nodes = [line.split() for line in trace.read_text().splitlines()]
    evaluated = [entry for entry in nodes if entry[7] != "P"]
    result.update(new_nodes=len(evaluated), popped_fathoms=len(nodes) - len(evaluated),
                  new_node_cpu=sum(float(entry[11]) for entry in evaluated))
    if evaluated:
        result.update(root_outcome=evaluated[0][7], root_qpb_bound=evaluated[0][3], root_node_cpu=float(evaluated[0][11]))
    dnn_trace = folder / "dnn.txt"
    if dnn_trace.exists():
        solves = [dnn(line) for line in dnn_trace.read_text().splitlines()]
        if solves:
            for key in ("first_iters", "first_cpu", "phase_iters", "phase_cpu"):
                result[f"total_{key}"] = sum(solve.get(key, 0) for solve in solves)
            result.update(solves[0])
    log = (folder / "stdout.txt").read_text()
    if objective := re.findall(r"(?im)^\s*(?:new objective value|optimal objective(?: value)?|objective value|best objective(?: value)?)\s*[:=]\s*([-+.eE\d]+)", log):
        result["final_objective"] = float(objective[-1])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, nargs="+", help="run folders (each containing run.json)")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    rows = [summarize(folder) for folder in args.folder]
    if args.json:
        print(json.dumps(rows, indent=2))
        return
    writer = csv.DictWriter(sys.stdout, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
