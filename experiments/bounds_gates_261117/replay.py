#!/usr/bin/env python3
"""Capture baseline parent states and replay individual QAP nodes."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time


HERE = Path(os.environ.get("QAP_EXPERIMENT_WORK", Path(__file__).resolve().parent)).resolve()
REPO = Path(os.environ.get("QAP_EXPERIMENT_SOURCE", HERE.parent / "qap")).resolve()
SEED = Path(os.environ.get("QAP_EXPERIMENT_SEED", HERE.parent / "qap-experiment-seed")).resolve()
STATES = HERE / "states"
RECIPES = REPO / "infra/phase1/recipes"
PURPOSES = {
    "tho40-d3-hard": ("tho40", "26:10,8:0,2:32"),
    "tho40-d3-mid": ("tho40", "26:10,8:0,3:32"),
    "tho40-d4-hard-parent": ("tho40", "26:10,8:0,2:32,6:39"),
    "tho40-d4-subtree": ("tho40", "26:10,8:0,2:32,33:24"),
    "tho40-d5-hard": ("tho40", "26:10,8:0,2:32,6:39,33:24"),
    "tho40-d5-tight": ("tho40", "26:10,8:0,2:32,33:24,39:8"),
    "tai60b-d2": ("tai60b", "44:45,43:58"),
    "tai60b-d3-tight": ("tai60b", "44:45,43:58,1:35"),
    "nug12": ("nug12", ""),
}


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def settings(instance, folder, cap=None):
    if instance == "tai60b":
        instance_path = REPO / "qaplib/tai60b.dat"
        recipe_path = RECIPES / "tai60b.par"
        extra = ["-u", "608215055"]
    else:
        instance_path = REPO / f"qaplib/{'tho40t' if instance == 'tho40' else 'nug12'}.dat"
        recipe_path = RECIPES / "tho40c.par"
        extra = ["-u", "240516", "-s", str(REPO / "latest/SYM/tho40t.sym")] if instance == "tho40" else ["-u", "580"]

    source = recipe_path.read_text()
    # Preserve the recipe's relative cluster-region path and its settings hash.
    overrides = {}
    if instance == "nug12":
        overrides.update({"cut.grid": "3:4", "dnn.depths": "0:11", "cut.zone": "1", "dnn.giveup": "40"})
    if cap is not None:
        overrides["cut.cp_iters"] = str(cap)
    if instance == "nug12" and cap is None:
        overrides["cut.cp_iters"] = "400"
    for key, value in overrides.items():
        pattern = rf"^{re.escape(key)}\s*=[^\n]*$"
        spec = f"{key} = {value}"
        source, matches = re.subn(pattern, spec, source, flags=re.M)
        if not matches:
            source = f"{source.rstrip()}\n{spec}\n"
    if overrides:
        recipe_path = folder / "params.par"
        recipe_path.write_text(source)
    return ["-p", str(recipe_path), *extra], instance_path, recipe_path


def depth(fixes):
    return len(fixes.split(",")) if fixes else 0


def assignments(fixes):
    return {int(i): int(j) for i, j in (pair.split(":") for pair in fixes.split(",") if pair)}


def state_content(path):
    lines = path.read_text().splitlines(keepends=True)
    dimension = int(next(line for line in lines if line.startswith("n ")).split()[1])
    return lines, dimension


def slice_state(source, target, fixes):
    lines, n = state_content(source)
    wanted = assignments(fixes)
    matching = []
    for line in lines:
        fields = line.split()
        if not fields or fields[0] != "node" or int(fields[1]) != len(wanted):
            continue
        actual = {i: int(j) for i, j in enumerate(fields[5:5 + n]) if int(j) >= 0}
        if actual == wanted:
            matching.append(line)
    if len(matching) != 1:
        raise ValueError(f"expected one matching node {fixes!r} in {source}, found {len(matching)}")
    start = next(i for i, line in enumerate(lines) if line.startswith("nodes "))
    warm_id = int(matching[0].split()[-1])
    warms = [line for line in lines if line.startswith("warm ")]
    if not 0 <= warm_id < len(warms):
        raise ValueError(f"missing parent warm line for {fixes!r}")
    pool, token = warms[warm_id].split()[-2:]
    if pool == "-" or token == "-":
        raise ValueError(f"missing inherited cut pool or DNN state for {fixes!r}")
    warm_blob = (source.parent / token[1:]).resolve() if token.startswith("@") else None
    if source.parent != target.parent:
        raise ValueError("state and its slice must share their directory for scripts using relative warm references")
    temporary = target.with_suffix(".tmp")
    temporary.write_text("".join([*lines[:start], "nodes 1\n", matching[0], "end\n"]))
    temporary.replace(target)
    facts = {"fixes": fixes, "state": str(target), "state_sha256": digest(target),
             "source_state": str(source), "warm_id": warm_id, "pool_bytes": len(pool)}
    if warm_blob:
        facts.update(warm_blob=str(warm_blob), warm_blob_sha256=digest(warm_blob))
    target.with_suffix(".json").write_text(json.dumps(facts, indent=2) + "\n")


def wait_for_solver(solver, folder, start):
    snapshot = folder / "progress.json"
    stat = Path(f"/proc/{solver.pid}/stat")
    ticks = os.sysconf("SC_CLK_TCK")
    page_kib = os.sysconf("SC_PAGE_SIZE") / 1024
    while True:
        try:
            fields = stat.read_text().rsplit(") ", 1)[1].split()
            info = {"status": fields[0], "cpu_seconds": (int(fields[11]) + int(fields[12])) / ticks,
                    "rss_kib": int(fields[21]) * page_kib, "elapsed_seconds": time.monotonic() - start}
            snapshot.write_text(json.dumps(info, indent=2) + "\n")
        except (OSError, IndexError, ValueError):
            pass
        try:
            return solver.wait(timeout=15)
        except subprocess.TimeoutExpired:
            continue


def process(binary, instance, folder, arguments, changes, cap=None, parent_state=None):
    folder.mkdir(parents=True)
    bin_path = binary.resolve()
    base_arguments, instance_path, recipe_path = settings(instance, folder, cap)
    command = [str(bin_path), *base_arguments, *arguments, str(instance_path)]
    environment = {key: val for key, val in os.environ.items() if not key.startswith(("QAP_", "DNN_"))}
    seeds = {key: "1" for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS")}
    seeds.update(QAP_DNN_TRACE=str(folder / "dnn.txt"), QAP_NODE_TRACE=str(folder / "nodes.txt"), QAP_DNN_TRACE_HIST="1")
    seeds.update(changes)
    environment.update(seeds)
    manifest = {"command": command, "cwd": str(bin_path.parent), "environment": seeds,
                "commit": subprocess.check_output(["git", "-C", str(bin_path.parent), "rev-parse", "HEAD"], text=True).strip(),
                "binary_sha256": digest(bin_path), "instance_sha256": digest(instance_path),
                "recipe_sha256": digest(recipe_path), "input_state_sha256": digest(parent_state) if parent_state else None}
    if parent_state:
        source = json.loads(parent_state.with_suffix(".json").read_text())
        if source["state_sha256"] != manifest["input_state_sha256"]:
            raise ValueError(f"saved state changed: {parent_state}")
        if warm_blob := source.get("warm_blob"):
            manifest["input_warm_blob_sha256"] = digest(Path(warm_blob))
            if source["warm_blob_sha256"] != manifest["input_warm_blob_sha256"]:
                raise ValueError(f"warm state changed: {warm_blob}")
    if instance == "tai60b":
        manifest["regions_sha256"] = digest(RECIPES / "tai60b.regions")
        if digest(bin_path.parent / "../infra/phase1/recipes/tai60b.regions") != manifest["regions_sha256"]:
            raise ValueError("solver tree cluster regions differ from baseline")
    artifact = folder / "run.json"
    artifact.write_text(json.dumps(manifest, indent=2) + "\n")
    start = time.monotonic()
    initial = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (folder / "stdout.txt").open("w") as stdout, (folder / "stderr.txt").open("w") as stderr:
        with subprocess.Popen(command, cwd=bin_path.parent, env=environment, stdout=stdout, stderr=stderr) as solver:
            exit_code = wait_for_solver(solver, folder, start)
    final = resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest.update(wall_seconds=time.monotonic() - start, user_seconds=final.ru_utime - initial.ru_utime,
                    system_seconds=final.ru_stime - initial.ru_stime, max_rss_kib=final.ru_maxrss, returncode=exit_code)
    artifact.write_text(json.dumps(manifest, indent=2) + "\n")
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, command)


def capture(args):
    STATES.mkdir(exist_ok=True)
    warm = STATES / "warm"
    warm.mkdir(exist_ok=True)
    folder = HERE / "capture" / args.tag
    out = STATES / f"{args.tag}-all.state"
    parent = STATES / f"{args.parent}.state" if args.parent else None
    opening = ["-r", str(parent)] if parent else ["-F", args.fixes]
    options = {"QAP_FORCE_LINE": f"{depth(args.fixes)}:{args.line}", "QAP_DNN_WARM_DIR": str(warm)}
    process(args.binary, args.instance, folder, [*opening, "-d", str(depth(args.fixes)), "-o", str(out)], options, parent_state=parent)
    print(out)


def extract(args):
    _, fixes = PURPOSES[args.node]
    target = STATES / f"{args.node}.state"
    slice_state(STATES / f"{args.source}-all.state", target, fixes)
    print(target)


def run(args):
    instance, fixes = PURPOSES[args.node]
    parent = STATES / f"{args.node}.state" if fixes else None
    opening = ["--allow-settings-change", "-r", str(parent)] if parent else []
    depth_limit = [] if args.complete else ["-d", str(depth(fixes))]
    settings_env = dict(value.split("=", 1) for value in args.env)
    out = HERE / "runs" / args.arm / args.node
    process(args.binary, instance, out, [*opening, *depth_limit], settings_env, args.cap, parent)
    print(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    seed = commands.add_parser("capture")
    seed.add_argument("tag")
    seed.add_argument("instance", choices=["tho40", "tai60b"])
    seed.add_argument("fixes")
    seed.add_argument("line", help="forced branching line, e.g. c32 or r43")
    seed.add_argument("--parent", help="existing one-node warm-state tag under states/")
    seed.add_argument("--binary", type=Path, default=SEED / "latest/bnbdrv")
    seed.set_defaults(callback=capture)
    select = commands.add_parser("extract")
    select.add_argument("source", help="capture tag without the -all.state suffix")
    select.add_argument("node", choices=PURPOSES)
    select.set_defaults(callback=extract)
    arm = commands.add_parser("run")
    arm.add_argument("arm")
    arm.add_argument("node", choices=PURPOSES)
    arm.add_argument("binary", type=Path)
    arm.add_argument("--cap", type=int, help="override the cut phase iteration cap")
    arm.add_argument("--env", action="append", default=[], metavar="KEY=VALUE")
    arm.add_argument("--complete", action="store_true", help="solve the whole saved subtree")
    arm.set_defaults(callback=run)
    args = parser.parse_args()
    args.callback(args)


if __name__ == "__main__":
    main()
