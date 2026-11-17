#!/usr/bin/env python3
"""Compare fixed coarse cuts and regional CP on a shared same-node checkpoint."""

import argparse
import json
import os
from pathlib import Path
import re
import resource
import shlex
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = {
    "tho40_d5": ("tho40t", "tho40c", "26:10,8:0,2:32,6:39,33:24", 240516, 1, 300, 400, 500),
    "tho40_d3": ("tho40t", "tho40c", "26:10,8:0,2:32", 240516, 1, 300, 400, 500),
    "tai60b_d2": ("tai60b", "tai60b", "44:45,43:58", 608215055, 0, 100, 200, 300),
}
SOLVER_CHECK = re.compile(r"^it\s+(\d+) lb .*? cpu (\S+)(?: extra)?$", re.M)
CP_CHECK = re.compile(r"^cp it\s+(\d+) cert (\S+)", re.M)
SOLVE = re.compile(r"^cp solve  lb (\S+).*?status (\S+)  iters (\d+)  cpu (\S+)  certs (\d+) failures (\d+)", re.M)
FINAL = re.compile(r"^FINAL best_lb (\S+) first_lb (\S+)", re.M)


def run(command, log):
    start = resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open("w") as out:
        print(shlex.join(map(str, command)), file=out, flush=True)
        subprocess.run(command, cwd=ROOT, stdout=out, stderr=subprocess.STDOUT,
                       env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}, check=True)
    finish = resource.getrusage(resource.RUSAGE_CHILDREN)
    return finish.ru_utime + finish.ru_stime - start.ru_utime - start.ru_stime


def parse(log, process_cpu=None):
    text = log.read_text()
    lb, status, iters, cpu, certs, failures = SOLVE.search(text).groups()
    best, first = map(float, FINAL.search(text).groups())
    cp_cpu = {int(i): float(t) for i, t in SOLVER_CHECK.findall(text)}
    trajectory = []
    running_best = first
    for it, cert in CP_CHECK.findall(text):
        running_best = max(running_best, float(cert))
        trajectory.append({"it": int(it), "best": running_best, "cpu": cp_cpu[int(it)]})
    # A lazy exit can recover a stronger certified value than any check's
    # trace reported. Attribute it to the end of the solve, at actual CPU.
    if trajectory and float(lb) > running_best:
        terminal = {"it": int(iters), "best": float(lb), "cpu": float(cpu)}
        if terminal["it"] == trajectory[-1]["it"]:
            trajectory[-1] = terminal
        else:
            trajectory.append(terminal)
    coarse = re.findall(r"coarse (\d+) active (\d+) sum_nu (\S+) maxviol (\S+)/", text)
    setup = re.search(r"setup cpu (\S+)", text)
    return {"lb": float(lb), "best": best, "first": first, "status": status, "iters": int(iters),
            "cpu": float(cpu), "process_cpu": process_cpu, "certs": int(certs), "failures": int(failures),
            "coarse_setup_cpu": float(setup[1]) if setup else 0.0,
            "coarse_active_last": int(coarse[-1][1]) if coarse else 0,
            "coarse_active_max": max(int(r[1]) for r in coarse) if coarse else 0,
            "coarse_maxviol_last": float(coarse[-1][3]) if coarse else None,
            "trajectory": trajectory}


def compare(regional, coarse, separation_cpu):
    extra_cpu = separation_cpu + coarse["coarse_setup_cpu"]
    base = {p["it"]: p for p in regional["trajectory"]}
    matched = [{"it": p["it"], "gain": p["best"] - base[p["it"]]["best"],
                "base_cpu": base[p["it"]]["cpu"], "coarse_cpu_with_separation": p["cpu"] + extra_cpu}
               for p in coarse["trajectory"] if p["it"] in base]
    budget = regional["cpu"]
    eligible = [p for p in coarse["trajectory"] if p["cpu"] + extra_cpu <= budget]
    at_budget = eligible[-1]["best"] if eligible else coarse["first"]
    return {"gain_at_final_check": coarse["best"] - regional["best"],
            "max_same_iteration_gain": max(p["gain"] for p in matched),
            "gain_at_regional_cpu_budget": at_budget - regional["best"],
            "matched_checks": [p for p in matched if p["it"] % 100 == 0]}


def run_case(name, work, args):
    inst, recipe, fixes, ub, side, ordinary, warm_iters, step_iters = CASES[name]
    instance = f"qaplib/{inst}.dat"
    checkpoint = work / "seed.cp"
    prefix = work / "seed"
    cuts = work / "coarse.pool"
    recipe = f"infra/phase1/recipes/{recipe}.par"
    command = [str(ROOT / "experiments/dnn_cuts/dnncuts"), "-x", fixes, "-u", str(ub), "-H", "0",
               "-k", "1000", "-cpsep", "100:1000", "-p", recipe, "-side", str(side), "-warmY"]
    if side == 1:
        command += ["-grid", "5:8"]
    else:
        command += ["-regions", "infra/phase1/recipes/tai60b.regions"]
    if not args.reuse_seed:
        seed = command + ["-I", str(ordinary), "-cpiters", str(args.seed_iters or warm_iters), "-savecp", checkpoint,
                          "-dumpY", prefix, instance]
        seed_cpu = run(seed, work / "seed.log")
    else:
        seed_cpu = None
    sep = [sys.executable, str(HERE / "separate.py"), instance, work / "seed.Y.bin", cuts,
           "--side", str(side), "--fixes", fixes, "--samples", str(args.samples),
           "--cut-count", str(args.cut_count), "--sizes", *map(str, args.sizes)]
    if side == 1:
        sep += ["--grid", "5:8"]
    separation_cpu = run(sep, work / "separate.log")
    cut_data = json.loads(cuts.with_suffix(".json").read_text())
    screened = None
    if args.binary_filter:
        filtered = work / "filtered.pool"
        control = work / "local-controls.json"
        screening_cpu = run([sys.executable, str(HERE / "local_control.py"), cuts.with_suffix(".json"), control,
                             "--filter-pool", filtered], work / "local-controls.log")
        separation_cpu += screening_cpu
        screened = json.loads(control.read_text())
        cuts = filtered
    phase = command + ["-cpiters", str(args.continuation_iters or step_iters), "-startcp", checkpoint]
    results = {}
    for arm, flags in (("regional", []), ("coarse", ["-coarse", cuts])):
        log = work / f"{arm}.log"
        process_cpu = None if arm == "regional" and args.reuse_regional else run(phase + flags + [instance], log)
        results[arm] = parse(log, process_cpu)
    summary = {"case": name, "fixes": fixes, "ub": ub, "seed_process_cpu": seed_cpu,
               "seed": parse(work / "seed.log") if not args.reuse_seed else None,
               "n_cuts": len(screened["retained_ranks"]) if screened else cut_data["n_cuts"],
               "capacity": cut_data["capacity"], "binary_screen": screened,
               "sizes": cut_data["sizes"], "samples": cut_data["samples"],
               "separation_cpu": separation_cpu, "separator_internal_cpu": cut_data["cpu_s"], **results,
               "comparison": compare(results["regional"], results["coarse"], separation_cpu)}
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(name, json.dumps(summary["comparison"]), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("cases", nargs="+", choices=CASES)
    parser.add_argument("--reuse-seed", action="store_true")
    parser.add_argument("--reuse-regional", action="store_true", help="use the existing regional.log from this checkpoint")
    parser.add_argument("--seed-iters", type=int)
    parser.add_argument("--continuation-iters", type=int)
    parser.add_argument("--sizes", type=int, nargs="+", default=[4, 5])
    parser.add_argument("--samples", type=int, default=35)
    parser.add_argument("--cut-count", type=int, default=40)
    parser.add_argument("--binary-filter", action="store_true")
    args = parser.parse_args()
    for name in args.cases:
        work = args.output.resolve() / name
        work.mkdir(parents=True, exist_ok=True)
        run_case(name, work, args)


if __name__ == "__main__":
    main()
