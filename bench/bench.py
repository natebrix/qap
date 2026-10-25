#!/usr/bin/env python3
"""Run the QAP branch-and-bound solver over a suite of instances and compare runs.

  python3 bench/bench.py run --suite quick --rev master --label base
  python3 bench/bench.py run --suite quick --build --label mine
  python3 bench/bench.py compare base mine --log --note "what changed"

Run it from anywhere; paths are resolved relative to the repository root. See
RUNNING.md for what "apples to apples" means here and which settings are fixed,
and CONTRIBUTING.md for how collaborators record results in bench/LOG.csv.
Uses only the Python standard library.
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import platform
import re
import socket
import statistics
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER_DIR = os.path.join(ROOT, "latest")
SOLVER = os.path.join(SOLVER_DIR, "bnbdrv")
MESCHACH_DIR = os.path.join(ROOT, "src", "MESCHACH")
QAPLIB = os.path.join(ROOT, "qaplib")
QAPSOLN = os.path.join(ROOT, "QAPSOLN")
SUITES = os.path.join(ROOT, "bench", "suites")
RESULTS = os.path.join(ROOT, "bench", "results")
WORKTREES = os.path.join(ROOT, "bench", ".worktrees")
LOG = os.path.join(ROOT, "bench", "LOG.csv")

DEFAULT_PARAMS = "PARAM/default.par"

CSV_FIELDS = [
    "instance", "n", "status", "optimum", "objective", "upper_bound",
    "nodes", "fathomed", "enumerated", "root_bound", "fw_iterations",
    "solver_time", "wall_time", "wall_time_min", "repeats", "params", "sym_file",
]


# ---------------------------------------------------------------- helpers

def sh(cmd, cwd=ROOT):
    try:
        return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def sha1_file(path):
    with open(path, "rb") as f:
        return hashlib.sha1(f.read()).hexdigest()[:12]


def read_optimum(name):
    """Return (n, optimal objective) from QAPSOLN/<name>.sln, or (None, None)."""
    path = os.path.join(QAPSOLN, name + ".sln")
    if not os.path.exists(path):
        return None, None
    with open(path) as f:
        tok = f.read().split()
    return int(tok[0]), int(float(tok[1]))


def read_suite(suite):
    path = suite if os.path.exists(suite) else os.path.join(SUITES, suite + ".txt")
    if not os.path.exists(path):
        sys.exit(f"suite not found: {suite} (looked in {SUITES})")
    names = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(line)
    return names


def resolve_params(spec, name):
    """Return params file path (relative to latest/) for an instance.

    'auto' means the solver's own default (PARAM/<name>.par), falling back to
    DEFAULT_PARAMS when that file does not exist. The solver does not check
    that a params file opened, so we always pass an existing one explicitly.
    """
    if spec == "auto":
        own = f"PARAM/{name}.par"
        spec = own if os.path.exists(os.path.join(SOLVER_DIR, own)) else DEFAULT_PARAMS
    full = spec if os.path.isabs(spec) else os.path.join(SOLVER_DIR, spec)
    if not os.path.exists(full):
        sys.exit(f"params file not found: {full}")
    return spec


def use_tree(root):
    """Point the solver paths at the source tree rooted at `root`."""
    global SOLVER_DIR, SOLVER, MESCHACH_DIR
    SOLVER_DIR = os.path.join(root, "latest")
    SOLVER = os.path.join(SOLVER_DIR, "bnbdrv")
    MESCHACH_DIR = os.path.join(root, "src", "MESCHACH")


def checkout_rev(rev):
    """Check out `rev` in a cached worktree under bench/.worktrees/; return (path, sha)."""
    sha = sh(["git", "rev-parse", "--verify", rev + "^{commit}"])
    if not sha:
        sys.exit(f"unknown git revision: {rev}")
    path = os.path.join(WORKTREES, sha[:12])
    if not os.path.isdir(path):
        os.makedirs(WORKTREES, exist_ok=True)
        r = subprocess.run(["git", "worktree", "add", "--detach", path, sha],
                           cwd=ROOT, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"git worktree add failed:\n{r.stderr}")
    return path, sha


def build():
    print(f"building Meschach and bnbdrv in {os.path.relpath(os.path.dirname(SOLVER_DIR), ROOT) or '.'} ...",
          flush=True)
    for d, target in ((MESCHACH_DIR, "basic"), (SOLVER_DIR, "all")):
        r = subprocess.run(["make", target], cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            sys.stderr.write(r.stdout + r.stderr)
            sys.exit(f"build failed in {d}\n(commits before the Meschach isatty fix, "
                     f"e.g. master at 16a72e8, do not build with modern g++)")


def machine_info():
    cpu = ""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except OSError:
        cpu = platform.processor()
    return {
        "hostname": socket.gethostname(),
        "cpu": cpu,
        "cpus": os.cpu_count(),
        "platform": platform.platform(),
        "compiler": sh(["g++", "--version"]).splitlines()[0] if sh(["g++", "--version"]) else "",
    }


def git_info(rev_sha=None, rev=None):
    if rev_sha:
        return {"commit": rev_sha, "branch": rev, "dirty": False}
    return {
        "commit": sh(["git", "rev-parse", "HEAD"]),
        "branch": sh(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(sh(["git", "status", "--porcelain", "--untracked-files=no"])),
    }


# ---------------------------------------------------------------- parsing

PATTERNS = {
    "solver_time": r"Total Time:\s*([0-9.eE+-]+)",
    "nodes_line": r"Nodes:\s*(\d+)\s*/\s*(\d+)\s*/\s*(\d+)",
    "root_bound": r"root bound:\s*([0-9.eE+-]+)",
    "fw_iterations": r"fw iterations:\s*(\d+)",
}


def parse_output(text):
    out = {}
    m = re.search(PATTERNS["solver_time"], text)
    out["solver_time"] = float(m.group(1)) if m else None
    m = re.search(PATTERNS["nodes_line"], text)
    if m:
        out["nodes"], out["fathomed"], out["enumerated"] = map(int, m.groups())
    m = re.search(PATTERNS["root_bound"], text)
    out["root_bound"] = float(m.group(1)) if m else None
    m = re.search(PATTERNS["fw_iterations"], text)
    out["fw_iterations"] = int(m.group(1)) if m else None
    found = re.findall(r"New objective value:\s*([0-9.eE+-]+)", text)
    out["objective"] = int(float(found[-1])) if found else None
    out["finished"] = "---- RESULTS ----" in text
    return out


# ---------------------------------------------------------------- run

def run_one(name, params, ub, timeout, log_path):
    dat = os.path.join(QAPLIB, name + ".dat")
    cmd = [SOLVER, "-p", params]
    if ub is not None:
        cmd += ["-u", str(ub)]
    cmd.append(dat)
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, cwd=SOLVER_DIR, capture_output=True, text=True,
                           timeout=timeout)
        text, rc, timed_out = r.stdout + r.stderr, r.returncode, False
    except subprocess.TimeoutExpired as e:
        out = e.stdout or ""
        text = out.decode(errors="replace") if isinstance(out, bytes) else out
        rc, timed_out = None, True
    wall = time.perf_counter() - t0
    with open(log_path, "w") as f:
        f.write("$ (cd latest && " + " ".join(cmd) + ")\n")
        f.write(text)
    res = parse_output(text)
    res["wall_time"] = wall
    res["status"] = ("timeout" if timed_out else
                     "ok" if rc == 0 and res["finished"] else f"error(rc={rc})")
    return res


def cmd_run(args):
    rev_sha = None
    if args.rev:
        path, rev_sha = checkout_rev(args.rev)
        use_tree(path)
        if not os.path.exists(SOLVER):
            build()  # first use of this worktree; later runs reuse the binary
    if args.build:
        build()
    if not os.path.exists(SOLVER):
        sys.exit(f"{SOLVER} not found; build first (see RUNNING.md) or pass --build")

    names = args.instances or read_suite(args.suite)
    label = args.label or datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(RESULTS, label)
    if os.path.exists(outdir) and not args.force:
        sys.exit(f"{outdir} exists; choose another --label or pass --force")
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)

    params_used = {}
    meta = {
        "label": label,
        "started": datetime.datetime.now().isoformat(timespec="seconds"),
        "suite": args.suite if not args.instances else None,
        "instances": names,
        "params": args.params,
        "upper_bound": args.ub,
        "repeats": args.repeats,
        "timeout": args.timeout,
        "binary_sha1": sha1_file(SOLVER),
        "git": git_info(rev_sha, args.rev),
        "machine": machine_info(),
        "params_sha1": params_used,
        "notes": args.notes,
    }

    rows = []
    for name in names:
        n, opt = read_optimum(name)
        if not os.path.exists(os.path.join(QAPLIB, name + ".dat")):
            print(f"{name:>10}: missing qaplib/{name}.dat, skipped")
            continue
        params = resolve_params(args.params, name)
        params_used[params] = sha1_file(os.path.join(SOLVER_DIR, params))
        if args.ub == "opt":
            if opt is None:
                print(f"{name:>10}: no QAPSOLN/{name}.sln, cannot use --ub opt, skipped")
                continue
            ub = opt
        elif args.ub == "none":
            ub = None
        else:
            ub = float(args.ub)

        walls, res = [], None
        for rep in range(args.repeats):
            log = os.path.join(outdir, "logs", f"{name}.log" if rep == 0 else f"{name}.{rep}.log")
            res_rep = run_one(name, params, ub, args.timeout, log)
            walls.append(res_rep["wall_time"])
            if res is None:
                res = res_rep
            elif res_rep.get("nodes") != res.get("nodes"):
                print(f"{name:>10}: WARNING node count changed between repeats "
                      f"({res.get('nodes')} vs {res_rep.get('nodes')}); run is not deterministic")
            if res_rep["status"] != "ok":
                break

        status = res["status"]
        # Correctness: with --ub opt nothing better than the optimum may appear;
        # without it, the search must end on the optimum.
        if status == "ok" and opt is not None:
            obj = res["objective"]
            if obj is not None and obj < opt:
                status = "wrong(below-opt)"
            elif ub is None and obj != opt:
                status = "wrong(not-opt)"

        row = {
            "instance": name, "n": n, "status": status, "optimum": opt,
            "objective": res["objective"], "upper_bound": ub,
            "nodes": res.get("nodes"), "fathomed": res.get("fathomed"),
            "enumerated": res.get("enumerated"), "root_bound": res["root_bound"],
            "fw_iterations": res["fw_iterations"], "solver_time": res["solver_time"],
            "wall_time": round(statistics.median(walls), 4),
            "wall_time_min": round(min(walls), 4), "repeats": len(walls),
            "params": params, "sym_file": "auto",
        }
        rows.append(row)
        t = f"{row['solver_time']:.3f}s" if row["solver_time"] is not None else "-"
        print(f"{name:>10}: {status:<8} nodes={row['nodes'] if row['nodes'] is not None else '-':>10}  "
              f"time={t:>10}  wall={row['wall_time']:.2f}s", flush=True)

    with open(os.path.join(outdir, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    meta["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    with open(os.path.join(outdir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    bad = [r["instance"] for r in rows if r["status"] != "ok"]
    print(f"\nwrote {os.path.relpath(outdir, ROOT)}/  ({len(rows)} instances, "
          f"{len(bad)} not ok{': ' + ', '.join(bad) if bad else ''})")
    return 1 if any(r["status"].startswith(("wrong", "error")) for r in rows) else 0


# ---------------------------------------------------------------- compare

def load(label):
    d = label if os.path.isdir(label) else os.path.join(RESULTS, label)
    with open(os.path.join(d, "meta.json")) as f:
        meta = json.load(f)
    with open(os.path.join(d, "results.csv")) as f:
        rows = {r["instance"]: r for r in csv.DictReader(f)}
    return meta, rows


def geomean(xs):
    xs = [x for x in xs if x > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def cmd_compare(args):
    meta_a, a = load(args.a)
    meta_b, b = load(args.b)

    # Apples-to-apples guard: say loudly when the runs differ in ways that
    # make time (or node) ratios meaningless.
    warnings = []
    for key in ("upper_bound", "params"):
        if meta_a.get(key) != meta_b.get(key):
            warnings.append(f"{key}: {meta_a.get(key)} vs {meta_b.get(key)}")
    if meta_a.get("params_sha1") != meta_b.get("params_sha1"):
        warnings.append("params file contents differ")
    for key in ("hostname", "cpu"):
        if meta_a["machine"].get(key) != meta_b["machine"].get(key):
            warnings.append(f"machine {key}: {meta_a['machine'].get(key)} vs "
                            f"{meta_b['machine'].get(key)} (times not comparable)")

    print(f"A = {args.a}  ({meta_a['git']['commit'][:8]}{' dirty' if meta_a['git']['dirty'] else ''}, "
          f"params={meta_a['params']}, ub={meta_a['upper_bound']})")
    print(f"B = {args.b}  ({meta_b['git']['commit'][:8]}{' dirty' if meta_b['git']['dirty'] else ''}, "
          f"params={meta_b['params']}, ub={meta_b['upper_bound']})")
    for w in warnings:
        print(f"  WARNING: {w}")
    print()

    tcol = "solver_time" if args.time == "solver" else "wall_time"
    hdr = f"{'instance':>10} {'nodes A':>11} {'nodes B':>11} {'B/A':>7}   {'time A':>9} {'time B':>9} {'B/A':>7}  status"
    print(hdr)
    print("-" * len(hdr))
    node_ratios, time_ratios = [], []
    for name in [n for n in a if n in b]:
        ra, rb = a[name], b[name]
        st = ra["status"] if ra["status"] == rb["status"] else f"{ra['status']}->{rb['status']}"
        both_ok = ra["status"] == "ok" and rb["status"] == "ok"
        na, nb = ra["nodes"], rb["nodes"]
        ta, tb = ra[tcol], rb[tcol]
        nr = tr = ""
        if both_ok and na and nb:
            node_ratios.append(int(nb) / int(na))
            nr = f"{node_ratios[-1]:7.3f}"
        if both_ok and ta and tb and float(ta) > 0:
            time_ratios.append(float(tb) / float(ta))
            tr = f"{time_ratios[-1]:7.3f}"
        print(f"{name:>10} {na or '-':>11} {nb or '-':>11} {nr:>7}   "
              f"{float(ta) if ta else '-':>9} {float(tb) if tb else '-':>9} {tr:>7}  {st}")
    only = sorted(set(a) ^ set(b))
    if only:
        print(f"\nnot in both runs (ignored): {', '.join(only)}")
    print(f"\ngeometric mean B/A over {len(node_ratios)} instances solved in both:"
          f"  nodes {geomean(node_ratios):.3f}   {args.time} time {geomean(time_ratios):.3f}"
          f"   (< 1 means B is better)")

    if args.log:
        if meta_a["instances"] != meta_b["instances"]:
            warnings.append("the two runs used different instance lists")
        for side, meta in (("A", meta_a), ("B", meta_b)):
            if meta["git"]["dirty"]:
                warnings.append(f"run {side} was made with uncommitted changes")
        if warnings and not args.force_log:
            print("\nnot logged:")
            for w in warnings:
                print(f"  {w}")
            print("fix these (commit your change, re-run) or pass --force-log to record them")
            return 1
        ok = lambda rows: sum(r["status"] == "ok" for r in rows.values())
        bad = [n for n in b if b[n]["status"].startswith(("wrong", "error"))]
        entry = {
            "date": datetime.date.today().isoformat(),
            "author": sh(["git", "config", "user.name"]),
            "suite": meta_b.get("suite") or "custom",
            "params": meta_b["params"],
            "upper_bound": meta_b["upper_bound"],
            "base_commit": meta_a["git"]["commit"][:10],
            "base_label": meta_a["label"],
            "cand_commit": meta_b["git"]["commit"][:10],
            "cand_branch": meta_b["git"]["branch"],
            "cand_label": meta_b["label"],
            "instances": len(meta_b["instances"]),
            "ok_base": ok(a),
            "ok_cand": ok(b),
            "nodes_ratio": f"{geomean(node_ratios):.4f}",
            "time_ratio": f"{geomean(time_ratios):.4f}",
            "worst_nodes_ratio": f"{max(node_ratios):.4f}" if node_ratios else "",
            "wrong_or_error": " ".join(bad),
            "host": meta_b["machine"]["hostname"],
            "cpu": meta_b["machine"]["cpu"],
            "warnings": "; ".join(warnings),
            "note": args.note,
        }
        new = not os.path.exists(LOG)
        with open(LOG, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(entry))
            if new:
                w.writeheader()
            w.writerow(entry)
        print(f"\nappended to {os.path.relpath(LOG, ROOT)}; commit it with your change")
    return 0


# ---------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run a suite and record results")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--suite", default="quick",
                   help="suite name in bench/suites/ or a path (default: quick)")
    g.add_argument("--instances", nargs="+", help="instance names instead of a suite")
    r.add_argument("--params", default=DEFAULT_PARAMS,
                   help=f"params file relative to latest/, used for every instance "
                        f"(default: {DEFAULT_PARAMS}); 'auto' = per-instance PARAM/<name>.par")
    r.add_argument("--ub", default="opt",
                   help="initial upper bound: 'opt' = known optimum from QAPSOLN (default), "
                        "'none' = solver default, or a number")
    r.add_argument("--repeats", type=int, default=1,
                   help="runs per instance; wall time is the median (default: 1)")
    r.add_argument("--timeout", type=float, default=600, help="seconds per run (default: 600)")
    r.add_argument("--label", help="results go to bench/results/<label>/ (default: timestamp)")
    r.add_argument("--notes", default="", help="free text saved in meta.json")
    r.add_argument("--build", action="store_true", help="run make before benchmarking")
    r.add_argument("--rev", help="benchmark this git commit/branch instead of the working "
                                 "tree (checked out and built under bench/.worktrees/)")
    r.add_argument("--force", action="store_true", help="overwrite an existing label")

    c = sub.add_parser("compare", help="compare two recorded runs (A = baseline)")
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("--time", choices=["solver", "wall"], default="solver",
                   help="which time to compare (default: solver-reported)")
    c.add_argument("--log", action="store_true", help="append a summary row to bench/LOG.csv")
    c.add_argument("--note", default="", help="what changed (saved with --log)")
    c.add_argument("--force-log", action="store_true",
                   help="log even when the runs are not comparable (warnings are recorded)")

    args = p.parse_args()
    sys.exit(cmd_run(args) if args.cmd == "run" else cmd_compare(args))


if __name__ == "__main__":
    main()
