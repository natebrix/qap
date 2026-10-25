#!/usr/bin/env python3
"""Check the branch-and-bound solver against full enumeration on small random QAPs.

  python3 tests/brute_force.py            # default grid (~1 min)
  python3 tests/brute_force.py --quick    # fewer instances, for a fast check
  python3 tests/brute_force.py --solver bench/.worktrees/<sha>/latest/bnbdrv

For each random instance (n = 2..8) the optimum is found by enumerating every
permutation. The solver is then run under several params files and initial
upper bounds, and must:

  * report a root lower bound no greater than the optimum,
  * never report an objective below the optimum,
  * finish on exactly the optimum when started with no upper bound or with
    upper bound optimum + 2,
  * report an objective that some permutation actually attains.

Instances with n = 2..8 are drawn from a fixed seed, so a run is repeatable.

Instance families mimic QAPLIB structure and deliberately contain many ties:
nug-like grids (symmetric, even objectives), small-integer symmetric matrices,
tai-b-like (A symmetric, B asymmetric: odd objectives possible), and symmetric
data plus a linear cost matrix C. Uses only the Python standard library.
"""

import argparse
import concurrent.futures
import os
import random
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SOLVER = os.path.join(ROOT, "latest", "bnbdrv")

PARAMS = ["PARAM/default.par", "PARAM/rule1.par", "PARAM/rule2.par",
          "PARAM/rule3.par", "PARAM/rule4.par", "PARAM/glb.par",
          # plain QP bound (bound type 0), which the PARAM/ files don't cover
          os.path.join(ROOT, "tests", "params", "qp0.par")]
# initial upper bound modes: None = solver default (no incumbent)
UB_MODES = ["none", "opt", "opt+2"]
TOL = 1e-6


# ---------------------------------------------------------------- instances

def sym_matrix(rng, n, lo, hi, zero_prob=0.0):
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = 0 if rng.random() < zero_prob else rng.randint(lo, hi)
            M[i][j] = M[j][i] = v
    return M


def grid_distance(rng, n):
    rows = rng.randint(1, n)
    cols = -(-n // rows)
    cells = [(r, c) for r in range(rows) for c in range(cols)][:n]
    return [[abs(a[0] - b[0]) + abs(a[1] - b[1]) for b in cells] for a in cells]


def make_instance(family, n, rng):
    """Return (A, B, C or None)."""
    if family == "grid":       # nug-like: integer flows, Manhattan grid distances
        return sym_matrix(rng, n, 0, 9, zero_prob=0.3), grid_distance(rng, n), None
    if family == "ties":       # tiny value range: lots of equal-cost permutations
        return sym_matrix(rng, n, 0, 2), sym_matrix(rng, n, 0, 2), None
    if family == "asymB":      # tai-b-like: A symmetric, B asymmetric, zero diagonals
        B = [[0 if i == j else rng.randint(0, 9) for j in range(n)] for i in range(n)]
        return sym_matrix(rng, n, 0, 9), B, None
    if family == "linear":     # symmetric quadratic part plus linear costs
        C = [[rng.randint(0, 20) for _ in range(n)] for _ in range(n)]
        return sym_matrix(rng, n, 0, 9), sym_matrix(rng, n, 0, 9), C
    raise ValueError(family)


FAMILIES = ["grid", "ties", "asymB", "linear"]


def write_dat(path, A, B, C):
    n = len(A)
    with open(path, "w") as f:
        f.write(f"{n}\n\n")
        for M in (A, B) + ((C,) if C else ()):
            for row in M:
                f.write(" ".join(map(str, row)) + "\n")
            f.write("\n")


def brute_force(A, B, C):
    """Minimum of sum_ij A[i][j] B[p(i)][p(j)] + sum_i C[i][p(i)], plus all attained values."""
    n = len(A)
    best = [None]
    values = set()
    p = [0] * n
    used = [False] * n

    def rec(k, cost):
        if k == n:
            values.add(cost)
            if best[0] is None or cost < best[0]:
                best[0] = cost
            return
        for loc in range(n):
            if used[loc]:
                continue
            # cost of placing facility k at loc, against facilities 0..k-1
            d = A[k][k] * B[loc][loc] + (C[k][loc] if C else 0)
            for i in range(k):
                d += A[i][k] * B[p[i]][loc] + A[k][i] * B[loc][p[i]]
            used[loc] = True
            p[k] = loc
            rec(k + 1, cost + d)
            used[loc] = False

    rec(0, 0)
    return best[0], values


# ---------------------------------------------------------------- solver

def run_solver(solver, dat, params, ub):
    cmd = [solver, "-p", params, "-s", "/nonexistent.sym"]
    if ub is not None:
        cmd += ["-u", str(ub)]
    cmd.append(dat)
    try:
        r = subprocess.run(cmd, cwd=os.path.dirname(solver), capture_output=True,
                           text=True, timeout=60)
    except subprocess.TimeoutExpired:
        return None, "timeout"
    text = r.stdout + r.stderr
    if r.returncode != 0 or "---- RESULTS ----" not in text:
        return None, f"crashed (rc={r.returncode})"
    objs = [float(x) for x in re.findall(r"New objective value:\s*([-0-9.eE+]+)", text)]
    m = re.search(r"root bound:\s*([-0-9.eE+]+|nan|-?inf)", text)
    root = float(m.group(1)) if m else None
    return {"objs": objs, "root": root, "cmd": " ".join(cmd)}, None


def check(case, solver):
    """Run one (instance, params, ub mode) combination; return an error string or None."""
    name, dat, opt, values, params, mode = case
    ub = {"none": None, "opt": opt, "opt+2": opt + 2}[mode]
    res, err = run_solver(solver, dat, params, ub)
    where = f"{name} params={os.path.basename(params)} ub={mode} (opt={opt})"
    if err:
        return f"{where}: {err}"
    final = res["objs"][-1] if res["objs"] else None
    if res["root"] is None or res["root"] != res["root"]:
        return f"{where}: no root bound reported\n    {res['cmd']}"
    if res["root"] > opt + TOL:
        return f"{where}: root bound {res['root']} exceeds optimum\n    {res['cmd']}"
    for o in res["objs"]:
        if o < opt - TOL:
            return f"{where}: reported objective {o} below optimum\n    {res['cmd']}"
        if not any(abs(o - v) < TOL for v in values):
            return f"{where}: reported objective {o} is not attained by any permutation\n    {res['cmd']}"
    if mode in ("none", "opt+2") and (final is None or abs(final - opt) > TOL):
        return f"{where}: finished with {final}, expected optimum {opt}\n    {res['cmd']}"
    return None


# ---------------------------------------------------------------- main

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=DEFAULT_SOLVER, help="bnbdrv binary to test")
    p.add_argument("--sizes", default="2-8", help="instance sizes, e.g. 3-7 (default 2-8)")
    p.add_argument("--count", type=int, default=6, help="instances per family and size (default 6)")
    p.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    p.add_argument("--families", nargs="+", default=FAMILIES, choices=FAMILIES)
    p.add_argument("--quick", action="store_true", help="2 instances per family and size")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    p.add_argument("--keep", action="store_true", help="keep generated .dat files")
    p.add_argument("--show", type=int, default=25, help="failures to print (default 25)")
    args = p.parse_args()

    solver = os.path.abspath(args.solver)
    if not os.path.exists(solver):
        sys.exit(f"{solver} not found; build first (see RUNNING.md)")
    lo, _, hi = args.sizes.partition("-")
    sizes = range(int(lo), int(hi or lo) + 1)
    count = 2 if args.quick else args.count
    rng = random.Random(args.seed)

    tmp = tempfile.mkdtemp(prefix="qap-bf-")
    cases = []
    for family in args.families:
        for n in sizes:
            for k in range(count):
                A, B, C = make_instance(family, n, rng)
                name = f"{family}_n{n}_{k}"
                dat = os.path.join(tmp, name + ".dat")
                write_dat(dat, A, B, C)
                opt, values = brute_force(A, B, C)
                for params in PARAMS:
                    for mode in UB_MODES:
                        cases.append((name, dat, opt, values, params, mode))

    n_inst = len(cases) // (len(PARAMS) * len(UB_MODES))
    print(f"{n_inst} instances x {len(PARAMS)} params x {len(UB_MODES)} upper-bound modes "
          f"= {len(cases)} solver runs (seed {args.seed}, files in {tmp})", flush=True)

    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for err in ex.map(lambda c: check(c, solver), cases):
            if err:
                failures.append(err)

    by_family = {}
    for f in failures:
        fam = f.split("_n", 1)[0]
        by_family[fam] = by_family.get(fam, 0) + 1
    for f in failures[:args.show]:
        print("FAIL " + f)
    if len(failures) > args.show:
        print(f"... and {len(failures) - args.show} more (use --show N)")
    if not args.keep and not failures:
        for name in os.listdir(tmp):
            os.remove(os.path.join(tmp, name))
        os.rmdir(tmp)
    summary = ", ".join(f"{k}: {v}" for k, v in sorted(by_family.items()))
    print(f"\n{len(cases) - len(failures)}/{len(cases)} passed"
          + (f"  (failures by family: {summary}; instances kept in {tmp})" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
