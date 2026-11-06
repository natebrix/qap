#!/usr/bin/env python3
"""Run two bnbdrv binaries on every case of tests/brute_force.py (same
instances, params files and upper-bound modes for a seed) and require the
whole output to be identical, apart from CPU times: incumbents, node /
fathomed / enumerated counts, root bound, Frank-Wolfe iterations, the
per-level table, strategy breakdown and gap statistics. brute_force.py only
checks the results against enumeration; this checks that the search itself
did not change, on tiny instances (n = 2 to 8) where sizes 0 and 1 occur.

  python3 bf_compare.py OLD_BNBDRV NEW_BNBDRV [--seed 1] [--jobs 2] [--sizes 2-8]
                        [--solver-args "-x 0:20"]
"""
import argparse
import concurrent.futures
import os
import random
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import brute_force as bf  # noqa: E402

ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
EXTRA = []   # --solver-args


def output(solver, dat, params, ub):
    cmd = [solver, "-p", params, "-s", "/nonexistent.sym"] + EXTRA
    if ub is not None:
        cmd += ["-u", str(ub)]
    cmd.append(dat)
    r = subprocess.run(cmd, cwd=os.path.dirname(solver), capture_output=True, text=True,
                       timeout=120, env=ENV)
    lines = []
    for line in (r.stdout + r.stderr).splitlines():
        if "Time" in line or line.startswith("time "):
            continue                                   # CPU times
        if re.match(r"^\s*\d+\s+\d+\t", line):          # per-level table: drop the time column
            line = line.rsplit(None, 1)[0]
        line = re.sub(r"\b[0-9.]+(e-?\d+)? s$", "<t> s", line)  # "DNN bound: ..., 0.0121 s"
        lines.append(line)
    return r.returncode, "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("old")
    p.add_argument("new")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--sizes", default="2-8")
    p.add_argument("--count", type=int, default=6)
    p.add_argument("--solver-args", default="", help="extra solver flags, e.g. '-x 0:20'")
    a = p.parse_args()
    EXTRA.extend(a.solver_args.split())
    old, new = os.path.abspath(a.old), os.path.abspath(a.new)
    lo, _, hi = a.sizes.partition("-")
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp(prefix="qap-bfcmp-")
    cases = []
    for family in bf.FAMILIES:
        for n in range(int(lo), int(hi or lo) + 1):
            for k in range(a.count):
                A, B, C = bf.make_instance(family, n, rng)
                dat = os.path.join(tmp, f"{family}_n{n}_{k}.dat")
                bf.write_dat(dat, A, B, C)
                opt, _ = bf.brute_force(A, B, C)
                for params in bf.PARAMS:
                    for mode in bf.UB_MODES:
                        ub = {"none": None, "opt": opt, "opt+2": opt + 2}[mode]
                        cases.append((dat, params, ub))

    def one(case):
        ra, oa = output(old, *case)
        rb, ob = output(new, *case)
        return case, ra == rb and oa == ob, ra, rb

    bad = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for case, same, ra, rb in ex.map(one, cases):
            if not same:
                bad += 1
                if bad <= 10:
                    print(f"DIFFERENT: {case} rc {ra} / {rb}")
    print(f"seed {a.seed}{' ' + a.solver_args if a.solver_args else ''}: {len(cases)} runs, {len(cases) - bad} with identical output, {bad} different")
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
