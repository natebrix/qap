#!/usr/bin/env python3
"""Check enabled CP and frozen-equality certificates against enumeration."""

import itertools
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
from brute_force import make_instance, write_dat


def fixed_optimum(a, b, c, fixed_i, fixed_j):
    n = len(a)
    values = []
    for free in itertools.permutations(j for j in range(n) if j != fixed_j):
        perm = free[:fixed_i] + (fixed_j,) + free[fixed_i:]
        value = sum(a[i][j] * b[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if c is not None:
            value += sum(c[i][perm[i]] for i in range(n))
        values.append(value)
    return min(values)


def verify():
    solver = ROOT / "experiments/dnn_cuts/dnncuts"
    env = dict(os.environ, DNN_CP_REGION_EQ="1", OPENBLAS_NUM_THREADS="1")
    with_initial_cuts = 0
    with tempfile.TemporaryDirectory() as directory:
        for seed in (1, 2, 3):
            rng = random.Random(seed)
            for family in ("grid", "asymB", "linear"):
                a, b, c = make_instance(family, 7, rng)
                path = Path(directory) / f"{family}.dat"
                write_dat(path, a, b, c)
                fixed_i, fixed_j = seed, 6 - seed
                optimum = fixed_optimum(a, b, c, fixed_i, fixed_j)
                for side in (0, 1):
                    check = seed == 1 and family == "grid" and side == 1
                    command = [str(solver), "-side", str(side), "-prefix", "-x", f"{fixed_i}:{fixed_j}",
                               "-I", "7", "-c", str(optimum + 1), "-cpsep", "10:35:100:20",
                               "-cpiters", "55", "-cpdiag", "20:10", "-cpquiet", "25", "-v"]
                    if check:
                        command.append("-cpcheck")
                    command.append(str(path))
                    result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
                    if check:
                        assert "MISMATCH" not in result.stdout and "failed" not in result.stderr, result.stdout + result.stderr
                        equations = [line for line in result.stdout.splitlines() if line.startswith("cp-region-eq ")]
                        assert len(equations) == 2 and all(line.endswith("  ok") for line in equations), result.stdout
                        print("\n".join(equations))
                    cp = re.search(r"^cp solve  lb ([0-9.e+-]+).* certs ([0-9]+) failures ([0-9]+)  kept ([0-9]+)",
                                   result.stdout, re.MULTILINE)
                    quiet = re.search(r"^cp quiet round  lb ([0-9.e+-]+)", result.stdout, re.MULTILINE)
                    assert cp and quiet, result.stdout + result.stderr
                    lb, quiet_lb = float(cp[1]), float(quiet[1])
                    assert int(cp[2]) > 0 and int(cp[3]) == 0, cp[0]
                    with_initial_cuts += int(cp[4]) > 0
                    assert max(lb, quiet_lb) <= optimum + 1e-6, (family, seed, side, lb, quiet_lb, optimum)
                    print(f"{family:6} seed {seed} side {side}  CP {lb:.7f}  quiet {quiet_lb:.7f}  exact {optimum}")
    assert with_initial_cuts > 0


if __name__ == "__main__":
    verify()
