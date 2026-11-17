#!/usr/bin/env python3
"""Run feature-on CP certificate checks on small QAPs with enumerated optima."""

import argparse
import array
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import re
import shlex
import subprocess
import tempfile


TRACKS = ("coarse", "scaling", "equalities")
SOURCE_FILES = (
    "latest/dnn.C", "latest/dnn.h", "latest/dnncut.C", "latest/dnncut.h",
    "latest/settings.C", "experiments/dnn_cuts/dnncuts.C",
)
OBJECTS = (
    "qpbnd qpbnd2 dlap qap assign link1 solution branch bnb util linalg fw "
    "estimate mat sa_qap stats console lgb_train state symeig settings dnn "
    "dnncut amx_gemm"
).split()
DRIVER = "experiments/dnn_cuts/dnncuts"
CERTIFICATE = re.compile(r"^cp it +[0-9]+ cert (\S+)", re.MULTILINE)
SUMMARY = re.compile(r"^cp solve  lb (\S+).*?status (\S+)  iters (\d+)", re.MULTILINE)
QUIET = re.compile(r"^cp quiet round  lb (\S+)", re.MULTILINE)
LAZY_EXIT = re.compile(r"^lazycert exit [^:]+: cert (\S+)", re.MULTILINE)


def run(command, **kwargs):
    result = subprocess.run(command, capture_output=True, text=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{shlex.join(map(str, command))}\n{result.stdout}{result.stderr}")
    return result.stdout


def required(condition, details):
    if not condition:
        raise RuntimeError(details)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identify(root):
    return {
        "commit": run(["git", "rev-parse", "HEAD"], cwd=root).strip(),
        "source_sha256": {name: sha256(root / name) for name in SOURCE_FILES},
        "linked_object_sha256": {
            f"{name}.o": sha256(root / "latest" / f"{name}.o")
            for name in ("dnn", "dnncut", "settings")
        },
        "driver_sha256": sha256(root / DRIVER),
    }


def make_driver(track, root, temp, cxx, blas_flags):
    # Use the same driver and solver objects as the experiment. The temporary
    # copy sets policy keys which the driver does not expose on all branches.
    source = (root / f"{DRIVER}.C").read_text()
    main = "int main(int argc, char **argv)\n{"
    iteration = "p.max_iter = cp_iters;"
    required(source.count(main) == source.count(iteration) == 1, "Unexpected driver layout")
    source = '#include "settings.h"\n' + source.replace(
        main, main + '\n  settings_set("dnn.dscale0", "2.5", SRC_FLAG);'
        '\n  settings_set("dnn.lazycert", "1", SRC_FLAG);'
        '\n  settings_set("cut.cp_lazycert", "1", SRC_FLAG);',
    )
    source = source.replace(iteration, iteration + " p.tol = 0.005;")
    driver = temp / f"driver_{track}.C"
    driver.write_text(source)
    exe = temp / f"verify_{track}"
    if not blas_flags:
        dependencies = run(["ldd", root / DRIVER])
        openblas = re.search(r"lib(?:scipy_)?openblas\S* => (\S+)", dependencies)
        if openblas:
            lib = Path(openblas.group(1))
            blas_flags = [str(lib), f"-Wl,-rpath,{lib.parent}"]
        else:
            blas_flags = ["-l:liblapack.so.3", "-l:libblas.so.3"]
    command = [
        cxx, "-O2", "-std=gnu++17", "-Wno-register", "-ffp-contract=off",
        "-DHAVE_DNN", "-I.", driver, *[f"{name}.o" for name in OBJECTS],
        "-o", exe, "-lm", *blas_flags,
    ]
    run(command, cwd=root / "latest")
    return exe, {"command": shlex.join(map(str, command)), "temporary_driver_sha256": sha256(driver)}


def fixtures(temp):
    n = 7
    for sample in range(4):
        rng = random.Random(10519 + sample)
        a = [[rng.randint(-20, 31) for _ in range(n)] for _ in range(n)]
        b = [[rng.randint(-9, 43) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            a[2][i] = a[i][2] = 0
        fixes = {0: 4} if sample % 2 else {}
        instance = temp / f"random_{sample}.dat"
        instance.write_text(f"{n}\n" + "\n".join(" ".join(map(str, row)) for row in a + b) + "\n")
        optimum = min(
            sum(a[i][j] * b[p[i]][p[j]] for i in range(n) for j in range(n))
            for p in itertools.permutations(range(n))
            if all(p[i] == j for i, j in fixes.items())
        )
        yield sample, instance, fixes, n - len(fixes), optimum


def coarse_pool(temp, sample, m, side):
    # Each signed pair encodes an exact permutation identity. The first
    # rows sum all coarse pair masses (RHS ±1), and the remaining rows
    # compare a cross-label marginal with a diagonal label mass (RHS 0).
    labels, n_labels, rows = [a % 3 for a in range(m)], 3, []
    for i, j in ((0, 1), (2, 3), (4, 5)):
        for sign in (1, -1):
            prefix = [2, i, j]
            rows.append(prefix + [0] * 6 + [sign * 1024] * 9)
            for g in range(n_labels):
                diagonal = [0] * 6
                diagonal[g] = -sign * 1024
                pair = [sign * 1024 if a == g else 0
                        for a in range(n_labels) for _ in range(n_labels)]
                rows.append(prefix + diagonal + pair)
    path = temp / f"pool_{sample}_{side}.txt"
    with path.open("w") as out:
        print("QAP_COARSE_POOL", 1, side, m, n_labels, len(rows), file=out)
        print(*labels, file=out)
        for row in rows:
            print(*row, file=out)
    return path


def verify_case(track, exe, temp, fixture, side):
    sample, instance, fixes, m, optimum = fixture
    saved = temp / f"state_{track}"
    command = [
        exe, "-prefix", "-side", str(side), "-I", "25", "-k", "20", "-cpsep",
        "10:10:100:20", "-cpiters", "500", "-cpquiet", "12", "-c",
        str(optimum + 1000), "-dumpY", saved,
    ]
    environment = {}
    if track == "scaling":
        environment = {"DNN_CP_DSCALE": ".5", "DNN_CP_DSCALE_UNTIL": "22" if fixes else "0"}
    elif track == "equalities":
        environment = {"DNN_CP_REGION_EQ": "1"}
    else:
        command += ["-coarse", coarse_pool(temp, sample, m, side)]
    if fixes:
        command += ["-x", ",".join(f"{i}:{j}" for i, j in fixes.items())]
    command += [instance]
    log = run(command, env={**os.environ, **environment})
    summary = SUMMARY.search(log)
    quiet = QUIET.search(log)
    required(summary and quiet, log)
    lower, status, iterations = summary.groups()
    cp_checks = list(CERTIFICATE.finditer(log))
    required(cp_checks, log)
    after_cp = log[cp_checks[-1].end():]
    before_quiet_checks = re.split(r"\nit +\d+", after_cp, maxsplit=1)[0]
    values = [float(match.group(1)) for match in cp_checks]
    recovery = [float(value) for value in LAZY_EXIT.findall(before_quiet_checks)]
    values += [float(lower), float(quiet.group(1))]
    precise = [value for value in values if math.isfinite(value)]
    finite = precise + [value for value in recovery if math.isfinite(value)]
    required(status not in ("error", "fathomed") and int(iterations) > 10, log)
    required(precise and max(finite) <= optimum + 1e-8, log)

    activation = {}
    if track == "scaling":
        start = "dscale p 0.5 d0 2.5 on 1"
        short_window = "dscale off at iteration 23" in log
        required(start in log and short_window == bool(fixes), log)
        activation = {"scale_on": True, "transition_to_plain": short_window}
    elif track == "equalities":
        pattern = re.escape("eq max |C-s|/||C|| ") + r"(\S+)" + re.escape(" |G|/||G|| ") + r"(\S+)"
        residuals = [(float(c), float(g)) for c, g in re.findall(pattern, log)]
        required(any(c > 1e-6 and g > 1e-6 for c, g in residuals), log)
        activation = {"max_cardinality_residual": max(c for c, _ in residuals),
                      "max_row_residual": max(g for _, g in residuals)}
    else:
        active = [int(value) for value in re.findall(r"coarse \d+ active (\d+)", log)]
        required(active and max(active) > 0, log)
        activation = {"max_active_coarse_rows": max(active)}

    with saved.with_suffix(".Y.bin").open("rb") as stream:
        y = array.array("d")
        y.fromfile(stream, (1 + m * m) ** 2)
    required(
        all(map(math.isfinite, y)) and abs(y[0] - 1) < 1e-12
        and min(y) >= -1e-12 and max(y) <= 1 + 1e-12,
        "Saved Y is outside the original lifted box",
    )
    return {
        "sample": sample, "side": side, "free_order": m, "exact_optimum": optimum,
        "status": status, "iterations": int(iterations), "finite_values": len(finite),
        "best_value": max(finite), "best_value_minus_optimum": max(finite) - optimum,
        "best_precise_value_minus_optimum": max(precise) - optimum,
        "quiet_value": float(quiet.group(1)),
        "lazy_checks": len(re.findall(r", lazy\)[^\n]*\ncp it ", log)),
        "lazy_exit_values": recovery, "activation": activation,
        "fixture_sha256": sha256(instance), "overridden_environment": environment,
        "command": shlex.join(map(str, command)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for track in TRACKS:
        parser.add_argument(f"--{track}", type=Path, help=f"built {track} checkout")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--blas-link-flags", help="override link flags inferred from the built driver")
    parser.add_argument("--out", type=Path, help="also save results to this JSON file")
    args = parser.parse_args()
    if not any(getattr(args, t) for t in TRACKS):
        parser.error("specify at least one experiment checkout")
    results = {}
    with tempfile.TemporaryDirectory(prefix="verify-qap-cp-") as folder:
        temp = Path(folder)
        inputs = list(fixtures(temp))
        for track in TRACKS:
            path = getattr(args, track)
            if path is None:
                continue
            root = path.resolve()
            build = identify(root)
            exe, compile_info = make_driver(
                track, root, temp, args.cxx, shlex.split(args.blas_link_flags or ""),
            )
            cases = [verify_case(track, exe, temp, fixture, side)
                     for fixture in inputs for side in (0, 1)]
            lazy_exits = sum(bool(case["lazy_exit_values"]) for case in cases)
            required(lazy_exits > 0, f"{track} did not exercise lazy certificate recovery")
            results[track] = {
                "build": build | compile_info, "cases": cases,
                "case_count": len(cases), "lazy_exit_cases": lazy_exits,
                "lazy_checks": sum(case["lazy_checks"] for case in cases),
                "finite_values": sum(case["finite_values"] for case in cases),
                "worst_value_minus_optimum": max(case["best_value_minus_optimum"] for case in cases),
                "worst_precise_value_minus_optimum": max(
                    case["best_precise_value_minus_optimum"] for case in cases
                ),
            }
    result = json.dumps(results, indent=2, allow_nan=False) + "\n"
    if args.out:
        args.out.write_text(result)
    print(result, end="")


if __name__ == "__main__":
    main()
