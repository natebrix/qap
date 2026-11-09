#!/usr/bin/env python3
"""Check that the parent's DNN state survives a stop and resume (state files
version 6, ship.bits; latest/state.C).

  python3 tests/ship_test.py [--solver PATH] [--slice SECONDS]

The case is a DNN search with every warm-start mechanism of the wil50
configuration (parent-Y warm start, cut phase, pool inheritance), on nug20.
It is solved straight through, then in CPU slices (-t SLICE -o STATE, then
-r STATE) with ship.bits = 64: every node's DNN trace (bound, primal, status,
iterations, warm start, pool, cut phase) must equal the straight run's, so a
resumed search is the same search, and so must the per-level DNN, cut-phase and
fathom-cause tables the last slice prints (state version 7 carries them
through every resume; CPU columns aside). Slices with ship.bits = 32 and 16 and
without shipping are reported for comparison: their warm-started node count
and DNN iterations.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import envpar  # noqa: E402  (a par file per ship.bits value)
# the solver's settings are the [settings] section of tests/par/ship.par (the
# DNN bound at depths 0-4, the cut phase, the pool); only the thread variables
# and the pool diagnostic stay in the environment
ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
       "QAP_DNN_POOL_CHECK": "20"}
CASE = ("nug20", 2570)
PAR = os.path.join(ROOT, "tests", "par", "ship.par")
FLAGS = ["-p", PAR]


def node_key(line):
    return " ".join(sorted(re.search(r"fixes([\d: ]*)\|", line).group(1).split()))


def node_fields(line):
    # everything but the pid and the CPU times
    line = re.sub(r"^pid \d+ ", "", line)
    return re.sub(r" cpu \S+", "", line).strip()


def read_trace(paths):
    rows = {}
    for p in paths:
        if os.path.exists(p):
            for line in open(p):
                if " fixes " in line:
                    rows[node_key(line)] = node_fields(line)
    return rows


def run(solver, args, env):
    e = dict(os.environ)
    e.update(ENV)
    e.update(env)
    r = subprocess.run([solver] + args, cwd=os.path.dirname(solver), capture_output=True, text=True,
                       timeout=3600, env=e)
    return r.returncode, r.stdout + r.stderr


def sliced(solver, dat, ub, tmp, tag, ship, slice_s):
    state = os.path.join(tmp, f"{tag}.state")
    traces, k, rc = [], 0, 3
    flags = FLAGS
    if ship:   # the case's par file with ship.bits set
        par = os.path.join(tmp, f"ship{ship}.par")
        envpar.write_par(PAR, {"ship.bits": str(ship)}, par, solver=solver)
        flags = ["-p", par]
    args = flags + ["-u", str(ub), "-o", state, "-t", str(slice_s), dat]
    while rc == 3:
        traces.append(os.path.join(tmp, f"{tag}.{k}.trace"))
        rc, out = run(solver, args, {"QAP_DNN_TRACE": traces[-1]})
        args = flags + ["-o", state, "-t", str(slice_s), "-r", state, dat]
        k += 1
    sizes = os.path.getsize(state) if os.path.exists(state) else 0
    return rc, k, read_trace(traces), out, sizes


def level_tables(out):
    """the DNN, cut-phase and fathoms-by-cause tables of a run's output,
    without their CPU columns (the last column of the first two)"""
    tables = {}
    for title in ("DNN bound by level", "Cut phase by level", "Fathoms by cause"):
        m = re.search(title + r":\n[^\n]*\n((?: ?\d+ [^\n]*\n)*)", out)
        if m:
            rows = [l.split() for l in m.group(1).splitlines()]
            tables[title] = [r[:-1] if title != "Fathoms by cause" else r for r in rows]
    return tables


def totals(rows):
    it = sum(int(re.search(r"iters (\d+)", r).group(1)) for r in rows.values())
    cut = sum(int(m.group(1)) for r in rows.values() for m in [re.search(r"rounds \d+ iters (\d+)", r)] if m)
    warm = sum(1 for r in rows.values() if " warm 1 " in r)
    return it + cut, warm


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    p.add_argument("--slice", type=float, default=3.0, help="CPU seconds per slice")
    a = p.parse_args()
    solver = os.path.abspath(a.solver)
    tmp = tempfile.mkdtemp(prefix="qap-ship-")
    name, ub = CASE
    dat = os.path.join(ROOT, "qaplib", name + ".dat")

    straight_trace = os.path.join(tmp, "straight.trace")
    rc, out = run(solver, FLAGS + ["-u", str(ub), dat], {"QAP_DNN_TRACE": straight_trace})
    if rc != 0:
        print(f"FAIL {name}: straight run exited {rc}\n{out[-2000:]}")
        return 1
    ref = read_trace([straight_trace])
    ref_tables = level_tables(out)
    if "DNN bound by level" not in ref_tables or "Fathoms by cause" not in ref_tables:
        print(f"FAIL {name}: the straight run printed no DNN or fathom-cause table")
        return 1
    it0, w0 = totals(ref)
    print(f"straight: {len(ref)} DNN nodes, {w0} warm-started, {it0:,} DNN iterations")

    failures = 0
    for tag, ship in (("ship64", 64), ("ship32", 32), ("ship16", 16), ("noship", 0)):
        rc, k, rows, out, size = sliced(solver, dat, ub, tmp, tag, ship, a.slice)
        if rc != 0:
            print(f"FAIL {tag}: slice {k} exited {rc}\n{out[-2000:]}")
            failures += 1
            continue
        it, w = totals(rows)
        line = f"{tag}: {k} slices, {len(rows)} DNN nodes, {w} warm-started, {it:,} DNN iterations"
        if ship == 64:
            diff = [n for n in set(ref) | set(rows) if ref.get(n) != rows.get(n)]
            tables = level_tables(out)
            bad_tables = [t for t in set(ref_tables) | set(tables) if ref_tables.get(t) != tables.get(t)]
            if bad_tables:
                failures += 1
                print(f"FAIL {tag}: the last slice's {', '.join(sorted(bad_tables))} differ from the straight run's")
                for t in sorted(bad_tables):
                    print(f"  {t}:\n    straight: {ref_tables.get(t)}\n    sliced:   {tables.get(t)}")
            if diff or k < 3:
                failures += 1
                print(f"FAIL {line}; {len(diff)} nodes differ from the straight run"
                      + (" (fewer than 3 slices: lower --slice)" if k < 3 else ""))
                for n in sorted(diff)[:5]:
                    print(f"  [{n}]\n    straight: {ref.get(n)}\n    sliced:   {rows.get(n)}")
            else:
                print(f"ok   {line}: every node identical to the straight run"
                      + ("" if bad_tables else f", and the {len(tables)} level tables"))
        else:
            print(f"     {line}")
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
