#!/usr/bin/env python3
"""Check that the parent's DNN state survives a stop and resume (state files
version 6, QAP_DNN_SHIP; latest/state.C).

  python3 tests/ship_test.py [--solver PATH] [--slice SECONDS]

The case is a DNN search with every warm-start mechanism of the wil50
configuration (parent-Y warm start, cut phase, pool inheritance), on nug20.
It is solved straight through, then in CPU slices (-t SLICE -o STATE, then
-r STATE) with QAP_DNN_SHIP=64: every node's DNN trace (bound, primal, status,
iterations, warm start, pool, cut phase) must equal the straight run's, so a
resumed search is the same search. Slices with QAP_DNN_SHIP=32 and 16 and
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
ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
       "QAP_DNN_GAMMA1": "1.0", "QAP_DNN_GAMMA": "1.0", "QAP_DNN_CHEB": "3", "QAP_DNN_TIGHT_LAM": "2",
       "QAP_DNN_REFRESH": "40", "DNN_CHOLQR_PASSES": "1", "QAP_DNN_WARM_Y": "2", "QAP_DNN_EXT": "20",
       "QAP_DNN_CUT_ZONE": "0.005", "QAP_DNN_CUT_GRID": "4:5", "QAP_DNN_CUT_SIDE": "0",
       "QAP_DNN_CUT_FREEZE": "1500", "QAP_DNN_POOL": "1", "QAP_DNN_POOL_CHECK": "20"}
CASE = ("nug20", 2570)
FLAGS = ["-p", "PARAM/nug_deep_xbar9_all.par", "-x", "0:4:1:300:1"]


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
    args = FLAGS + ["-u", str(ub), "-o", state, "-t", str(slice_s), dat]
    while rc == 3:
        traces.append(os.path.join(tmp, f"{tag}.{k}.trace"))
        env = {"QAP_DNN_TRACE": traces[-1]}
        if ship:
            env["QAP_DNN_SHIP"] = str(ship)
        rc, out = run(solver, args, env)
        args = FLAGS + ["-o", state, "-t", str(slice_s), "-r", state, dat]
        k += 1
    sizes = os.path.getsize(state) if os.path.exists(state) else 0
    return rc, k, read_trace(traces), out, sizes


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
            if diff or k < 3:
                failures += 1
                print(f"FAIL {line}; {len(diff)} nodes differ from the straight run"
                      + (" (fewer than 3 slices: lower --slice)" if k < 3 else ""))
                for n in sorted(diff)[:5]:
                    print(f"  [{n}]\n    straight: {ref.get(n)}\n    sliced:   {rows.get(n)}")
            else:
                print(f"ok   {line}: every node identical to the straight run")
        else:
            print(f"     {line}")
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
