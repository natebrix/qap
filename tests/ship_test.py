#!/usr/bin/env python3
"""Check that the parent's DNN state survives a stop and resume (state files
version 6, QAP_DNN_SHIP; latest/state.C).

  python3 tests/ship_test.py [--solver PATH] [--slice SECONDS] [--basis]

The case is a DNN search with every warm-start mechanism of the wil50
configuration (parent-Y warm start, cut phase, pool inheritance), on nug20.
It is solved straight through, then in CPU slices (-t SLICE -o STATE, then
-r STATE) with QAP_DNN_SHIP=64: every node's DNN trace (bound, primal, status,
iterations, warm start, pool, cut phase) must equal the straight run's, so a
resumed search is the same search. Slices with QAP_DNN_SHIP=32 and 16 and
without shipping are reported for comparison: their warm-started node count
and DNN iterations.

--basis: every run with DNN_SHIP_BASIS=1, so state files are version 7 and carry
the parents' tracked eigenbases (latest/state.C, dnn.C's seed_block); the sliced
run with QAP_DNN_SHIP=64 (the basis in fp32, its own precision) must again equal
the straight run node for node, i.e. a child whose parent's basis came from a
state file solves as one whose parent was bounded in the same process. The
nodes whose first solve was seeded (" seed 1 res R" in the trace) are counted, and
the basis tokens in the state files measured. A bf16 basis
(QAP_DNN_SHIP_BASIS_BITS=16) is reported for comparison.
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


def basis_tokens(state):
    """(version, [(m, len of the dnn token, len of the basis token)]) of a state file's warm lines"""
    if not os.path.exists(state):
        return 0, []
    with open(state) as f:
        version = int(f.readline().split()[1])
        rows = []
        for line in f:
            if line.startswith("warm "):
                t = line.split()
                if version >= 7:
                    rows.append((int(t[1]), len(t[-2]) if t[-2] != "-" else 0, len(t[-1]) if t[-1] != "-" else 0))
    return version, rows


def sliced(solver, dat, ub, tmp, tag, ship, slice_s, extra_env=None):
    state = os.path.join(tmp, f"{tag}.state")
    traces, k, rc = [], 0, 3
    args = FLAGS + ["-u", str(ub), "-o", state, "-t", str(slice_s), dat]
    toks = []
    while rc == 3:
        traces.append(os.path.join(tmp, f"{tag}.{k}.trace"))
        env = {"QAP_DNN_TRACE": traces[-1]}
        env.update(extra_env or {})
        if ship:
            env["QAP_DNN_SHIP"] = str(ship)
        rc, out = run(solver, args, env)
        if rc == 3:
            toks.append(basis_tokens(state))
        args = FLAGS + ["-o", state, "-t", str(slice_s), "-r", state, dat]
        k += 1
    sizes = os.path.getsize(state) if os.path.exists(state) else 0
    return rc, k, read_trace(traces), out, sizes, toks


def totals(rows):
    it = sum(int(re.search(r"iters (\d+)", r).group(1)) for r in rows.values())
    cut = sum(int(m.group(1)) for r in rows.values() for m in [re.search(r"rounds \d+ iters (\d+)", r)] if m)
    warm = sum(1 for r in rows.values() if " warm 1 " in r)
    return it + cut, warm


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    p.add_argument("--slice", type=float, default=3.0, help="CPU seconds per slice")
    p.add_argument("--basis", action="store_true", help="every run with DNN_SHIP_BASIS=1 (state version 7)")
    a = p.parse_args()
    if a.basis:
        ENV["DNN_SHIP_BASIS"] = "1"
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
    # bnb.C ends a node's trace line with " seed S res R" when DNN_SHIP_BASIS=1
    seeded = lambda rows: sum(1 for r in rows.values() if re.search(r" seed 1 res \S+$", r))
    print(f"straight: {len(ref)} DNN nodes, {w0} warm-started, {it0:,} DNN iterations"
          + (f", {seeded(ref)} seeded from the parent's basis" if a.basis else ""))
    if a.basis and seeded(ref) == 0:
        print("FAIL: no node was seeded from a basis in the straight run")
        return 1

    failures = 0
    cases = [("ship64", 64, {}), ("ship32", 32, {}), ("ship16", 16, {}), ("noship", 0, {})]
    if a.basis:
        cases.insert(1, ("ship64b16", 64, {"QAP_DNN_SHIP_BASIS_BITS": "16"}))
    for tag, ship, xenv in cases:
        rc, k, rows, out, size, toks = sliced(solver, dat, ub, tmp, tag, ship, a.slice, xenv)
        if rc != 0:
            print(f"FAIL {tag}: slice {k} exited {rc}\n{out[-2000:]}")
            failures += 1
            continue
        it, w = totals(rows)
        line = f"{tag}: {k} slices, {len(rows)} DNN nodes, {w} warm-started, {it:,} DNN iterations"
        if a.basis:
            line += f", {seeded(rows)} seeded"
            vers = sorted({v for v, _ in toks})
            bas = [(m, d, b) for _, r in toks for (m, d, b) in r]
            nb = sum(1 for _, _, b in bas if b)
            if bas:
                line += (f"; state version {'/'.join(map(str, vers))}, {nb} of {len(bas)} warm lines with a basis"
                         + (f" (basis token mean {sum(b for _, _, b in bas if b) / max(1, nb) / 1e3:.0f} kB,"
                            f" DNN state token mean {sum(d for _, d, b in bas if b) / max(1, nb) / 1e3:.0f} kB)" if nb else ""))
            if ship and (vers != [7] or nb == 0) and k > 1:
                failures += 1
                print(f"FAIL {tag}: the state files carry no basis (versions {vers}, {nb} bases)")
        # the exact case: the DNN state in doubles and the basis in fp32 (its own
        # precision); a bf16 basis (ship64b16) is reported only, like ship32 and ship16
        if ship == 64 and not xenv:
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
