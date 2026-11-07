#!/usr/bin/env python3
"""Validity of the CP cut phase in the search (bnbdrv, QAP_DNN_CUT_CP), by brute
force: small instances searched with phases forced (QAP_DNN_CUT_ZONE=1: every DNN
miss runs one), pools on (QAP_DNN_POOL=1, QAP_DNN_POOL_CHECK=20), the per-node trace
(QAP_DNN_TRACE). The upper bound is the optimum + 2, so the search has to find
the optimum itself: a bound above a node's optimum would show as a wrong final
objective or a node bound above that node's optimum.

Checks, for every arm and instance:
- the final objective equals the optimum (enumeration);
- every traced node's first-solve bound and phase bound are at most that node's
  optimum (by enumeration when at most 9 elements are free, else by a plain QPB
  search from the node, bnbdrv -F), within the trace's 4 decimals (1e-4).

Instances: the four families of tests/brute_force.py at n = 8 and 9 (prefix
regions), and nug12 with its matrices swapped so that the 3 x 4 grid is the first
matrix (grid 3:4 regions, side 0). Arms:
  cp       QAP_DNN_CUT_CP=50:100 (Y restarted at 0)
  cpwarmy  the same with QAP_DNN_CUT_WARMY=1 (the phase starts from the first solve's Y)
  poolcp   cpwarmy plus QAP_DNN_POOL_CP=1 (a pooled child's first solve in CP form)
  cplazy   poolcp plus the lazy certificate, in the plain solves (DNN_CORE_LAZYCERT=1)
           and in the CP ones (DNN_CP_LAZYCERT=1)

  python3 validate_search.py [--core 0] [--count 2] [--out DIR] [--arms cp cpwarmy poolcp cplazy]
"""
import argparse
import importlib.util
import itertools
import os
import random
import re
import subprocess
import sys

import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(D, "..", "..", "..", "..", ".."))
BNB = os.path.join(ROOT, "latest", "bnbdrv")
BASE_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "QAP_DNN_GAMMA1": "1.0", "QAP_DNN_GAMMA": "1.0",
            "QAP_DNN_CHEB": "3", "QAP_DNN_WARM_Y": "2", "QAP_DNN_CUT_ZONE": "1.0", "QAP_DNN_CUT_SIDE": "0",
            "QAP_DNN_CUT_K": "100", "QAP_DNN_CUT_ROUNDS": "4", "QAP_DNN_CUT_ITERS": "100", "QAP_DNN_CUT_FREEZE": "200",
            "QAP_DNN_POOL": "1", "QAP_DNN_POOL_CHECK": "20", "QAP_DNN_CUT_CP": "50:100"}
ARMS = {"cp": {}, "cpwarmy": {"QAP_DNN_CUT_WARMY": "1"},
        "poolcp": {"QAP_DNN_CUT_WARMY": "1", "QAP_DNN_POOL_CP": "1"},
        "cplazy": {"QAP_DNN_CUT_WARMY": "1", "QAP_DNN_POOL_CP": "1", "DNN_CORE_LAZYCERT": "1", "DNN_CP_LAZYCERT": "1"}}
FLAGS = ["-p", "PARAM/nug_deep_xbar9_all.par", "-s", "/nonexistent.sym", "-x", "0:30:1:0:1"]
TOL = 1e-4


def load_bf():
    spec = importlib.util.spec_from_file_location("brute_force", os.path.join(ROOT, "tests", "brute_force.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def objective(A, B, C, p):
    v = float((A * B[np.ix_(p, p)]).sum())
    if C is not None:
        v += float(C[np.arange(len(p)), p].sum())
    return v


def node_opt_enum(A, B, C, fixes):
    n = len(A)
    fac = [i for i, _ in fixes]; loc = [j for _, j in fixes]
    ff = [i for i in range(n) if i not in fac]; fl = [j for j in range(n) if j not in loc]
    best = None
    for perm in itertools.permutations(fl):
        p = np.zeros(n, dtype=int)
        for i, j in fixes: p[i] = j
        for i, j in zip(ff, perm): p[i] = j
        v = objective(A, B, C, p)
        best = v if best is None else min(best, v)
    return best


def run(core, args, env, timeout=1800):
    e = dict(os.environ); e.update(BASE_ENV); e.update(env)
    r = subprocess.run(["taskset", "-c", str(core), BNB] + args, cwd=os.path.dirname(BNB), capture_output=True,
                       text=True, timeout=timeout, env=e)
    return r.returncode, r.stdout + r.stderr


def node_opt_search(core, dat, fixes):
    """the node's optimum by a plain QPB search from it (no DNN)"""
    x = ",".join("%d:%d" % f for f in fixes)
    rc, out = run(core, ["-p", "PARAM/nug_deep_xbar9_all.par", "-s", "/nonexistent.sym", "-F", x, dat],
                  {"QAP_DNN_CUT_CP": "", "QAP_DNN_POOL": "0", "QAP_DNN_CUT_ZONE": "0"})
    objs = [float(v) for v in re.findall(r"New objective value:\s*([-0-9.eE+]+)", out)]
    return objs[-1] if objs else None


def read_trace(path):
    rows = []
    if not os.path.exists(path):
        return rows
    for line in open(path):
        if " fixes " not in line:
            continue
        fx = re.search(r"fixes([\d: ]*)\|", line).group(1).split()
        fixes = tuple(sorted(tuple(int(v) for v in f.split(":")) for f in fx))
        lb = float(re.search(r"first lb (\S+)", line).group(1))
        m = re.search(r"\| cut 1 best (\S+) rounds (\d+) iters (\d+)", line)
        rows.append(dict(fixes=fixes, lb=lb, cut=bool(m), best=float(m.group(1)) if m else None,
                         seps=int(m.group(2)) if m else 0, cut_iters=int(m.group(3)) if m else 0,
                         cp=" cp 1 " in line, pcp=" pcp 1 " in line, pool=int(re.search(r" pool (\d)", line).group(1)),
                         pool_out=int(re.search(r"pool_out (\d+)", line).group(1))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", type=int, default=0)
    ap.add_argument("--count", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(D, "logs", "validate_search"))
    ap.add_argument("--arms", nargs="+", default=["cp", "cpwarmy"])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bf = load_bf()
    rng = random.Random(1107)
    insts = []
    for family in bf.FAMILIES:
        for n in (8, 9):
            for k in range(a.count):
                A, B, C = bf.make_instance(family, n, rng)
                dat = os.path.join(a.out, "%s_n%d_%d.dat" % (family, n, k))
                bf.write_dat(dat, A, B, C)
                insts.append((os.path.basename(dat)[:-4], dat, np.array(A, float), np.array(B, float),
                              None if C is None else np.array(C, float), {}))
    lines = open(os.path.join(ROOT, "qaplib", "nug12.dat")).read().split()
    n = int(lines[0]); nums = list(map(int, lines[1:1 + 2 * n * n]))
    F = [nums[i * n:(i + 1) * n] for i in range(n)]; Dm = [nums[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]
    dat = os.path.join(a.out, "nug12_swapped.dat")
    bf.write_dat(dat, Dm, F, None)
    insts.append(("nug12_swapped", dat, np.array(Dm, float), np.array(F, float), None, {"QAP_DNN_CUT_GRID": "3:4"}))
    failures = []
    stats = {arm: dict(runs=0, nodes=0, phases=0, cp_phases=0, seps=0, pooled=0, pcp=0, pool_out=0, worst=-np.inf,
                       checked=0) for arm in a.arms}
    opt_cache = {}
    for name, dat, A, B, C, env in insts:
        nn = len(A)
        if (name, ()) not in opt_cache:
            opt_cache[(name, ())] = 578.0 if name == "nug12_swapped" else \
                min(objective(A, B, C, np.array(p)) for p in itertools.permutations(range(nn)))
        opt = opt_cache[(name, ())]
        for arm in a.arms:
            trace = os.path.join(a.out, "%s.%s.trace" % (name, arm))
            if os.path.exists(trace): os.unlink(trace)
            e = dict(env); e.update(ARMS[arm]); e["QAP_DNN_TRACE"] = trace
            rc, out = run(a.core, FLAGS + ["-u", repr(opt + 2), dat], e)
            open(os.path.join(a.out, "%s.%s.out" % (name, arm)), "w").write(out)
            st = stats[arm]; st["runs"] += 1
            objs = [float(v) for v in re.findall(r"New objective value:\s*([-0-9.eE+]+)", out)]
            if rc != 0 or not objs or abs(objs[-1] - opt) > 1e-6:
                failures.append("%s %s: rc %d, final %s, optimum %g" % (name, arm, rc, objs[-1] if objs else None, opt))
            for r in read_trace(trace):
                st["nodes"] += 1; st["phases"] += r["cut"]; st["cp_phases"] += r["cut"] and r["cp"]
                st["seps"] += r["seps"]; st["pooled"] += r["pool"]; st["pcp"] += r["pcp"]; st["pool_out"] += r["pool_out"] > 0
                key = (name, r["fixes"])
                if key not in opt_cache:
                    free = nn - len(r["fixes"])
                    opt_cache[key] = node_opt_enum(A, B, C, list(r["fixes"])) if free <= 9 else \
                        node_opt_search(a.core, dat, list(r["fixes"]))
                no = opt_cache[key]
                if no is None:
                    failures.append("%s %s node %s: no node optimum" % (name, arm, r["fixes"])); continue
                st["checked"] += 1
                for what, v in (("first", r["lb"]), ("phase", r["best"])):
                    if v is None: continue
                    st["worst"] = max(st["worst"], v - no)
                    if v > no + TOL:
                        failures.append("%s %s node %s: %s bound %.4f > node optimum %.4f" % (name, arm, r["fixes"], what, v, no))
    print("CP cut phase in the search: %d instances (n = 8, 9 and nug12), ub = optimum + 2" % len(insts))
    for arm, st in stats.items():
        print("  %-8s %d runs, %d traced DNN nodes (%d checked), %d phases (%d in CP form, %d separations), %d pooled "
              "nodes (%d first solves in CP form), %d nodes handed on a pool; max(bound - node optimum) = %.4g"
              % (arm, st["runs"], st["nodes"], st["checked"], st["phases"], st["cp_phases"], st["seps"], st["pooled"],
                 st["pcp"], st["pool_out"], st["worst"]))
    if failures:
        print("FAILURES (%d):" % len(failures))
        for f in failures[:40]: print("  " + f)
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
