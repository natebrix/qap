#!/usr/bin/env python3
"""Validity of the CP-ADMM certificate on small instances, by brute force
(as tests/dnn_test.py part A does for the plain bound).

For each instance and node: the node's optimum by enumeration; 3 Polyak rounds
(dnncuts, cutoff = the optimum) build a regional triangle pool, written whole
after round 3 (-savepoolall); then dnncuts -cpadmm runs the uncut round and the
CP-ADMM solve with that pool frozen in the Y-step, in three arms:
  base     1 sweep, 1 pass, cutoff = optimum (-cpcheck: the extraction self-test)
  hard     3 sweeps of 5 passes (multipliers pushed harder), cutoff = optimum
  fathom   1 sweep, cutoff = optimum - 1: a fathom must come with lb > cutoff
Every certified value of every check (the "cert" field of each "it" line) must
be <= optimum + 1e-9 |optimum|, no run may report "fathomed" with the cutoff at
the optimum, and every self-test line must end in "ok".

Instances: the four families of tests/brute_force.py at n = 7, 8, 9 (prefix
regions of the free elements) and nug12 with its matrices swapped so that the
3 x 4 grid distance is the first matrix (grid 3:4 regions, side 0): its root
(optimum 578) and nodes with three facilities fixed (m = 9, enumerated).

  python3 validate_cp.py [--core 2] [--count 2] [--out logs/validate]
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
DNNCUTS = os.path.join(ROOT, "experiments", "dnn_cuts", "dnncuts")


def load_bf():
    spec = importlib.util.spec_from_file_location("brute_force", os.path.join(ROOT, "tests", "brute_force.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def values(A, B, C):
    n = len(A)
    P = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    A = np.asarray(A, float); B = np.asarray(B, float)
    v = np.einsum("ik,pik->p", A, B[P[:, :, None], P[:, None, :]])
    if C is not None:
        v = v + np.asarray(C, float)[np.arange(n), P].sum(axis=1)
    return P, v


def node_opt_fixed(A, B, C, fixes):
    """optimum with facility i at location j for (i, j) in fixes, by enumerating the rest"""
    n = len(A)
    A = np.asarray(A, float); B = np.asarray(B, float)
    fac = [i for i, _ in fixes]; loc = [j for _, j in fixes]
    free_f = [i for i in range(n) if i not in fac]
    free_l = [j for j in range(n) if j not in loc]
    best = None
    for perm in itertools.permutations(free_l):
        p = [0] * n
        for i, j in fixes: p[i] = j
        for i, j in zip(free_f, perm): p[i] = j
        p = np.array(p)
        v = float((A * B[np.ix_(p, p)]).sum())
        if C is not None: v += float(sum(C[i][p[i]] for i in range(n)))
        best = v if best is None else min(best, v)
    return best


def run(args, log):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    out = subprocess.run(["taskset", "-c", str(args_core)] + args, capture_output=True, text=True, env=env, timeout=900)
    with open(log, "w") as f: f.write(out.stdout + out.stderr)
    return out.returncode, out.stdout


def check_node(dat, fixes, opt, regions, tag, outdir, failures, stats):
    x = ",".join("%d:%d" % fx for fx in fixes) if fixes else "-"
    pool = os.path.join(outdir, tag + ".pool")
    rc, _ = run([DNNCUTS] + regions + ["-c", repr(opt), "-T", "1", "-k", "150", "-v", "-I", "200", "-i", "200",
                 "-r", "3", "-delta", "20", "-nocap", "-savepoolall", "3:" + pool, "-x", x, dat],
                os.path.join(outdir, tag + ".rounds.log"))
    if rc != 0:
        failures.append("%s: rounds failed rc=%d" % (tag, rc)); return
    npool = 0
    with open(pool) as f:
        tok = f.read().strip()
        if tok != "-": npool = int(tok.split("/")[4])
    stats["pool"] += npool
    if npool == 0:
        stats["nopool"] += 1
    arms = [("base", ["-cpsweeps", "1", "-cppasses", "1", "-cpcheck"], opt),
            ("hard", ["-cpsweeps", "3", "-cppasses", "5"], opt),
            ("fathom", ["-cpsweeps", "1"], opt - 1)]
    tol = 1e-9 * max(1.0, abs(opt))
    for arm, extra, cut in arms:
        rc, txt = run([DNNCUTS] + regions + ["-c", repr(cut), "-T", "1", "-v", "-I", "200", "-cpadmm", pool,
                      "-cpiters", "800", "-x", x] + extra + [dat], os.path.join(outdir, "%s.%s.log" % (tag, arm)))
        if rc != 0:
            failures.append("%s %s: rc=%d" % (tag, arm, rc)); continue
        for line in txt.splitlines():
            if line.startswith("cp-selftest"):
                stats["selftests"] += 1
                if not line.rstrip().endswith("ok"): failures.append("%s %s: %s" % (tag, arm, line))
        # the full-precision certified value of every check (the "cp it" lines)
        certs = [float(c) for c in re.findall(r"^cp it +\d+ cert (\S+) ", txt, re.M) if c not in ("-inf", "inf")]
        m = re.search(r"^cp solve  lb (\S+) .*status (\S+)", txt, re.M)
        if not m:
            failures.append("%s %s: no cp solve line" % (tag, arm)); continue
        lb, status = float(m.group(1)), m.group(2)
        stats["certs"] += len(certs); stats["runs"] += 1
        worst = max(certs + [lb])
        stats["worst"] = max(stats["worst"], worst - opt)
        r0 = re.search(r"^round  0  lb (\S+)", txt, re.M)
        lb0 = float(r0.group(1)) if r0 else float("-inf")
        if arm == "hard": stats["gain"].append((lb - lb0, opt - lb0))
        if worst > opt + tol:
            failures.append("%s %s: certified %.9g > optimum %.9g" % (tag, arm, worst, opt))
        if arm != "fathom" and status == "fathomed":
            failures.append("%s %s: fathomed with the cutoff at the optimum" % (tag, arm))
        if arm == "fathom" and status == "fathomed":
            stats["fathomed"] += 1
            if not lb > cut: failures.append("%s fathom: fathomed but lb %.9g <= cutoff" % (tag, lb))


def main():
    global args_core
    p = argparse.ArgumentParser()
    p.add_argument("--core", type=int, default=2)
    p.add_argument("--count", type=int, default=2)
    p.add_argument("--out", default=os.path.join(D, "logs", "validate"))
    a = p.parse_args()
    args_core = a.core
    os.makedirs(a.out, exist_ok=True)
    bf = load_bf()
    rng = random.Random(1107)
    failures = []
    stats = dict(runs=0, certs=0, worst=-np.inf, fathomed=0, selftests=0, pool=0, nopool=0, gain=[])
    nodes = 0
    for family in bf.FAMILIES:
        for n in (7, 8, 9):
            for k in range(a.count):
                A, B, C = bf.make_instance(family, n, rng)
                dat = os.path.join(a.out, "%s_n%d_%d.dat" % (family, n, k))
                bf.write_dat(dat, A, B, C)
                P, v = values(A, B, C)
                specs = [[]]
                for nfix in (1, 2):
                    fac = rng.sample(range(n), nfix); loc = rng.sample(range(n), nfix)
                    specs.append(list(zip(fac, loc)))
                for s, fixes in enumerate(specs):
                    if fixes:
                        mask = np.all(P[:, [i for i, _ in fixes]] == np.array([j for _, j in fixes]), axis=1)
                        opt = float(v[mask].min())
                    else:
                        opt = float(v.min())
                    check_node(dat, fixes, opt, ["-prefix"], "%s_n%d_%d_s%d" % (family, n, k, s), a.out, failures, stats)
                    nodes += 1
    # nug12, matrices swapped: the grid distance first (side 0 regions on it)
    lines = open(os.path.join(ROOT, "qaplib", "nug12.dat")).read().split()
    n = int(lines[0]); nums = list(map(int, lines[1:1 + 2 * n * n]))
    F = [nums[i * n:(i + 1) * n] for i in range(n)]
    Dm = [nums[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]
    dat = os.path.join(a.out, "nug12_swapped.dat")
    bf.write_dat(dat, Dm, F, None)
    check_node(dat, [], 578.0, ["-grid", "3:4"], "nug12s_root", a.out, failures, stats)
    nodes += 1
    for k in range(4):
        fac = rng.sample(range(n), 3); loc = rng.sample(range(n), 3)
        fixes = list(zip(fac, loc))
        opt = node_opt_fixed(Dm, F, None, fixes)
        check_node(dat, fixes, opt, ["-grid", "3:4"], "nug12s_node%d" % k, a.out, failures, stats)
        nodes += 1
    g = stats["gain"]
    closed = [x / y for x, y in g if y > 1e-6]
    print("CP-ADMM validity: %d nodes, %d dnncuts -cpadmm runs, %d certified values, %d self-test lines"
          % (nodes, stats["runs"], stats["certs"], stats["selftests"]))
    print("  pools: %d cuts in all, %d nodes without a violated cut" % (stats["pool"], stats["nopool"]))
    print("  max(certified - optimum) = %.3g (must be <= 1e-9 |optimum|); fathomed at cutoff opt-1: %d"
          % (stats["worst"], stats["fathomed"]))
    if closed:
        print("  hard arm: share of the uncut gap closed: median %.3f, max %.3f over %d nodes with a gap"
              % (float(np.median(closed)), max(closed), len(closed)))
    if failures:
        print("FAILURES (%d):" % len(failures))
        for f in failures[:40]: print("  " + f)
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
