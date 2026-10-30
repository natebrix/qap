#!/usr/bin/env python3
"""Validate the certified DNN bound (latest/dnn.C) through latest/dnntest.

  make -C latest dnntest
  python3 tests/dnn_test.py            # all five parts, a few minutes on one core
  python3 tests/dnn_test.py --quick    # fewer instances, no nug20 nodes

A. Brute force. Small random instances (n = 3..8, the four families of
   tests/brute_force.py) and random subproblems of them (0 to n-1 facilities
   fixed). The optimum of each node is found by enumeration. Every bound
   variant must give lb <= optimum, and must not report "fathomed" with the
   cutoff set to the optimum. A second run with cutoff optimum - 1 checks that
   a fathomed status always comes with lb > cutoff.
B. Agreement with the Python prototype (experiments/dnn_bound/nug20/dnn.py).
   dnntest -P computes the prototype's certificate; its value at every check
   must match dnn.py's within 1e-6 relative (nug12 root and subproblems, and,
   unless --quick, nug20 depth-2/3 nodes).
C. Nodes of nug20 with known subproblem optima (solved exactly by bnbdrv in
   the nug20 study): the default bound must stay below them, and fathom
   exactly the nodes the study fathomed. A second run must reproduce every
   result bit for bit.
E. Eliminated assignments (dnntest -E): random small instances and nodes with
   random forbidden assignments; the optimum over the allowed permutations is
   found by enumeration. The bound must stay at or below it and not fathom
   with the cutoff at it; a node with no allowed permutation must fathom.
D. Cancellation regressions: data whose products are large and cancel
   (against the linear term, or between asymmetric A and B), where every
   rounding step before the certificate matters. Each node's exact optimum
   is computed with Python integers; the bound must not exceed it and must
   not fathom with the cutoff at the optimum.

Needs numpy. Runs on one core (OPENBLAS_NUM_THREADS=1 is set for children).
"""

import argparse
import importlib.util
import itertools
import os
import re
import subprocess
import sys
import tempfile
import time

# one core: also for numpy in this process (the prototype runs here)
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DNNTEST = os.path.join(ROOT, "latest", "dnntest")
QAPLIB = os.path.join(ROOT, "qaplib")
PROTO = os.path.join(ROOT, "experiments", "dnn_bound", "nug20")
ENV = dict(os.environ)

# bound variants checked in part A: (name, dnntest flags)
MODES = [("default", []), ("ssyevd", ["-e", "1"]), ("dsyevd", ["-e", "0"]),
         ("lap-reduce", ["-l", "0"]), ("lap-exact", ["-l", "1"]), ("box", ["-b"]),
         ("prototype", ["-P", "-s", "0"])]


def load_brute_force():
    spec = importlib.util.spec_from_file_location(
        "brute_force", os.path.join(ROOT, "tests", "brute_force.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_dnntest(dat, nodes, flags, iters):
    """nodes: list of (fixes string, cutoff). Returns list of result dicts."""
    with tempfile.NamedTemporaryFile("w", suffix=".nodes", delete=False) as f:
        for fx, c in nodes:
            f.write("%s %s\n" % (fx, c if isinstance(c, str) else repr(c)))
        nf = f.name
    try:
        out = subprocess.run([DNNTEST, "-N", nf, "-i", str(iters)] + flags + [dat],
                             capture_output=True, text=True, env=ENV, timeout=600)
    finally:
        os.unlink(nf)
    res = []
    for line in out.stdout.splitlines():
        if line.startswith("RESULT"):
            d = dict(kv.split("=", 1) for kv in line.split()[1:])
            res.append(d)
    if out.returncode != 0 or len(res) != len(nodes):
        raise RuntimeError("dnntest failed (rc=%d):\n%s\n%s" % (out.returncode, out.stdout,
                                                                 out.stderr))
    return res


def all_values(A, B, C):
    """Objective of every permutation (n <= 8), as (perms, values)."""
    n = len(A)
    P = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    v = np.einsum("ik,pik->p", A, B[P[:, :, None], P[:, None, :]])
    if C is not None:
        v = v + np.asarray(C, float)[np.arange(n), P].sum(axis=1)
    return P, v


def part_a(args, rng):
    bf = load_brute_force()
    tmp = tempfile.mkdtemp(prefix="dnn-test-")
    count = 1 if args.quick else args.count
    worst = worst_dnn = -np.inf
    nnodes = nruns = ntight = nfathom = 0
    failures = []
    t0 = time.time()
    for family in bf.FAMILIES:
        for n in range(3, 9):
            for k in range(count):
                A, B, C = bf.make_instance(family, n, rng)
                dat = os.path.join(tmp, "%s_n%d_%d.dat" % (family, n, k))
                bf.write_dat(dat, A, B, C)
                P, v = all_values(A, B, C)
                # nodes: root, then random partial assignments with 1..n-1 fixed
                nodes = [("-", float(v.min()))]
                for nfix in sorted(set([1, n // 2, n - 2, n - 1]) - {0}):
                    fac = rng.sample(range(n), nfix)
                    loc = rng.sample(range(n), nfix)
                    mask = np.all(P[:, fac] == np.array(loc), axis=1)
                    fx = ",".join("%d:%d" % (i, l) for i, l in zip(fac, loc))
                    nodes.append((fx, float(v[mask].min())))
                for name, flags in MODES:
                    if args.quick and name in ("lap-exact", "box", "dsyevd"):
                        continue
                    # cutoff = optimum: fathoming here would be wrong
                    res = run_dnntest(dat, nodes, flags, args.iters)
                    # cutoff = optimum - 1: fathoming allowed, must have lb > cutoff
                    res2 = run_dnntest(dat, [(fx, opt - 1) for fx, opt in nodes], flags,
                                       args.iters)
                    nruns += 2
                    for (fx, opt), r, r2 in zip(nodes, res, res2):
                        nnodes += 1
                        lb = float(r["lb"])
                        tolv = 1e-9 * max(1.0, abs(opt))
                        worst = max(worst, lb - opt)
                        if int(r["m"]) >= 2:   # m <= 1 is solved exactly, lb = optimum
                            worst_dnn = max(worst_dnn, lb - opt)
                        if lb >= opt - 1e-3 * max(1.0, abs(opt)):
                            ntight += 1
                        where = "%s %s node [%s] opt=%g" % (name, os.path.basename(dat), fx, opt)
                        if lb > opt + tolv:
                            failures.append("%s: lb %.9g > optimum" % (where, lb))
                        if r["status"] == "fathomed":
                            failures.append("%s: fathomed with cutoff = optimum (lb %.9g)" % (where, lb))
                        if r["status"] == "error" or r2["status"] == "error":
                            failures.append("%s: error status" % where)
                        lb2 = float(r2["lb"])
                        if lb2 > opt + tolv:
                            failures.append("%s: lb %.9g > optimum (cutoff opt-1)" % (where, lb2))
                        if r2["status"] == "fathomed":
                            nfathom += 1
                            if not lb2 > opt - 1:
                                failures.append("%s: fathomed but lb %.9g <= cutoff" % (where, lb2))
    print("A. brute force: %d node bounds in %d dnntest runs, %.0f s" % (nnodes, nruns,
                                                                        time.time() - t0))
    print("   max(lb - optimum) = %.3g, over nodes with m >= 2: %.3g (must be <= 0 up to"
          " 1e-9 relative)" % (worst, worst_dnn))
    print("   bound within 1e-3 of the optimum: %d of %d; fathomed at cutoff opt-1: %d"
          % (ntight, nnodes, nfathom))
    return failures


def proto_nodes(quick):
    """(instance, fixes, iterations) for part B."""
    nodes = [("nug12", "-", 900), ("nug12", "0:3,5:7", 300), ("nug12", "2:2,7:0,9:11", 300)]
    if not quick:
        nodes += [("nug20", "0:19,1:0", 300), ("nug20", "1:7,19:11", 200),
                  ("nug20", "1:1,12:19,16:15", 200)]
    return nodes


def part_b(args):
    sys.path.insert(0, PROTO)
    import qaputil
    import dnn as proto
    failures = []
    t0 = time.time()
    worst = 0.0
    nvals = 0
    for inst, fx, iters in proto_nodes(args.quick):
        dat = os.path.join(QAPLIB, inst + ".dat")
        n, a, b, c = qaputil.read_qaplib(dat)
        p = [-1] * n
        if fx != "-":
            for t in fx.split(","):
                i, l = map(int, t.split(":"))
                p[i] = l
        F, L, AFF, BLL, Cp, const = qaputil.reduce_subproblem(a, b, c, p)
        py = proto.dnn_admm(AFF, BLL, Cp, const, maxit=iters, check_every=50)
        pyhist = {h[0]: h[1] for h in py["hist"]}
        out = subprocess.run([DNNTEST, "-P", "-s", "0", "-v", "-i", str(iters), "-x", fx, dat],
                             capture_output=True, text=True, env=ENV)
        chist = {}
        for line in out.stdout.splitlines():
            mt = re.match(r"it\s+(\d+) lb \S+ \(cert (\S+),", line)
            if mt:
                chist[int(mt.group(1))] = float(mt.group(2))
        if set(chist) != set(pyhist):
            failures.append("%s [%s]: checks differ: C %s python %s" % (inst, fx, sorted(chist),
                                                                       sorted(pyhist)))
            continue
        for it in sorted(pyhist):
            d = abs(chist[it] - pyhist[it]) / max(1.0, abs(pyhist[it]))
            worst = max(worst, d)
            nvals += 1
            if d > 1e-6:
                failures.append("%s [%s] it %d: C %.6f python %.6f" % (inst, fx, it, chist[it],
                                                                        pyhist[it]))
        print("   %s [%s] m=%d: certified lb after %d iterations C %.4f, python %.4f"
              % (inst, fx, len(F), iters, chist[max(chist)], pyhist[max(pyhist)]))
    print("B. prototype agreement: %d checkpoint values, max relative difference %.2g, %.0f s"
          % (nvals, worst, time.time() - t0))
    return failures


# nug20 nodes from the nug20 study (verify.log): subproblem optimum from bnbdrv,
# and whether the study's certified bound fathomed them at incumbent 2570
NUG20_EXACT = [("1:6,19:12", 2570, False), ("1:6,19:13", 2588, False),
               ("1:7,19:11", 2594, True), ("1:6,19:8", 2594, True), ("0:19,1:0", 2646, True),
               ("1:1,12:19,16:15", 2584, True), ("1:6,15:4,19:10", 2644, True)]


def part_c(args):
    failures = []
    t0 = time.time()
    dat = os.path.join(QAPLIB, "nug20.dat")
    nodes = [(fx, 2570 - 1.95) for fx, _, _ in NUG20_EXACT]
    res = run_dnntest(dat, nodes, [], 1200)
    again = run_dnntest(dat, nodes, [], 1200)
    for r, r2 in zip(res, again):
        if any(r[k] != r2[k] for k in r if k not in ("cpu", "ms_per_iter")):
            failures.append("not reproducible: %s vs %s" % (r, r2))
    tight = run_dnntest(dat, [(fx, opt) for fx, opt, _ in NUG20_EXACT], ["-i", "400"], 400)
    for (fx, opt, fath), r, r2 in zip(NUG20_EXACT, res, tight):
        lb = float(r["lb"])
        ok = lb <= opt and float(r2["lb"]) <= opt and r2["status"] != "fathomed"
        print("   nug20 [%s] optimum %d: lb %.3f %s in %s iterations (%.1f s), study %s"
              % (fx, opt, lb, r["status"], r["iters"], float(r["cpu"]),
                 "fathomed" if fath else "not fathomed"))
        if not ok:
            failures.append("nug20 [%s]: lb %s / %s vs optimum %d" % (fx, r["lb"], r2["lb"], opt))
        if (r["status"] == "fathomed") != fath:
            failures.append("nug20 [%s]: fathomed=%s, study %s" % (fx, r["status"], fath))
    print("C. nug20 nodes with known optima: %.0f s" % (time.time() - t0))
    return failures


def part_d(args):
    """Cancellation regressions (exact optima by integer enumeration)."""
    t0 = time.time()
    failures = []
    cases = []
    # Astra's review: a b + c = 3 with a b = 2^64 + 2^34 + 3, at m = 1 and
    # after reducing a 3 x 3 instance with facility 0 fixed (shift = 3 exactly)
    a, b, c = 2**32 + 1, 2**32 + 3, -(2**64 + 2**34)
    cases.append(("m1", [[a]], [[b]], [[c]], "-"))
    Z = lambda v: [[v, 0, 0], [0, 0, 0], [0, 0, 0]]
    cases.append(("reduced", Z(a), Z(b), Z(c), "0:0"))
    # products near 2^62 cancelling in the shift and in the linear terms
    r = __import__("random").Random(5)
    for k in range(4):
        n = 5
        big = 2**31 + r.randrange(2**20)
        A = [[r.randrange(-3, 4) for _ in range(n)] for _ in range(n)]
        B = [[r.randrange(-3, 4) for _ in range(n)] for _ in range(n)]
        for M in (A, B):
            for i in range(n):
                for j in range(i):
                    M[i][j] = M[j][i]
        A[0][1] = A[1][0] = big
        B[0][1] = B[1][0] = big + 1
        C = [[r.randrange(-9, 10) for _ in range(n)] for _ in range(n)]
        C[1][1] = -2 * big * (big + 1) + r.randrange(-5, 6)
        cases.append(("big%d" % k, A, B, C, "0:0"))
    # neither matrix symmetric: large antisymmetric parts cancel in L
    for k in range(3):
        n = 5
        S = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                v = r.randrange(1, 2**26)
                S[i][j], S[j][i] = v, -v
        P = [[r.randrange(0, 9) for _ in range(n)] for _ in range(n)]
        A = [[P[i][j] + P[j][i] + S[i][j] for j in range(n)] for i in range(n)]
        B = [[r.randrange(0, 2**26) for _ in range(n)] for _ in range(n)]
        C = [[0] * n for _ in range(n)]
        cases.append(("asym%d" % k, A, B, C, "-"))
    # found by adversarial testing (tests/data/dnn_cancel): asymmetric A and B
    # cancelling in L, and products cancelling in a subproblem's shift and C
    from fractions import Fraction
    for name, fx in [("ac_nonint_3", "-"), ("ac_int_4", "-"),
                     ("rc_nonint_1", "0:0,1:1,2:2"), ("rc_int_3", "0:0,1:1,2:2")]:
        t = open(os.path.join(ROOT, "tests", "data", "dnn_cancel", name + ".dat")).read().split()
        n = int(t[0])
        v = [Fraction(x) for x in t[1:1 + 3 * n * n]]
        v += [Fraction(0)] * (3 * n * n - len(v))
        M = [[v[k * n * n + i * n:k * n * n + (i + 1) * n] for i in range(n)] for k in range(3)]
        cases.append((name, M[0], M[1], M[2], fx))
    with tempfile.TemporaryDirectory() as d:
        for name, A, B, C, fx in cases:
            n = len(A)
            fixed = {}
            if fx != "-":
                for t in fx.split(","):
                    i, j = map(int, t.split(":"))
                    fixed[i] = j
            best = None
            for p in itertools.permutations(range(n)):
                if any(p[i] != j for i, j in fixed.items()):
                    continue
                v = sum(A[i][k] * B[p[i]][p[k]] for i in range(n) for k in range(n)) + \
                    sum(C[i][p[i]] for i in range(n))
                best = v if best is None else min(best, v)
            dat = os.path.join(ROOT, "tests", "data", "dnn_cancel", name + ".dat")
            if not os.path.exists(dat):
                dat = os.path.join(d, name + ".dat")
                with open(dat, "w") as f:
                    f.write("%d\n\n" % n)
                    for M in (A, B, C):
                        f.write("\n".join(" ".join(str(x) for x in row) for row in M) + "\n\n")
            # cutoff: the optimum, rounded down to a double, so that a
            # fathom at it is wrong
            cut = float(best)
            if Fraction(cut) > best:
                cut = np.nextafter(cut, -np.inf)
            res = run_dnntest(dat, [(fx, repr(float(cut)))], [], 300)[0]
            lb = float(res["lb"])
            # compare exactly: the double lb against the exact optimum
            if Fraction(lb) > best or res["fathomed"] == "1":
                failures.append("D %s: lb %r > optimum %s or fathomed (%s)"
                                % (name, lb, float(best), res["fathomed"]))
    print("D. cancellation regressions: %d cases, %.0f s" % (len(cases), time.time() - t0))
    return failures


def part_e(args, rng):
    """Eliminated assignments."""
    t0 = time.time()
    failures = []
    stronger = same = cases = 0
    with tempfile.TemporaryDirectory() as d:
        for inst in range(16 if args.quick else 40):
            n = rng.choice([5, 6, 7, 8])
            A = [[0] * n for _ in range(n)]
            B = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    A[i][j] = A[j][i] = rng.randrange(0, 20)
                    B[i][j] = B[j][i] = rng.randrange(1, 20)
            C = [[rng.randrange(0, 30) if rng.random() < 0.3 else 0 for _ in range(n)]
                 for _ in range(n)]
            dat = os.path.join(d, "e%d.dat" % inst)
            with open(dat, "w") as f:
                f.write("%d\n\n" % n)
                for M in (A, B, C):
                    f.write("\n".join(" ".join(map(str, r)) for r in M) + "\n\n")
            P, vals = all_values(A, B, C)
            nodes = []
            for _ in range(4):
                nfix = rng.randrange(0, n - 2)
                perm = list(range(n)); rng.shuffle(perm)
                rows = rng.sample(range(n), nfix)
                fixed = {i: perm[i] for i in rows}
                free_r = [i for i in range(n) if i not in fixed]
                free_c = [j for j in range(n) if j not in fixed.values()]
                prob = rng.choice([0.2, 0.4, 0.6])
                elim = [(i, j) for i in free_r for j in free_c if rng.random() < prob]
                ok = np.ones(len(P), bool)
                for i, j in fixed.items():
                    ok &= P[:, i] == j
                for i, j in elim:
                    ok &= P[:, i] != j
                fx = ",".join("%d:%d" % (i, j) for i, j in sorted(fixed.items())) or "-"
                el = ",".join("%d:%d" % e for e in elim) or "-"
                if ok.any():
                    opt = float(vals[ok].min())
                    nodes.append((fx, opt, el, False))
                else:
                    nodes.append((fx, 1e9, el, True))
            res = run_dnntest_e(dat, nodes, 600)
            plain = run_dnntest(dat, [(fx, repr(opt)) for fx, opt, el, inf in nodes], [], 600)
            for (fx, opt, el, infeasible), r, rp in zip(nodes, res, plain):
                cases += 1
                lb = float(r["lb"])
                if infeasible:
                    if r["fathomed"] != "1":
                        failures.append("E %s fixes %s elim %s: infeasible node not fathomed" % (dat, fx, el))
                    continue
                if lb > opt + 1e-9 or r["fathomed"] == "1":
                    failures.append("E %s fixes %s elim %s: lb %r > optimum %r or fathomed"
                                    % (dat, fx, el, lb, opt))
                if lb > float(rp["lb"]) + 1e-6:
                    stronger += 1
                else:
                    same += 1
    print("E. eliminated assignments: %d nodes (%d with a stronger bound than without "
          "eliminations), %.0f s" % (cases, stronger, time.time() - t0))
    return failures


def run_dnntest_e(dat, nodes, iters):
    """nodes: (fixes, cutoff, eliminations, infeasible); returns RESULT dicts."""
    with tempfile.NamedTemporaryFile("w", suffix=".nodes", delete=False) as f:
        for fx, c, el, _ in nodes:
            f.write("%s %r %s\n" % (fx, c, el))
        nf = f.name
    try:
        out = subprocess.run([DNNTEST, "-N", nf, "-i", str(iters)] + [dat],
                             capture_output=True, text=True, env=ENV, timeout=600)
    finally:
        os.unlink(nf)
    res = [dict(kv.split("=", 1) for kv in l.split()[1:])
           for l in out.stdout.splitlines() if l.startswith("RESULT")]
    if out.returncode != 0 or len(res) != len(nodes):
        raise RuntimeError("dnntest failed (rc=%d):\n%s\n%s" % (out.returncode, out.stdout,
                                                                 out.stderr))
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--count", type=int, default=2, help="instances per family and size (A)")
    p.add_argument("--iters", type=int, default=1500, help="iteration limit in part A")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--parts", default="ABCDE")
    args = p.parse_args()
    if not os.path.exists(DNNTEST):
        sys.exit("%s not found; run make -C latest dnntest" % DNNTEST)
    import random
    rng = random.Random(args.seed)
    failures = []
    if "A" in args.parts:
        failures += part_a(args, rng)
    if "B" in args.parts:
        failures += part_b(args)
    if "C" in args.parts and not args.quick:
        failures += part_c(args)
    if "D" in args.parts:
        failures += part_d(args)
    if "E" in args.parts:
        failures += part_e(args, rng)
    if failures:
        print("FAILED: %d problems" % len(failures))
        for f in failures[:30]:
            print("  " + f)
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
