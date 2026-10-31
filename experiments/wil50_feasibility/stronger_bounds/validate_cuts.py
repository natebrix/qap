#!/usr/bin/env python3
"""Validity check of dnncuts' regional cuts (-prefix) and atomic cuts against
enumeration on random small instances: every certified value printed by every
round must be at or below the optimum of the node, and the final bound must not
exceed it. Also records how often the cuts strengthen the bound.

  python3 experiments/wil50_feasibility/stronger_bounds/validate_cuts.py [-n 7] [-N 12] [--seed 1]
"""
import argparse, importlib.util, os, random, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DNNCUTS = os.path.join(ROOT, "experiments", "dnn_cuts", "dnncuts")
spec = importlib.util.spec_from_file_location("brute_force", os.path.join(ROOT, "tests", "brute_force.py"))
bf = importlib.util.module_from_spec(spec); spec.loader.exec_module(bf)


def node_optimum(A, B, C, fixes):
    """optimum over permutations with the fixed assignments (A-index i -> B-index j)"""
    n = len(A)
    fixed = dict(fixes)
    best = None
    def rec(k, used, p, cost):
        nonlocal best
        if k == n:
            if best is None or cost < best: best = cost
            return
        locs = [fixed[k]] if k in fixed else [l for l in range(n) if l not in used and l not in fixed.values()]
        for loc in locs:
            d = A[k][k] * B[loc][loc] + (C[k][loc] if C else 0)
            for i in range(k):
                d += A[i][k] * B[p[i]][loc] + A[k][i] * B[loc][p[i]]
            p.append(loc); used.add(loc)
            rec(k + 1, used, p, cost + d)
            p.pop(); used.discard(loc)
    rec(0, set(), [], 0)
    return best


def run(dat, fixes, flags):
    cmd = [DNNCUTS, "-x", fixes, "-e", "0", "-g", "1.618", "-G", "0", "-r", "6", "-I", "1500",
           "-i", "800", "-k", "30", "-delta", "2"] + flags + [dat]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("dnncuts failed: %s\n%s" % (" ".join(cmd), out.stderr))
    lbs = [float(m.group(1)) for m in re.finditer(r"round\s+\d+\s+lb (\S+)", out.stdout)]
    fin = re.search(r"FINAL best_lb (\S+) first_lb (\S+)", out.stdout)
    return lbs, float(fin.group(1)), float(fin.group(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=7)
    ap.add_argument("-N", type=int, default=12, help="instances")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    tmp = tempfile.mkdtemp()
    checks = bad = 0
    gains = {"prefix": [], "atomic": []}
    for t in range(args.N):
        fam = bf.FAMILIES[t % len(bf.FAMILIES)]
        A, B, C = bf.make_instance(fam, args.n, rng)
        dat = os.path.join(tmp, "i%d.dat" % t)
        bf.write_dat(dat, A, B, C)
        # root and a node with one or two fixed assignments
        nodes = ["-"]
        i, j = rng.randrange(args.n), rng.randrange(args.n)
        nodes.append("%d:%d" % (i, j))
        i2 = (i + 1 + rng.randrange(args.n - 1)) % args.n; j2 = (j + 1 + rng.randrange(args.n - 1)) % args.n
        nodes.append("%d:%d,%d:%d" % (i, j, i2, j2))
        for fx in nodes:
            fixes = [] if fx == "-" else [tuple(map(int, s.split(":"))) for s in fx.split(",")]
            opt = node_optimum(A, B, C, fixes)
            for name, flags in (("prefix", ["-prefix"]), ("atomic", [])):
                lbs, best, first = run(dat, fx, flags)
                for lb in lbs:
                    checks += 1
                    if lb > opt + 1e-6:
                        bad += 1
                        print("VIOLATION: %s node %s %s: lb %.6f > optimum %d" % (dat, fx, name, lb, opt))
                gains[name].append((best - first, opt - first))
                print("%-7s n=%d node %-10s %-7s opt %5d  plain %9.4f  with cuts %9.4f  (closes %5.1f%% of %6.3f)" % (
                    fam, args.n, fx, name, opt, first, best,
                    100 * (best - first) / (opt - first) if opt - first > 1e-9 else 100.0, opt - first))
    print("\n%d certified values checked against enumeration, %d above the optimum" % (checks, bad))
    for name, g in gains.items():
        gapped = [(a, b) for a, b in g if b > 1e-6]
        print("%s cuts: %d nodes, %d with a DNN gap; mean fraction closed %.1f%%" % (
            name, len(g), len(gapped), 100 * sum(a / b for a, b in gapped) / max(1, len(gapped))))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
