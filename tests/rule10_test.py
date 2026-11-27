#!/usr/bin/env python3
"""Rules 10 and 11 (the lookahead rules 4 and 9 with the size score, branch.la_exp).

  python3 tests/rule10_test.py [--jobs 4]

 * with every line a candidate the root of nug30 branches on row 23 and the root of tho30 on row 5
   (experiments/branching_lab/ucheap: the rows the score names), and rule 4 on the same candidates takes
   row 0 at nug30, so the nug30 check tells the two scores apart (tho30's does not: both take row 5);
 * a params file that uses rule 10 without branch.la_exp is refused (exit code 2);
 * nug18 under the DNN bound with rule 10 at every level closes at the optimum, serially and through
   bench/parallel.py with the same totals;
 * with "best _ rows" 1 the rule has one row and one column to choose from and still closes;
 * rule 11 (rule 9's marginal lines, the same score) does the same.
"""
import argparse, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAR = os.path.join(ROOT, "bench", "parallel.py")
ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
BLOCK = """ {rule}   branching rule
 150  Frank-Wolfe iterations, for bound
 100  Frank-wolfe iterations, when we know we can't fathom
 25   Frank-Wolfe iterations, for row selection
 1    which bound?
 30   step for parametric bound
 {best}   compute B for best _ rows of U
"""


def par(tmp, name, best, settings, rule=10):
    path = os.path.join(tmp, name)
    with open(path, "w") as f:
        f.write("0    has the location of a symmetry file been provided?\n-1    restriction on tree depth\n"
                "1    number of branching rules\n\n" + BLOCK.format(rule=rule, best=best) + "\n[settings]\n"
                + "\n".join(settings) + "\n")
    return path


def run(solver, args, env=()):
    r = subprocess.run([solver] + args, cwd=os.path.dirname(solver), capture_output=True, text=True, env=dict(ENV, **dict(env)))
    return r.returncode, r.stdout + r.stderr


def nodes(text):
    m = re.search(r"Nodes:\s*(\d+)\s*/\s*(\d+)\s*/\s*(\d+)", text)
    return tuple(int(x) for x in m.groups()) if m else None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    solver = os.path.abspath(args.solver)
    tmp = tempfile.mkdtemp(prefix="qap-rule10-")
    fails = 0

    def check(ok, what):
        nonlocal fails
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {what}")

    allp = par(tmp, "all.par", 99, ["branch.la_exp = 32"])
    rule4 = par(tmp, "rule4.par", 99, ["ship.bits = 64"], rule=4)
    for name, ub, row, pf, what in (("nug30", 6126, "r23", allp, "rule 10"), ("tho30", 149938, "r5", allp, "rule 10"),
                                    ("nug30", 6126, "r0", rule4, "rule 4 (the same candidates and child bounds, its own score)")):
        trace = os.path.join(tmp, name + what[:7].replace(" ", "") + ".root.txt")
        rc, out = run(solver, ["-p", pf, "-s", f"SYM/{name}.sym", "-u", str(ub), "-d", "0",
                               os.path.join(ROOT, "qaplib", name + ".dat")], {"QAP_NODE_TRACE": trace})
        f = open(trace).read().split() if os.path.exists(trace) else []
        check(rc == 0 and len(f) > 9 and f[9] == row, f"{name}'s root, every row, {what}: branches on {f[9] if len(f) > 9 else None} (expected {row})")

    for rule in (10, 11):
        rc, out = run(solver, ["-p", par(tmp, "noexp.par", 99, ["ship.bits = 64"], rule=rule), "-u", "1930", os.path.join(ROOT, "qaplib", "nug18.dat")])
        check(rc == 2 and "branch.la_exp" in out, f"rule {rule} without branch.la_exp: exit code {rc}")

    dnn = ["ship.bits = 64", "dnn.depths = 1:4", "dnn.giveup = 20", "dnn.warm = 1", "branch.la_exp = 32"]
    for label, best, rule in (("rule 10, every line", 99, 10), ("rule 10, one row and one column", 1, 10),
                              ("rule 11 (rule 9's lines)", 5, 11)):
        pf = par(tmp, f"dnn{rule}_{best}.par", best, dnn, rule=rule)
        rc, out = run(solver, ["-p", pf, "-u", "1930", os.path.join(ROOT, "qaplib", "nug18.dat")])
        n1 = nodes(out)
        best_obj = re.findall(r"New objective value:\s*(\S+)", out)
        check(rc == 0 and n1 is not None and (not best_obj or float(best_obj[-1]) == 1930), f"nug18, {label}: nodes {n1}")
        r = subprocess.run([sys.executable, PAR, os.path.join(ROOT, "qaplib", "nug18.dat"), "--solver", solver, "--jobs", str(args.jobs),
                            "--ub", "1930", "--params", pf], capture_output=True, text=True, env=ENV)
        check(r.returncode == 0 and nodes(r.stdout) == n1, f"nug18, {label}, bench/parallel.py: nodes {nodes(r.stdout)}")

    print(f"\n{'all rule 10 checks passed' if not fails else f'{fails} rule 10 check(s) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
