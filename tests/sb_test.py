#!/usr/bin/env python3
"""Check DNN strong branching (sb.*, latest/bnb.C).

  python3 tests/sb_test.py [--solver PATH] [--jobs 4]

For each case (an instance, its optimum as the upper bound, a params file with sb.depths):
 * the run ends on the optimum and its statistics have the strong-branching line and table;
 * QAP_SB_TRACE has one line per strong-branching node, and its solves and iterations add up
   to the statistics';
 * a second run over the first's QAP_SB_PROGRESS file solves no child again and has the same tree;
 * bench/parallel.py's totals (state version 10: the strong-branching counts by level) equal
   the serial run's;
 * with sb.gap out of reach the tree is the tree without sb.depths.
And a value the solver cannot use (depths the DNN bound does not cover, a bad sb.gap, sb.gap
alone, a bad sb.power) stops it with exit code 2.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAR = os.path.join(ROOT, "bench", "parallel.py")
CASES = [  # instance, optimum, params file
    ("tai15b", 51765268, "sb_root.par"),
    ("nug18", 1930, "sb.par"),
]
ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")


def solve(solver, par, ub, name, env=()):
    r = subprocess.run([solver, "-p", par, "-u", str(ub), os.path.join(ROOT, "qaplib", name + ".dat")],
                       cwd=os.path.dirname(solver), capture_output=True, text=True, env=dict(ENV, **dict(env)))
    return r.returncode, r.stdout + r.stderr


def nodes(text):
    m = re.search(r"Nodes:\s*(\d+)\s*/\s*(\d+)\s*/\s*(\d+)", text)
    return tuple(int(x) for x in m.groups()) if m else None


def sb_line(text):
    """(nodes, child solves, iterations) of the statistics' strong-branching line, or None"""
    m = re.search(r"strong branch: (\d+) nodes, (\d+) child solves, (\d+) iterations", text)
    return tuple(int(x) for x in m.groups()) if m else None


def sb_table(text):
    """the strong-branching table's rows as (level, nodes, solves, iterations)"""
    m = re.search(r"Strong branching by level:\n.*?\n((?:\s*\d+\s.*\n)+)", text)
    return [tuple(int(x) for x in l.split()[:4]) for l in m.group(1).splitlines()] if m else None


def variant(par, tmp, name, drop=(), add=()):
    """a copy of a params file without the keys in drop and with the lines in add"""
    lines = [l for l in open(par).read().splitlines() if l.split("=")[0].strip() not in drop]
    path = os.path.join(tmp, name)
    with open(path, "w") as f:
        f.write("\n".join(lines + list(add)) + "\n")
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    solver = os.path.abspath(args.solver)
    tmp = tempfile.mkdtemp(prefix="qap-sb-")
    fails = 0

    def check(ok, what):
        nonlocal fails
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {what}")

    for name, opt, parname in CASES:
        par = os.path.join(ROOT, "tests", "par", parname)
        trace, prog = os.path.join(tmp, name + ".sb.txt"), os.path.join(tmp, name + ".progress.txt")
        rc, out = solve(solver, par, opt, name, {"QAP_SB_TRACE": trace, "QAP_SB_PROGRESS": prog})
        n1, s1, t1 = nodes(out), sb_line(out), sb_table(out)
        best = re.findall(r"New objective value:\s*(\S+)", out)
        check(rc == 0 and n1 is not None and s1 is not None and s1[0] > 0 and s1[1] > 0 and (not best or float(best[-1]) == opt),
              f"{name} {parname}: nodes {n1}, strong branching {s1}")
        check(t1 is not None and tuple(sum(r[k] for r in t1) for k in (1, 2, 3)) == s1,
              f"{name}: the table by level adds up to the line: {t1}")
        tl = [l for l in open(trace)] if os.path.exists(trace) else []
        got = (len(tl), sum(int(re.search(r"\| (\d+) solves", l).group(1)) for l in tl),
               sum(int(re.search(r"(\d+) iterations", l).group(1)) for l in tl))
        check(got == s1, f"{name}: QAP_SB_TRACE has a line per node, {got[1]} solves and {got[2]} iterations")
        nprog = sum(1 for l in open(prog)) if os.path.exists(prog) else 0
        check(nprog == s1[1], f"{name}: QAP_SB_PROGRESS has a line per child solve ({nprog})")

        # again over the progress file: nothing solved twice, the same tree
        trace2 = os.path.join(tmp, name + ".sb2.txt")
        rc, out = solve(solver, par, opt, name, {"QAP_SB_TRACE": trace2, "QAP_SB_PROGRESS": prog})
        n2, s2 = nodes(out), sb_line(out)
        reused = sum(int(re.search(r"(\d+) reused", l).group(1)) for l in open(trace2)) if os.path.exists(trace2) else -1
        lines = lambda path: [re.search(r"chose (\S+)", l).group(1) for l in open(path)]
        check(rc == 0 and n2 == n1 and s2 is not None and s2[0] == s1[0] and s2[1] == 0 and reused == s1[1]
              and lines(trace) == lines(trace2),
              f"{name} over its progress file: nodes {n2}, strong branching {s2}, {reused} children reused, the same lines")

        # bench/parallel.py: the summed statistics equal the serial run's
        r = subprocess.run([sys.executable, PAR, os.path.join(ROOT, "qaplib", name + ".dat"), "--solver", solver,
                            "--jobs", str(args.jobs), "--ub", str(opt), "--params", par],
                           capture_output=True, text=True, env=ENV)
        q = r.stdout
        check(r.returncode == 0 and nodes(q) == n1 and sb_line(q) == s1 and sb_table(q) == t1,
              f"{name} bench/parallel.py: nodes {nodes(q)}, strong branching {sb_line(q)}, table {'equal' if sb_table(q) == t1 else sb_table(q)}")

        # a gap no node has: the tree without strong branching
        off = variant(par, tmp, name + ".off.par", drop=("sb.depths", "sb.all", "sb.k", "sb.iters"))
        gated = variant(par, tmp, name + ".gated.par", add=("sb.gap = 0:1e9",))
        (rc0, o0), (rc1, o1) = solve(solver, off, opt, name), solve(solver, gated, opt, name)
        check(rc0 == 0 and rc1 == 0 and nodes(o0) == nodes(o1) and sb_line(o1) is None and sb_line(o0) is None,
              f"{name} with sb.gap out of reach: nodes {nodes(o1)}, without sb.depths {nodes(o0)}")

    # values the solver refuses
    par = os.path.join(ROOT, "tests", "par", "sb.par")
    for what, kw in (("sb.depths above dnn.depths", dict(drop=("sb.depths",), add=("sb.depths = 0:2",))),
                     ("sb.depths below dnn.depths", dict(drop=("sb.depths",), add=("sb.depths = 1:9",))),
                     ("sb.depths without dnn.warm", dict(drop=("dnn.warm",))),
                     ("a bad sb.gap", dict(add=("sb.gap = 2:x",))),
                     ("sb.gap with decreasing depths", dict(add=("sb.gap = 2:0.1,1:0.2",))),
                     ("sb.gap without sb.depths", dict(drop=("sb.depths",), add=("sb.gap = 1:0.1",))),
                     ("a bad sb.power", dict(add=("sb.power = 0:x",))),
                     ("sb.power starting below sb.depths", dict(add=("sb.power = 2:8",)))):
        rc, out = solve(solver, variant(par, tmp, "bad.par", **kw), 1930, "nug18")
        check(rc == 2 and "sb." in out, f"{what}: exit code {rc}")

    print(f"\n{'all strong-branching checks passed' if not fails else f'{fails} strong-branching check(s) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
