#!/usr/bin/env python3
"""Check bench/parallel.py against the serial solver.

  python3 tests/parallel_test.py [--jobs 4]

 * With the optimum as upper bound, node / fathom / enumerate and FW
   iteration totals must equal a serial run exactly.
 * Without an upper bound, the run must finish on the known optimum.
 * A run with frequent checkpoints is interrupted (SIGTERM) and resumed with
   --resume; its totals must still equal the serial run.
 * With the DNN bound (state version 7: per-level DNN counts), the tree by
   level, the DNN table's nodes and fathom and give-up rates and the fathoms
   by cause, must equal a serial run's, and parallel.py must print those
   tables (the solver's own, from the summed statistics).
"""

import argparse
import os
import re
import signal
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "latest", "bnbdrv")
PAR = os.path.join(ROOT, "bench", "parallel.py")


def totals(text):
    m = re.search(r"Nodes:\s*(\d+\s*/\s*\d+\s*/\s*\d+)", text)
    f = re.search(r"fw iterations:\s*(\d+)", text)
    o = re.findall(r"New objective value:\s*(\S+)", text)
    return (m.group(1) if m else None, f.group(1) if f else None,
            float(o[-1]) if o else None)


def serial(name, ub):
    out = subprocess.run([SOLVER, "-p", "PARAM/default.par", "-u", str(ub),
                          os.path.join(ROOT, "qaplib", name + ".dat")],
                         cwd=os.path.dirname(SOLVER), capture_output=True, text=True).stdout
    return totals(out)


def parallel(name, jobs, *extra):
    out = subprocess.run([sys.executable, PAR, os.path.join(ROOT, "qaplib", name + ".dat"),
                          "--jobs", str(jobs)] + list(extra),
                         capture_output=True, text=True).stdout
    return totals(out)


def dnn_tree(text):
    """the totals, and the DNN and fathoms-by-cause tables without their
    iteration and CPU columns (the tree, not the work spent on it)"""
    d = re.search(r"DNN bound:\s*(\d+) nodes, (\d+) fathomed", text)
    tables = []
    for title, keep in (("DNN bound by level", 4), ("Fathoms by cause", 6)):
        m = re.search(title + r":\n[^\n]*\n((?: ?\d+ [^\n]*\n)*)", text)
        tables.append([l.split()[:keep] for l in m.group(1).splitlines()] if m else None)
    return totals(text)[:2], d.groups() if d else None, tables


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    fails = 0

    for name, ub in (("nug15", 1150), ("tai12b", 39464925), ("had14", 2724)):
        s, q = serial(name, ub), parallel(name, args.jobs, "--ub", str(ub))
        ok = s[:2] == q[:2]
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {name} --ub opt: serial {s[0]}, fw {s[1]}; "
              f"parallel {q[0]}, fw {q[1]}")

    for name, opt in (("nug14", 1014), ("tai12b", 39464925)):
        q = parallel(name, args.jobs)
        ok = q[2] == opt
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {name} no upper bound: found {q[2]}, optimum {opt}")

    # interrupt a checkpointed run, then resume it
    name, ub = "nug16a", 1610
    s = serial(name, ub)
    tmp = tempfile.mkdtemp(prefix="qap-ptest-")
    ck = os.path.join(tmp, "ck.state")
    proc = subprocess.Popen([sys.executable, PAR, os.path.join(ROOT, "qaplib", name + ".dat"),
                             "--jobs", str(args.jobs), "--ub", str(ub), "--checkpoint", ck,
                             "--checkpoint-every", "0.5"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deadline = time.time() + 30
    while not os.path.exists(ck) and time.time() < deadline:
        time.sleep(0.1)
    time.sleep(0.6)
    proc.send_signal(signal.SIGTERM)
    proc.wait()
    if not os.path.exists(ck):
        print(f"FAIL {name}: no checkpoint written")
        fails += 1
    else:
        q = parallel(name, args.jobs, "--ub", str(ub), "--resume", ck)
        ok = s[:2] == q[:2]
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {name} interrupted and resumed: serial {s[0]}, "
              f"fw {s[1]}; resumed {q[0]}, fw {q[1]}")
        os.remove(ck)
    os.rmdir(tmp)
    # the DNN bound, its per-level counts through state version 7
    name, ub, x = "nug16a", 1610, "2:4:1:0:1"
    env = dict(os.environ, QAP_DNN_SHIP="64", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    s_out = subprocess.run([SOLVER, "-p", "PARAM/default.par", "-u", str(ub), "-x", x,
                            os.path.join(ROOT, "qaplib", name + ".dat")],
                           cwd=os.path.dirname(SOLVER), capture_output=True, text=True, env=env).stdout
    q_out = subprocess.run([sys.executable, PAR, os.path.join(ROOT, "qaplib", name + ".dat"),
                            "--jobs", str(args.jobs), "--ub", str(ub), "--extra", f"-x {x}"],
                           capture_output=True, text=True, env=env).stdout
    s, q = dnn_tree(s_out), dnn_tree(q_out)
    ok = s == q and None not in s[2] and s[1] is not None
    fails += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {name} -x {x} --ub opt: serial {s[0][0]}, DNN {s[1]}; "
          f"parallel {q[0][0]}, DNN {q[1]}; per-level DNN and cause tables "
          + ("equal" if s[2] == q[2] else f"differ:\n  serial   {s[2]}\n  parallel {q[2]}"))

    print(f"\n{'all parallel checks passed' if not fails else f'{fails} parallel check(s) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
