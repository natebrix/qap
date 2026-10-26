#!/usr/bin/env python3
"""Check bench/parallel.py against the serial solver.

  python3 tests/parallel_test.py [--jobs 4]

 * With the optimum as upper bound, node / fathom / enumerate and FW
   iteration totals must equal a serial run exactly.
 * Without an upper bound, the run must finish on the known optimum.
 * A run with frequent checkpoints is interrupted (SIGTERM) and resumed with
   --resume; its totals must still equal the serial run.
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

    print(f"\n{'all parallel checks passed' if not fails else f'{fails} parallel check(s) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
