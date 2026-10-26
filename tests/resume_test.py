#!/usr/bin/env python3
"""Check that stopping and resuming the solver changes nothing.

  python3 tests/resume_test.py [--solver PATH] [--slice SECONDS]

Each case is solved once straight through and once in slices: run with
-t SLICE -o STATE, then -r STATE until the solver exits 0. The final node,
fathomed, enumerated and Frank-Wolfe counts, root bound, objective and
per-level node counts must be identical. One case runs without an upper
bound, so the incumbent has to survive the restarts.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CASES = [  # instance, upper bound (None = solver default)
    ("nug15", 1150), ("had14", 2724), ("tai12b", 39464925),
    ("chr12a", 9552), ("nug14", None),
]


def summary(text):
    keys = {
        "nodes": r"Nodes:\s*(\d+\s*/\s*\d+\s*/\s*\d+)",
        "root bound": r"root bound:\s*(\S+)",
        "fw iterations": r"fw iterations:\s*(\d+)",
    }
    out = {k: (re.search(p, text).group(1) if re.search(p, text) else None)
           for k, p in keys.items()}
    found = re.findall(r"New objective value:\s*(\S+)", text)
    out["objective"] = found[-1] if found else None
    m = re.search(r"\nlevel .*?\n(.*?)\n\n", text, re.S)
    out["per-level nodes"] = (" ".join(l.split()[1] for l in m.group(1).splitlines())
                              if m else None)
    return out


def run(solver, args):
    r = subprocess.run([solver] + args, cwd=os.path.dirname(solver),
                       capture_output=True, text=True, timeout=1800)
    return r.returncode, r.stdout + r.stderr


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    p.add_argument("--slice", type=float, default=0.2, help="CPU seconds per slice")
    args = p.parse_args()
    solver = os.path.abspath(args.solver)
    tmp = tempfile.mkdtemp(prefix="qap-resume-")

    failures = 0
    for name, ub in CASES:
        dat = os.path.join(ROOT, "qaplib", name + ".dat")
        base = ["-p", "PARAM/default.par"] + (["-u", str(ub)] if ub is not None else [])
        rc, full = run(solver, base + [dat])
        if rc != 0:
            print(f"FAIL {name}: straight run exited {rc}")
            failures += 1
            continue

        state = os.path.join(tmp, name + ".state")
        rc, out = run(solver, base + ["-o", state, "-t", str(args.slice), dat])
        slices, objectives = 1, re.findall(r"New objective value:\s*(\S+)", out)
        while rc == 3:
            rc, out = run(solver, ["-p", "PARAM/default.par", "-o", state,
                                   "-t", str(args.slice), "-r", state, dat])
            objectives += re.findall(r"New objective value:\s*(\S+)", out)
            slices += 1
        if rc != 0:
            print(f"FAIL {name}: slice {slices} exited {rc}")
            failures += 1
            continue

        a, b = summary(full), summary(out)
        b["objective"] = objectives[-1] if objectives else None
        diffs = [k for k in a if a[k] != b[k]]
        if diffs:
            failures += 1
            print(f"FAIL {name}: differs after {slices} slices in "
                  + ", ".join(f"{k} ({a[k]} vs {b[k]})" for k in diffs))
        else:
            print(f"ok   {name}: {slices} slices, nodes {a['nodes']}, "
                  f"objective {a['objective'] or ub}")
        if slices == 1:
            print(f"     (note: {name} finished within one slice; lower --slice to test resuming)")

    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    print(f"\n{len(CASES) - failures}/{len(CASES)} resume cases identical")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
