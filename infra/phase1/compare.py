#!/usr/bin/env python3
"""Phase 1 (README.md): the acceptance check of a finished run.

  compare.py RESULTS [REFERENCE] [--fixed-ub] [--opt X]

RESULTS is the coordinator's results.txt, REFERENCE the log of
bench/parallel.py run with the same settings on the same machine. With a
fixed incumbent (--fixed-ub) the two searches are the same tree, so nodes,
fathomed, enumerated, FW iterations and the DNN counts must be equal, and the
best objective must be the optimum. Without one, incumbents arrive in an
order that depends on timing, so only the optimum is checked and the totals
are shown side by side. --opt is the known optimum (QAPSOLN). Prints a table
and exits 0 when the check passes, 1 when it fails.
"""
import argparse
import re
import sys


def totals(text):
    def num(pat, cast=float):
        m = re.search(pat, text)
        return cast(m.group(1)) if m else None
    m = re.search(r"Nodes: (\d+) / (\d+) / (\d+)", text)
    d = re.search(r"DNN bound: (\d+) nodes, (\d+) fathomed", text)
    nodes = tuple(map(int, m.groups())) if m else (None, None, None)
    dnn = tuple(map(int, d.groups())) if d else (0, 0)
    best = re.findall(r"best objective (\S+)", text)
    return {"nodes": nodes[0], "fathomed": nodes[1], "enumerated": nodes[2],
            "fw iterations": num(r"fw iterations: (\d+)", int),
            "DNN nodes": dnn[0], "DNN fathomed": dnn[1],
            "best objective": float(best[-1]) if best else None,
            "CPU s": num(r"Total Time: (\S+) seconds"),
            "wall s": num(r"wall (\S+) s,")}


EXACT = ("nodes", "fathomed", "enumerated", "fw iterations", "DNN nodes", "DNN fathomed")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("results")
    p.add_argument("reference", nargs="?")
    p.add_argument("--fixed-ub", action="store_true")
    p.add_argument("--opt", type=float)
    a = p.parse_args(argv)
    got = totals(open(a.results).read())
    ref = totals(open(a.reference).read()) if a.reference else None
    bad = []
    if got["nodes"] is None:
        bad.append("no totals in the coordinator's results")
    if a.opt is not None and got["best objective"] != a.opt:
        bad.append(f"best objective {got['best objective']}, the optimum is {a.opt:g}")
    if ref and a.fixed_ub:
        if ref["nodes"] is None:
            bad.append("no totals in the reference")
        bad += [f"{k}: {got[k]} (coordinator) != {ref[k]} (parallel.py)" for k in EXACT if got[k] != ref[k]]
    w = max(map(len, got))
    print(f"{'':{w}}  {'coordinator':>14}" + (f"  {'parallel.py':>14}" if ref else ""))
    for k in got:
        def f(v):
            return "-" if v is None else (f"{v:.1f}" if isinstance(v, float) and k.endswith(" s") else f"{v:g}"
                                          if isinstance(v, float) else str(v))
        mark = ""
        if ref and a.fixed_ub and k in EXACT:
            mark = "  =" if got[k] == ref[k] else "  DIFFERENT"
        print(f"{k:{w}}  {f(got[k]):>14}" + (f"  {f(ref[k]):>14}{mark}" if ref else ""))
    what = "totals equal to parallel.py and the optimum" if ref and a.fixed_ub else "the optimum"
    if a.opt is None and not (ref and a.fixed_ub):
        what = "totals present (no optimum known, nothing to compare)"
    print(("PASS: " + what) if not bad else "FAIL:\n  " + "\n  ".join(bad))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
