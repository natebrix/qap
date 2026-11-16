"""[e] What early entry and the near hand-over would save on the laptop's complete wil50 subtrees (no
histories there), from what they did on a220's matched nodes here:
 - early: a warm give-up's first solve ends at 110 iterations instead of 300+ (its phase no longer: a220);
 - near hand-over: a missed pooled phase ends by about CAP iterations (a220: 880->810, 1620->1020,
   3400->1100, 990->760; fathoming phases unchanged), its children's cost unchanged (1,228->1,142, 1,144->1,144).
  python3 extrap.py [CAP=1050]"""
import sys
from tracelib import load
CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 1050
print("| subtree | DNN nodes | iterations | first solves that fathom | give-ups' first solves | fathoming phases | missed phases | [e] early saves | [e] hand-over saves | [e] factor |")
print("|---|---|---|---|---|---|---|---|---|---|")
for name in ("a220", "a251", "a325"):
    n = load("../../estimates/wil50/subtrees/%s/trace.txt" % name)
    tot = sum(e["iters"] + e["cit"] for e in n)
    ff = sum(e["iters"] for e in n if e["status"] == "fathomed")
    gu = sum(e["iters"] for e in n if e["status"] != "fathomed")
    pf = sum(e["cit"] for e in n if e["cut"] and e["cstatus"] == "fathomed")
    pm = sum(e["cit"] for e in n if e["cut"] and e["cstatus"] != "fathomed")
    early = sum(e["iters"] - 110 for e in n if e["status"] != "fathomed" and e["warm"])
    ho = sum(max(0, e["cit"] - CAP) for e in n if e["cut"] and e["cstatus"] != "fathomed" and e["pool"])
    p = lambda x: "%s (%.0f%%)" % (format(x, ","), 100.0 * x / tot)
    print("| %s (laptop) | %d | %s | %s | %s | %s | %s | %s | %s | %.2fx |" % (name, len(n), format(tot, ","), p(ff), p(gu), p(pf), p(pm), p(early), p(ho), tot / (tot - early - ho)))
