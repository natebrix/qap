#!/usr/bin/env python3
"""Validation of --root-facility (README.md, "Root lines"): a run that forces the
rule's own root facility must reproduce the rule's run task for task.

  validate_force.py REF_RUN FORCED_RUN

For every task the forced run did, against the same task id in REF_RUN: the
node's assignments, the DNN trace (first lb to the printed digits, iterations,
status, warm flag, the cut phase's best and status), the children (count, line,
assignments, inherited and predicted bounds), the per-level statistics. The one
field allowed to differ is the depth-1 nodes' own predicted bound (README).
"""
import json
import os
import sys

TRACE_KEYS = ("fixes", "lb", "primal", "status", "iters", "warm", "pool", "cut", "fathomed", "cut_best", "cut_status",
              "cut_iters")


def load(run, tid):
    p = os.path.join(run, "tasks", tid, "done.json")
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    ref, forced = sys.argv[1], sys.argv[2]
    ok = True
    tids = sorted(t for t in os.listdir(os.path.join(forced, "tasks")) if load(forced, t))
    for t in tids:
        a, b = load(ref, t), load(forced, t)
        bad = []
        if a is None:
            print(f"{t}: not in {ref}")
            ok = False
            continue
        if a["spec"] != b["spec"]:
            bad.append(f"spec {a['spec']} / {b['spec']}")
        ta, tb = (a["trace"] or [{}])[0], (b["trace"] or [{}])[0]
        for k in TRACE_KEYS:
            if ta.get(k) != tb.get(k):
                bad.append(f"trace {k}: {ta.get(k)} / {tb.get(k)}")
        if a["n_children"] != b["n_children"] or a["line"] != b["line"]:
            bad.append(f"children {a['n_children']} {a['line']} / {b['n_children']} {b['line']}")
        ka = a.get("children") or a.get("frontier") or []
        kb = b.get("children") or b.get("frontier") or []
        for x, y in zip(ka, kb):
            if x["spec"] != y["spec"] or x["bound"] != y["bound"] or (x["pred"] != y["pred"] and t != "r"):
                bad.append(f"child {x['spec']} bound {x['bound']}/{y['bound']} pred {x['pred']}/{y['pred']}")
                break
        sa, sb = a["stats"], b["stats"]
        for k in ("n_level", "f_level", "l_level", "e_level", "dnn_counts"):
            if sa[k] != sb[k]:
                bad.append(f"stats {k}")
        if sa["counters"][:3] != sb["counters"][:3]:
            bad.append(f"counters {sa['counters']} / {sb['counters']}")
        what = (f"lb {tb.get('lb')} ({tb.get('status')}, {tb.get('iters')} it"
                + (f", phase {tb.get('cut_status')} {tb.get('cut_best')}" if tb.get("cut") else "") + ")"
                if tb else "no DNN")
        print(f"{t} [{b['spec'] or 'root'}]: {what}; {b['n_children']} children ({b['line']})"
              + (f"; forced from {b.get('rule_line')}" if b.get("rule_line") else "")
              + (": SAME as the rule's run" if not bad else ": DIFFERENT: " + "; ".join(bad)))
        ok = ok and not bad
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
