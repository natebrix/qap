#!/usr/bin/env python3
"""Compare each replay (runs/<side>/<task>) with stage A's record of the task
(tasks/<task>: trace.txt, solver.log, out.full.state / out.slim.state / out.state)."""
import json
import os
import re
import sys

W = "/home/user/qap/.claude/worktrees/agent-aebd90129de6196e4"
SA = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA"
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
sys.path.insert(0, os.path.join(W, "experiments", "wil50_feasibility", "sampling"))
from common import State, Warm, node_fields, parse_trace, strip_dnn  # noqa: E402
from stage_a import child_line  # noqa: E402

FIELDS = ["depth", "fixes", "lb", "primal", "status", "iters", "warm", "pool", "kept", "cut", "cut_best",
          "cut_rounds", "cut_iters", "cut_status", "cut_cp", "cutoff", "fathomed", "pool_out"]


def norm_trace(path):
    out = []
    for l in open(path):
        l = re.sub(r"^pid \d+ ", "", l.rstrip("\n"))
        l = re.sub(r" cpu [0-9.]+", " cpu *", l)
        out.append(l)
    return out


def log_counts(path):
    txt = open(path, errors="replace").read()
    g = {}
    for key, rx in (("nodes", r"Nodes: (.*)"), ("root_bound", r"root bound: (.*)"), ("fw_iters", r"fw iterations: (.*)"),
                    ("dnn", r"DNN bound: (\d+ nodes, \d+ fathomed)")):
        m = re.search(rx, txt)
        g[key] = m.group(1).strip() if m else None
    m = re.search(r"Strategy breakdown:\nnodes(.*)", txt)
    g["strategy_nodes"] = m.group(1).split() if m else None
    return g


def children(path, n_fix_own):
    st = State(path)
    kids = []
    for text, w in st.nodes:
        d, b, pb, f = node_fields(text, st.n)
        kids.append({"text": text, "depth": d, "bound": b, "pred": pb, "fixes": f, "warm": w.text if w else None})
    return st, kids


def old_children_path(tid):
    for name in ("out.full.state", "out.state", "out.slim.state"):
        p = os.path.join(SA, "tasks", tid, name)
        if os.path.exists(p):
            return p, name
    return None, None


def compare(side, tid):
    rd = os.path.join(G, "runs", side, tid)
    rec = json.load(open(os.path.join(SA, "tasks", tid, "done.json")))
    res = {"side": side, "task": tid, "spec": rec["spec"], "depth": rec["depth"], "diffs": []}
    done = json.load(open(os.path.join(rd, "done.json")))
    res["rc"] = done["rc"]
    res["in_state_same_as_a667727_writer"] = done["in_state_same_as_a667727_writer"]
    # 1. the trace line, pid and cpu aside
    to, tn = norm_trace(os.path.join(SA, "tasks", tid, "trace.txt")), norm_trace(os.path.join(rd, "trace.txt"))
    res["trace_identical"] = to == tn
    po, pn = parse_trace(os.path.join(SA, "tasks", tid, "trace.txt")), parse_trace(os.path.join(rd, "trace.txt"))
    if not res["trace_identical"]:
        if len(po) != len(pn):
            res["diffs"].append(f"trace: {len(po)} lines in the record, {len(pn)} in the replay")
        for a, b in zip(po, pn):
            for k in FIELDS:
                if a.get(k) != b.get(k):
                    res["diffs"].append(f"trace field {k}: record {a.get(k)!r}, replay {b.get(k)!r}")
                    break
        res["trace_record"], res["trace_replay"] = to, tn
    res["rec_trace"], res["rep_trace"] = (po[0] if po else None), (pn[0] if pn else None)
    # 2. the solver log's counters
    lo, ln = log_counts(os.path.join(SA, "tasks", tid, "solver.log")), log_counts(os.path.join(rd, "solver.log"))
    res["log_record"], res["log_replay"] = lo, ln
    for k in lo:
        if lo[k] != ln[k]:
            res["diffs"].append(f"log {k}: record {lo[k]!r}, replay {ln[k]!r}")
    # 3. the children
    own = rec["fixes"]
    own = [tuple(x) for x in own]
    sn, kn = children(os.path.join(rd, "out.state"), len(own))
    res["cpu_replay"] = sn.times[0]
    res["cpu_record"] = rec["stats"]["times"][0]
    op, oname = old_children_path(tid)
    res["old_children_file"] = oname
    if rec["n_children"] == 0 and op is None:
        ko = []
    else:
        so, ko = children(op, len(own))
        if so.root_bound != sn.root_bound or so.obj != sn.obj or so.perm != sn.perm:
            res["diffs"].append(f"header: record root {so.root_bound} obj {so.obj}, replay root {sn.root_bound} obj {sn.obj}")
    res["n_children_record"], res["n_children_replay"] = len(ko), len(kn)
    if len(ko) != len(kn):
        res["diffs"].append(f"children: {len(ko)} in the record, {len(kn)} in the replay")
    res["line_record"] = rec["line"]
    res["line_replay"] = child_line(own, [k["fixes"] for k in kn])
    if res["line_record"] != res["line_replay"]:
        res["diffs"].append(f"line: record {res['line_record']}, replay {res['line_replay']}")
    text_same = all(a["text"] == b["text"] for a, b in zip(ko, kn)) and len(ko) == len(kn)
    val_same = all(a["bound"] == b["bound"] and a["pred"] == b["pred"] and a["fixes"] == b["fixes"]
                   and a["depth"] == b["depth"] for a, b in zip(ko, kn)) and len(ko) == len(kn)
    res["children_text_identical"], res["children_values_identical"] = text_same, val_same
    if not val_same:
        for i, (a, b) in enumerate(zip(ko, kn)):
            if a["bound"] != b["bound"] or a["pred"] != b["pred"] or a["fixes"] != b["fixes"]:
                res["diffs"].append(f"child {i}: record {a['fixes']} {a['bound']!r} {a['pred']!r}, "
                                    f"replay {b['fixes']} {b['bound']!r} {b['pred']!r}")
                break
    if ko:
        res["children_bound"] = sorted({k["bound"] for k in kn})
        res["children_pred_min_max"] = [min(k["pred"] for k in kn), max(k["pred"] for k in kn)]
        res["children_pred_max_absdiff"] = max(abs(a["pred"] - b["pred"]) for a, b in zip(ko, kn)) if val_same or len(ko) == len(kn) else None
        res["children_bound_max_absdiff"] = max(abs(a["bound"] - b["bound"]) for a, b in zip(ko, kn)) if len(ko) == len(kn) else None
        # the warm starts the children carry (the QPB point, the cut pool and, with the full file, the shipped DNN state)
        wo = {k["warm"] for k in ko}
        if oname == "out.slim.state":
            wn = {strip_dnn(Warm(k["warm"])).text if k["warm"] else None for k in kn}
        else:
            wn = {k["warm"] for k in kn}
        res["warm_lines_compared"] = "without DNN state" if oname == "out.slim.state" else "with the shipped DNN state"
        res["warm_identical"] = wo == wn
        res["warm_bytes"] = sum(len(x) for x in wn if x)
        if not res["warm_identical"]:
            a, b = sorted(x or "" for x in wo), sorted(x or "" for x in wn)
            for x, y in zip(a, b):
                tx, ty = x.split(" "), y.split(" ")
                for i, (p, q) in enumerate(zip(tx, ty)):
                    if p != q:
                        res["diffs"].append(f"warm line field {i} (of {len(tx)}/{len(ty)}): record {p[:60]!r}, "
                                            f"replay {q[:60]!r}")
                        break
                break
    res["identical"] = not res["diffs"]
    return res


def main():
    out = []
    for side in ("new", "old", "pre", "fix"):
        base = os.path.join(G, "runs", side)
        if not os.path.isdir(base):
            continue
        for tid in sorted(os.listdir(base)):
            if not os.path.exists(os.path.join(base, tid, "done.json")):
                continue
            r = compare(side, tid)
            out.append(r)
            t = r["rep_trace"]
            tr = (f"{t['status']} lb {t['lb']:.4f} it {t['iters']}" + (f" cut {t.get('cut_status')} best "
                  f"{t.get('cut_best'):.4f} it {t.get('cut_iters')}" if t and t['cut'] else "")) if t else "no DNN"
            print(f"{side} {tid:8s} [{r['spec']}] {'IDENTICAL' if r['identical'] else 'DIFFERS'}: trace "
                  f"{'same' if r['trace_identical'] else 'DIFF'}; {tr}; children {r['n_children_record']}/"
                  f"{r['n_children_replay']} values {'same' if r['children_values_identical'] else 'DIFF'} text "
                  f"{'same' if r['children_text_identical'] else 'DIFF'}; line {r['line_record']}/{r['line_replay']}; "
                  f"warm {r.get('warm_identical')} ({r.get('warm_lines_compared')}, {r.get('warm_bytes', 0)} bytes); "
                  f"cpu record {r['cpu_record']:.1f} replay {r['cpu_replay']:.1f}; in.state same {r['in_state_same_as_a667727_writer']}")
            for d in r["diffs"]:
                print("    ", d)
    with open(os.path.join(G, "compare.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
