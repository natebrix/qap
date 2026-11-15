"""What the DNN give-up (dnn.giveup) costs and saves, from an export's two
traces (qapctl export: nodes.txt.gz and dnn.txt.gz; docs/tho40b_results_plan.md,
(a) 3).

A node whose DNN solve gives up (status hopeless: no fathom in sight after
dnn.giveup iterations) branches. If every one of its children then fathoms at
once, the node itself was close: more iterations might have fathomed it, for
fewer nodes. Per depth: the DNN nodes and their statuses, the give-ups, the
give-ups whose children all fathomed at once ("one level"), and what those
children cost: the most a later give-up could save there, against the
iterations it would add.

usage: python3 giveups.py EXPORT_DIR      (nodes.txt.gz and dnn.txt.gz)
"""
import argparse
import collections
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ntrace import load  # noqa: E402


def dnn_fields(line):
    """the fixes key and the fields of a QAP_DNN_TRACE line"""
    head, rest = line.split("|", 1)
    h = head.split()
    k = h.index("fixes")
    key = tuple(sorted(tuple(int(x) for x in p.split(":")) for p in h[k + 1:]))
    f = rest.split("|", 1)[0].split()   # first lb X primal Y status S iters I cpu C warm W ...
    at = lambda name: f[f.index(name) + 1]
    return key, {"status": at("status"), "iters": int(at("iters")), "cpu": float(at("cpu")), "warm": at("warm") == "1"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export")
    a = ap.parse_args()
    t = load(os.path.join(a.export, "nodes.txt.gz"))
    dnn = {}
    with gzip.open(os.path.join(a.export, "dnn.txt.gz"), "rt") as f:
        for ln in f:
            if ln.strip():
                k, v = dnn_fields(ln)
                dnn[k] = v
    key = lambda i: tuple(sorted(tuple(int(x) for x in p.split(":")) for p in t.fixes[i]))
    rows = collections.defaultdict(lambda: collections.Counter())
    for i in range(len(t)):
        d = dnn.get(key(i))
        if d is None:
            continue
        r = rows[t.depth[i]]
        r["dnn"] += 1
        r["iters"] += d["iters"]
        r[d["status"]] += 1
        r["warm"] += d["warm"]
        if d["status"] == "hopeless":
            r["giveup"] += 1
            r["giveup_iters"] += d["iters"]
            r["giveup_cpu"] += d["cpu"]
            ch = t.children[i]
            if t.outcome[i] == "B" and ch and len(ch) == t.kids[i] and all(t.fathomed(c) for c in ch):
                r["one_level"] += 1
                r["one_level_child_cpu"] += sum(t.cpu[c] for c in ch)
                r["one_level_children"] += len(ch)
            r["giveup_subtree"] += t.size[i]
    print(f"{a.export}: {len(t)} nodes, {len(dnn)} DNN lines")
    print(f"{'depth':>5} {'DNN':>8} {'warm':>6} {'it/node':>8} {'gave up':>8} {'one level':>10} "
          f"{'their kids':>10} {'kids CPU s':>11} {'give-up CPU s':>14} {'subtree/give-up':>16}")
    for dep in sorted(rows):
        r = rows[dep]
        print(f"{dep:>5} {r['dnn']:>8,} {r['warm'] / r['dnn']:>6.0%} {r['iters'] / r['dnn']:>8.1f} {r['giveup']:>8,} "
              f"{r['one_level']:>10,} {r['one_level_children']:>10,} {r['one_level_child_cpu']:>11,.0f} "
              f"{r['giveup_cpu']:>14,.0f} {r['giveup_subtree'] / max(1, r['giveup']):>16,.1f}")
    st = collections.Counter()
    for r in rows.values():
        for s in ("fathomed", "hopeless", "converged", "stalled", "iter_limit", "time_limit", "error"):
            st[s] += r[s]
    print("statuses:", ", ".join(f"{s} {n:,}" for s, n in st.items() if n))


if __name__ == "__main__":
    main()
