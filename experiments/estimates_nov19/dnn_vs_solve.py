#!/usr/bin/env python3
"""The same nodes bounded by an estimate's dives and by a solve, from their DNN traces (QAP_DNN_TRACE lines,
latest/bnb.C: "depth D fixes a:p ... | first lb X ... iters N cpu C warm W pool P kept K dropped D | cut 0|1
[best .. rounds .. iters N cpu C ...] | cutoff ..."): matched by their fixes, the ADMM iterations and CPU of
each, the first solve and the cut phase apart, and how each started (warm, with an inherited pool).

    python3 experiments/estimates_nov19/dnn_vs_solve.py infra/phase1/runs/sko42-est2-261121 infra/phase1/runs/sko42-261120b

The estimate's trace is RUN/dives/dnn.txt.gz (a dives run with QAP_EST_DNN_TRACE=1), the solve's RUN/dnn.txt.gz.
A node the estimator bounded more than once (in several dives) counts once, its first line.
"""
import collections
import gzip
import os
import re
import statistics as st
import sys

NUM = r"(-?[0-9.eE+-]+|nan|-?inf)"


def parse(line):
    head, _, rest = line.partition("| first")
    fx = frozenset(head.split("fixes", 1)[1].split())
    g = lambda pat, s: (lambda m: float(m.group(1)) if m else None)(re.search(pat, s))
    first, _, after = rest.partition("| cut")
    cut, _, tail = after.partition("| cutoff")
    return fx, {
        "depth": int(re.search(r"depth (\d+)", head).group(1)),
        "it1": g(r"iters (\d+)", first), "cpu1": g(r"cpu " + NUM, first),
        "warm": g(r"warm (\d+)", first), "pool": g(r"pool (\d+)", first),
        "kept": g(r"kept (\d+)", first), "dropped": g(r"dropped (\d+)", first),
        "lb1": g(r"lb " + NUM, first), "status1": (re.search(r"status (\S+)", first) or [None, None])[1],
        "cut": cut.strip().startswith("1"), "itc": g(r"iters (\d+)", cut) or 0.0, "cpuc": g(r"cpu " + NUM, cut) or 0.0,
        "retry": "retry" in first,
    }


def load(path):
    out = {}
    for l in gzip.open(path, "rt"):
        if "fixes" in l and "| first" in l:
            k, v = parse(l)
            out.setdefault(k, v)
    return out


def main():
    est_dir, solve_dir = sys.argv[1:3]
    E = load(os.path.join(est_dir, "dives", "dnn.txt.gz"))
    S = load(os.path.join(solve_dir, "dnn.txt.gz"))
    both = [k for k in E if k in S]
    print(f"estimate: {len(E):,} nodes in its DNN trace; solve: {len(S):,}; in both: {len(both):,}")

    def row(name, ks):
        if not ks:
            return
        e_it = sum(E[k]["it1"] + E[k]["itc"] for k in ks); s_it = sum(S[k]["it1"] + S[k]["itc"] for k in ks)
        e_c = sum(E[k]["cpu1"] + E[k]["cpuc"] for k in ks); s_c = sum(S[k]["cpu1"] + S[k]["cpuc"] for k in ks)
        e1 = sum(E[k]["it1"] for k in ks); s1 = sum(S[k]["it1"] for k in ks)
        print(f"| {name} | {len(ks):,} | {e_it / s_it:.2f} | {e1 / s1:.2f} | {e_c / s_c:.2f} | "
              f"{(e_c / e_it) / (s_c / s_it):.2f} | {sum(1 for k in ks if E[k]['cut']):,} / {sum(1 for k in ks if S[k]['cut']):,} |")
    print("\n| nodes | n | iterations est/solve | first solve's iterations | CPU est/solve | CPU per iteration | cut phases est / solve |")
    print("|---|---|---|---|---|---|---|")
    row("all in both", both)
    for d in sorted({E[k]["depth"] for k in both}):
        row(f"depth {d}", [k for k in both if E[k]["depth"] == d])
    for name, f in (("estimate warm, solve warm", lambda k: E[k]["warm"] and S[k]["warm"]),
                    ("estimate with a pool, solve with a pool", lambda k: E[k]["pool"] and S[k]["pool"]),
                    ("estimate without a pool, solve with one", lambda k: not E[k]["pool"] and S[k]["pool"]),
                    ("estimate cold, solve warm", lambda k: not E[k]["warm"] and S[k]["warm"]),
                    ("the same outcome of the first solve", lambda k: E[k]["status1"] == S[k]["status1"])):
        row(name, [k for k in both if f(k)])
    kd = [(E[k]["kept"], S[k]["kept"]) for k in both if E[k]["pool"] and S[k]["pool"]]
    if kd:
        print(f"\ninherited cuts kept, median: estimate {st.median(a for a, _ in kd):.0f}, solve {st.median(b for _, b in kd):.0f}")
    lb = [E[k]["lb1"] - S[k]["lb1"] for k in both if E[k]["lb1"] is not None and S[k]["lb1"] is not None]
    if lb:
        print(f"first solve's bound, estimate - solve: median {st.median(lb):+.3f}")


if __name__ == "__main__":
    main()
