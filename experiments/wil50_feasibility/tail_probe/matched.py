#!/usr/bin/env python3
"""The matched tail probe (README.md, "The matched probe").

tail_probe.py drew its level-5 and level-6 samples per arm, round robin over
the nodes that arm branched, so the arms bounded mostly different nodes (5 of
30 shared at level 5) and their totals could not be compared. This probe
picks every sampled node by a hash of its assignment instead:

  level L0     the node                                    (reused: L3_0)
  level L0+1   all its children                            (reused: L4_j)
  deeper       of each branched sampled node's children, the k with the
               smallest hash (--k: 3 per parent), down to --max-level
               or until nothing is open

The arms branch on the same line almost everywhere (on 31 level-4 parents of
the first run, the sampled children of all four arms were on one line except
where the line itself differed), so the same parents give the same picks in
every arm. A node is bounded in an arm only if its parent branched there.

Estimates are Horvitz-Thompson: a sampled node stands for w nodes, w = its
parent's w times the parent's open children over the number picked. The
standard error is the ultimate-cluster one at the first sampled level
(stratified by level-(L0+1) parent, each sampled node carrying its own
subtree's estimate). Matched nodes (the same assignment bounded in two arms)
are compared directly.

  python3 matched.py ARM [--slots N] [--k 3] [--max-level 9]
  python3 matched.py --compare [ARM ...]

Tasks run in runs/ARM/, the same directory as tail_probe.py, so its level-3
and level-4 tasks are reused (and its level-5 and level-6 ones where they
picked the same node). Resumable. Results: runs/ARM/matched.json,
runs/ARM/matched.txt, runs/matched_compare.txt.
"""
import argparse
import concurrent.futures as cf
import glob
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import tail_probe as tp  # noqa: E402

ARM_ORDER = ("base", "baseship", "full", "fullwide")


def node_keys(path):
    """(key, node) for every open node of a state file, the key being the node's
    assignment as sorted 'facility:location' pairs (the trace's 'fixes')"""
    st = tp.State(path)
    st.zero_stats()
    out = []
    for text, warm in st.nodes:
        f = text.split()
        n = st.n
        p = [int(x) for x in f[5:5 + n]]
        key = " ".join(f"{i}:{j}" for i, j in enumerate(p) if 0 <= j < n)
        out.append((key, (text, warm)))
    return st, out


def trace_key(path):
    for line in open(path):
        m = re.search(r"fixes ([\d: ]+) \|", line)
        if m:
            return " ".join(sorted(m.group(1).split(), key=lambda s: int(s.split(":")[0])))
    return None


def hrank(seed, key):
    return hashlib.sha1(f"{seed}|{key}".encode()).hexdigest()


def cost(trace):
    rows = tp.trace_rows([trace])
    it = sum(r["iters"] + r["cut_iters"] for r in rows)
    cpu = sum(r["cpu"] + r["cut_cpu"] for r in rows)
    first = sum(r["iters"] for r in rows)
    return dict(it=it, cpu=cpu, first=first, phase=it - first, warm=sum(r["warm"] == 1 for r in rows),
                pool=sum(r["pool"] == 1 for r in rows), cut=sum(r["cut"] for r in rows),
                fathomed=all(r["fathomed"] == 1 for r in rows) if rows else True,
                short=[round(r["cutoff"] - (r["cut_best"] if r["cut"] else r["lb"]), 1) for r in rows])


def done(wd, tag):
    return os.path.exists(os.path.join(wd, tag + ".out.state")) and os.path.exists(os.path.join(wd, tag + ".trace"))


def shrink(path):
    """replace a parent's out.state, once every picked child is done, by its
    'nodes N' count alone (the shipped DNN state is tens of MB)"""
    n = tp.open_count(path)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write(f"shrunk\nnodes {n}\n")
    os.replace(tmp, path)


def run(args):
    env = dict(tp.BASE_ENV, **tp.ARMS[args.arm])
    env["QAP_DNN_CUT_GRID"] = args.grid
    wd = os.path.join(args.runs, args.arm)
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "matched_env.txt"), "w") as f:
        f.write(" ".join(f"{k}={v}" for k, v in sorted(env.items())) + "\n" + " ".join(sys.argv) + "\n")
    log = lambda m: print(f"{time.strftime('%H:%M:%S')} {args.arm}: {m}", flush=True)
    t0 = time.time()
    L0 = args.level
    # the first run's deeper tasks, by assignment, where their outputs are still there
    old = {}
    for tr in glob.glob(os.path.join(wd, "L*_*.trace")):
        tag = os.path.basename(tr)[:-6]
        lv = int(tag[1:].split("_")[0])
        if lv >= L0 + 2 and os.path.exists(os.path.join(wd, tag + ".out.state")):
            k = trace_key(tr)
            if k:
                old[k] = tag
    recs = {}                         # tag -> record
    pk_path = os.path.join(wd, "matched_picks.json")    # parent tag -> [[child tag, key]], for resuming
    saved = json.load(open(pk_path)) if os.path.exists(pk_path) else {}

    # level L0 and L0+1: as tail_probe.py (same tags, so its tasks are reused)
    out, tr = tp.run_task(args, env, None, L0, f"L{L0}_0", wd, spec=args.node)
    nk = tp.open_count(out)
    recs[f"L{L0}_0"] = dict(level=L0, w=1.0, parent=None, key=args.node, children=nk, picked=nk, **cost(tr))
    log(f"level {L0}: {nk} children open")
    tags = [f"L{L0 + 1}_{i}" for i in range(nk)]
    missing = [i for i, t in enumerate(tags) if not done(wd, t)]
    tasks = [(None, t) for t in tags]
    if missing:                       # inputs only for the children not yet bounded
        st = tp.State(out)
        st.zero_stats()
        for i in missing:
            p = os.path.join(wd, tags[i] + ".in.state")
            if not os.path.exists(p):
                st.write(p, [st.nodes[i]])
            tasks[i] = (p, tags[i])
        del st
    log(f"level {L0 + 1}: bounding {nk} on {args.slots} slots")
    with cf.ThreadPoolExecutor(max_workers=args.slots) as ex:
        res = dict(zip(tags, ex.map(lambda pt: tp.run_task(args, env, pt[0], L0 + 1, pt[1], wd), tasks)))
    for t in tags:
        o, tr = res[t]
        recs[t] = dict(level=L0 + 1, w=1.0, parent=f"L{L0}_0", key=trace_key(tr), children=tp.open_count(o), **cost(tr))
    frontier = [t for t in tags if recs[t]["children"] > 0]
    log(f"level {L0 + 1}: {len(frontier)} of {nk} branched")

    lv = L0 + 1
    while frontier and lv < args.max_level:
        lv += 1
        picks = []                    # (tag, parent tag, key)
        for pt in frontier:
            if pt not in saved:
                st, nodes = node_keys(os.path.join(wd, pt + ".out.state"))
                order = sorted(range(len(nodes)), key=lambda i: hrank(args.seed, nodes[i][0]))[:args.k]
                ch = []
                for i in order:
                    key = nodes[i][0]
                    t = f"M{lv}_{hashlib.sha1(key.encode()).hexdigest()[:12]}"
                    ch.append([t, key])
                    if done(wd, t):
                        continue
                    o = old.get(key)
                    if o:                 # the first run bounded this very node: reuse it
                        shutil.copyfile(os.path.join(wd, o + ".trace"), os.path.join(wd, t + ".trace"))
                        if not os.path.exists(os.path.join(wd, t + ".out.state")):
                            os.link(os.path.join(wd, o + ".out.state"), os.path.join(wd, t + ".out.state"))
                    elif not os.path.exists(os.path.join(wd, t + ".in.state")):
                        st.write(os.path.join(wd, t + ".in.state"), [nodes[i][1]])
                del st, nodes
                saved[pt] = ch
                json.dump(saved, open(pk_path + ".tmp", "w"))
                os.replace(pk_path + ".tmp", pk_path)
            picks += [(t, pt, key) for t, key in saved[pt]]
        todo = [t for t, *_ in picks if not done(wd, t)]
        log(f"level {lv}: {len(picks)} picked from {len(frontier)} parents, {len(todo)} to bound on {args.slots} slots")
        with cf.ThreadPoolExecutor(max_workers=args.slots) as ex:
            list(ex.map(lambda t: tp.run_task(args, env, os.path.join(wd, t + ".in.state"), lv, t, wd), todo))
        nxt = []
        for t, pt, key in picks:
            tr = os.path.join(wd, t + ".trace")
            o = os.path.join(wd, t + ".out.state")
            P = recs[pt]
            w = P["w"] * P["children"] / len(saved[pt])
            recs[t] = dict(level=lv, w=w, parent=pt, key=key, children=tp.open_count(o), **cost(tr))
            if recs[t]["children"] > 0:
                nxt.append(t)
        for pt in frontier:           # every picked child is done: the parent's shipped state is no longer needed
            if pt != f"L{L0}_0" and not pt.startswith(f"L{L0 + 1}_"):
                shrink(os.path.join(wd, pt + ".out.state"))
        frontier = nxt
        log(f"level {lv}: {len(nxt)} of {len(picks)} branched")
    meta = dict(arm=args.arm, node=args.node, level=L0, k=args.k, seed=args.seed, max_level=args.max_level,
                host=f"{platform.node()} {platform.machine()}", slots=args.slots, wall=time.time() - t0,
                open_below=sum(recs[t]["w"] * recs[t]["children"] for t in frontier), last_level=lv)
    json.dump(dict(meta=meta, nodes=recs), open(os.path.join(wd, "matched.json"), "w"), indent=0)
    txt = summary(meta, recs)
    open(os.path.join(wd, "matched.txt"), "w").write(txt)
    print(txt)


def estimate(meta, recs):
    """per level: est. nodes, iterations, CPU; total; and the ultimate-cluster SE"""
    L0 = meta["level"]
    lv = {}
    for r in recs.values():
        d = lv.setdefault(r["level"], dict(n=0.0, it=0.0, cpu=0.0, bounded=0, fathomed=0, branched=0, warm=0,
                                           pool=0, cut=0, first=0.0, phase=0.0))
        d["n"] += r["w"]; d["it"] += r["w"] * r["it"]; d["cpu"] += r["w"] * r["cpu"]
        d["bounded"] += 1; d["fathomed"] += r["children"] == 0; d["branched"] += r["children"] > 0
        d["warm"] += r["warm"]; d["pool"] += r["pool"]; d["cut"] += r["cut"]
        d["first"] += r["first"]; d["phase"] += r["phase"]
    # subtree totals relative to each node (y = own + descendants weighted by w_d / w_node)
    kids = {}
    for t, r in recs.items():
        if r["parent"]:
            kids.setdefault(r["parent"], []).append(t)

    def sub(t, f):
        r = recs[t]
        return r[f] + sum(recs[c]["w"] / r["w"] * sub(c, f) for c in kids.get(t, []))
    se = {}
    for f in ("it", "cpu"):
        var = 0.0
        for p, cs in kids.items():
            if recs[p]["level"] != L0 + 1:
                continue
            N, k = recs[p]["children"], len(cs)
            if k < 2:
                continue
            ys = [sub(c, f) for c in cs]
            m = sum(ys) / k
            s2 = sum((y - m) ** 2 for y in ys) / (k - 1)
            var += N * N * (1 - k / N) * s2 / k
        se[f] = math.sqrt(var)
    return lv, se


def summary(meta, recs):
    lv, se = estimate(meta, recs)
    lines = [f"# matched tail probe, arm {meta['arm']}: node {meta['node']} (depth {meta['level']}), "
             f"k = {meta['k']} per branched parent below level {meta['level'] + 1}, seed {meta['seed']}",
             f"# {meta['host']}, {meta['slots']} slots, wall {meta['wall']:.0f} s (this session)", "",
             "| level | est. nodes | bounded | fathomed | branched | iterations per bounded node (first + phase) | "
             "est. iterations | est. CPU s | warm | pooled | phases |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in sorted(lv):
        d = lv[k]
        b = d["bounded"]
        lines.append(f"| {k} | {d['n']:,.0f} | {b} | {d['fathomed']} | {d['branched']} | "
                     f"{(d['first'] + d['phase']) / b:,.0f} ({d['first'] / b:,.0f} + {d['phase'] / b:,.0f}) | "
                     f"{d['it']:,.0f} | {d['cpu']:,.0f} | {d['warm']} | {d['pool']} | {d['cut']} |")
    ti = sum(d["it"] for d in lv.values()); tc = sum(d["cpu"] for d in lv.values()); tn = sum(d["n"] for d in lv.values())
    lines += ["", f"Estimated through level {meta['last_level']}: {tn:,.0f} DNN nodes, {ti:,.0f} DNN iterations "
              f"(± {se['it']:,.0f}, 1 SE), {tc:,.0f} CPU s (± {se['cpu']:,.0f}). "
              f"Open below level {meta['last_level']} (estimated): {meta['open_below']:,.0f} nodes."]
    return "\n".join(lines) + "\n"


def compare(runs, arms):
    data = {}
    for a in arms:
        p = os.path.join(runs, a, "matched.json")
        if os.path.exists(p):
            data[a] = json.load(open(p))
    arms = [a for a in arms if a in data]
    out = ["| arm | est. DNN nodes | est. iterations (± 1 SE) | est. CPU s (± 1 SE) | open below | per level: nodes, fathomed/bounded |",
           "|---|---|---|---|---|---|"]
    tot = {}
    for a in arms:
        meta, recs = data[a]["meta"], data[a]["nodes"]
        lv, se = estimate(meta, recs)
        ti = sum(d["it"] for d in lv.values()); tc = sum(d["cpu"] for d in lv.values())
        tn = sum(d["n"] for d in lv.values())
        tot[a] = (ti, tc)
        cells = "; ".join(f"L{k}: {d['n']:,.0f}, {d['fathomed']}/{d['bounded']}" for k, d in sorted(lv.items()))
        out.append(f"| {a} | {tn:,.0f} | {ti:,.0f} (± {se['it']:,.0f}) | {tc:,.0f} (± {se['cpu']:,.0f}) | "
                   f"{meta['open_below']:,.0f} | {cells} |")
    out += ["", "Matched nodes (the same assignment bounded in both arms), by level: "
            "fathomed in each, DNN iterations in each and their ratio", "",
            "| arms | level | shared | fathomed (first / second) | iterations (first / second) | ratio |",
            "|---|---|---|---|---|---|"]
    for i, a in enumerate(arms):
        for b in arms[i + 1:]:
            ka = {(r["level"], r["key"]): r for r in data[a]["nodes"].values()}
            kb = {(r["level"], r["key"]): r for r in data[b]["nodes"].values()}
            for lvl in sorted({k[0] for k in ka} | {k[0] for k in kb}):
                sh = [k for k in ka if k[0] == lvl and k in kb]
                if not sh:
                    continue
                ia = sum(ka[k]["it"] for k in sh); ib = sum(kb[k]["it"] for k in sh)
                fa = sum(ka[k]["children"] == 0 for k in sh); fb = sum(kb[k]["children"] == 0 for k in sh)
                out.append(f"| {a} / {b} | {lvl} | {len(sh)} | {fa} / {fb} | {ia:,} / {ib:,} | "
                           f"{ia / ib if ib else float('nan'):.2f} |")
    out.append("")
    for ref in ("base", "baseship"):
        if ref in tot:
            for a in arms:
                if a != ref:
                    out.append(f"- {a} against {ref}: {tot[ref][0] / tot[a][0]:.2f}x fewer estimated iterations, "
                               f"{tot[ref][1] / tot[a][1]:.2f}x less estimated CPU")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arm", nargs="*")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--node", default="34:40,28:22,14:29")
    ap.add_argument("--level", type=int, default=3)
    ap.add_argument("--instance", default="qaplib/wil50.dat")
    ap.add_argument("--ub", type=float, default=48816)
    ap.add_argument("--params", default="PARAM/nug_deep_xbar9_all.par")
    ap.add_argument("--extra", default="-x 0:4:1:300:1")
    ap.add_argument("--grid", default="5:10")
    ap.add_argument("--slots", type=int, default=None)
    ap.add_argument("--k", type=int, default=3, help="children picked per branched sampled parent")
    ap.add_argument("--max-level", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--runs", default=os.path.join(D, "runs"))
    ap.add_argument("--solver", default=os.path.join(tp.R, "latest", "bnbdrv"))
    args = ap.parse_args()
    if args.compare:
        txt = compare(args.runs, args.arm or list(ARM_ORDER))
        open(os.path.join(args.runs, "matched_compare.txt"), "w").write(txt)
        print(txt)
        return
    if len(args.arm) != 1 or args.arm[0] not in tp.ARMS:
        ap.error(f"one arm of {sorted(tp.ARMS)}")
    args.arm = args.arm[0]
    if args.slots is None:
        args.slots = max(2, (os.cpu_count() or 2) // 2)
    run(args)


if __name__ == "__main__":
    main()
