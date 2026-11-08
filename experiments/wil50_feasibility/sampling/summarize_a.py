#!/usr/bin/env python3
"""Stage A's numbers (README.md) from a stage_a.py run directory: per depth,
the nodes bounded and fathomed, the DNN iterations and CPU, the cut phases,
the children kept and eliminated, the branching lines, the frontier and the
symmetry checks. Writes OUT/stageA_summary.md and OUT/frontier_stats.tsv.

  summarize_a.py [RUN] [--out results]
"""
import argparse
import collections
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import D  # noqa: E402


def load(run):
    docs = {}
    for t in os.listdir(os.path.join(run, "tasks")):
        p = os.path.join(run, "tasks", t, "done.json")
        if os.path.exists(p):
            docs[t] = json.load(open(p))
    return docs


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", nargs="?", default=os.path.join(D, "runs", "stageA"))
    ap.add_argument("--out", default=os.path.join(D, "results"))
    a = ap.parse_args()
    cfg = json.load(open(os.path.join(a.run, "config.json")))
    docs = load(a.run)
    Dmax = cfg["depth"]
    finished = os.path.exists(os.path.join(a.run, "frontier.tsv"))
    cutoff = None
    L = []
    w = L.append
    w(f"# Stage A: depths 0-{Dmax} of {os.path.basename(cfg['instance'])} "
      f"({'finished' if finished else 'IN PROGRESS'}; {len(docs)} tasks)")
    w("")
    c = cfg["commit"]
    w(f"- binary: commit `{c['commit'][:12]}`, latest/ modified: {c['latest_modified']}, sources newer than the "
      f"binary: {c['sources_newer_than_binary'] or 'none'} (built {c['built']})")
    w(f"- `bnbdrv -p {os.path.basename(cfg['params'])} -x {cfg['xspec']} -s {os.path.relpath(cfg['sym'], os.path.join(D, '../../..'))} "
      f"-u {cfg['ub']:g} -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh"
      + (f" plus {cfg['env_overrides']}" if cfg["env_overrides"] else " unchanged"))
    w(f"- stage-B sample seed {cfg['seed']}, parents kept below hash {cfg['tau']}")
    sym_ok = sum(1 for d in docs.values() if d["checks"]["sym_loaded"] and not d["checks"]["sym_missing"])
    w(f"- symmetry file loaded in {sym_ok} of {len(docs)} solver logs (\"read symmetry info from ...\")")
    w("")
    by = collections.defaultdict(list)
    for d in docs.values():
        by[d["depth"]].append(d)
    w("| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | "
      "CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | "
      "children kept | children eliminated in branching |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    tot = collections.Counter()
    phases = []
    for depth in sorted(by):
        ds = by[depth]
        nb = sum(d["stats"]["n_level"][depth] for d in ds)
        nf = sum(d["stats"]["f_level"][depth] for d in ds)
        tr = [r for d in ds for r in d["trace"]]
        if tr:
            cutoff = tr[0]["cutoff"]
        dn = len(tr)
        df = sum(r["fathomed"] for r in tr)
        it1 = sum(r["iters"] for r in tr)
        itc = sum(r.get("cut_iters", 0) for r in tr)
        it = it1 + itc
        dcpu = sum(d["stats"]["dnn_time"][0] for d in ds)       # the first solves (bnb.C: dnn_time)
        ccpu = sum(r.get("cut_cpu", 0.0) for r in tr)            # the cut phases (trace)
        cpu = sum(d["stats"]["times"][0] for d in ds)
        warm = sum(r["warm"] for r in tr)
        ph = [r for r in tr if r["cut"]]
        phases += ph
        kept = sum(d["n_children"] for d in ds)
        elim = sum(d["stats"]["e_level"][depth] for d in ds)
        w(f"| {depth} | {nb} | {nf} | {dn} | {df} | {it1:,} + {itc:,} = {it:,} | "
          f"{dcpu:,.0f} + {ccpu:,.0f} + {cpu - dcpu - ccpu:,.0f} = {cpu:,.0f} | {cpu / max(1, nb):,.1f} | "
          f"{warm} | {len(ph)} | {kept:,} | {elim:,} |")
        tot.update(nodes=nb, fathomed=nf, dnn=dn, dnn_f=df, it1=it1, itc=itc, dnn_cpu=dcpu, cut_cpu=ccpu, cpu=cpu)
    w(f"| all | {tot['nodes']} | {tot['fathomed']} | {tot['dnn']} | {tot['dnn_f']} | "
      f"{tot['it1']:,} + {tot['itc']:,} = {tot['it1'] + tot['itc']:,} | {tot['dnn_cpu']:,.0f} + {tot['cut_cpu']:,.0f} + "
      f"{tot['cpu'] - tot['dnn_cpu'] - tot['cut_cpu']:,.0f} = {tot['cpu']:,.0f} | | | {len(phases)} | | |")
    w("")
    wall = sum(d["wall"] for d in docs.values())
    w(f"Task wall time summed over slots: {wall:,.0f} s; CPU/wall {tot['cpu'] / wall if wall else 0:.3f}.")
    w("")
    # the DNN results per depth
    w(f"Cutoff {cutoff}. First DNN solve at each depth (status counts; certified bound at the end of the solve):")
    w("")
    w("| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |")
    w("|---|---|---|---|---|---|")
    for depth in sorted(by):
        tr = [r for d in by[depth] for r in d["trace"]]
        for st in sorted({r["status"] for r in tr}):
            rs = [r for r in tr if r["status"] == st]
            its = [r["iters"] for r in rs]
            sh = [r["cutoff"] - r["lb"] for r in rs]
            w(f"| {depth} | {st} | {len(rs)} | {min(its)} / {statistics.median(its):.0f} / {max(its)} | "
              f"{min(sh):.1f} / {statistics.median(sh):.1f} / {max(sh):.1f} | "
              f"{statistics.median(r['cpu'] for r in rs):.1f} |")
    w("")
    # shortfall histogram of the last depth's first solves (after the phase, if one ran)
    lastrows = [r for d in by.get(Dmax, []) for r in d["trace"]]
    if lastrows:
        edges = [-1e9, 0, 50, 100, 150, 200, 244, 300, 350, 400, 450, 500, 600, 1e9]
        names = ["fathomed", "0-50", "50-100", "100-150", "150-200", "200-244", "244-300", "300-350", "350-400",
                 "400-450", "450-500", "500-600", ">= 600"]
        cnt = [0] * len(names)
        for r in lastrows:
            s = -1.0 if r["fathomed"] else r["short"]
            for k in range(len(names)):
                if edges[k] < s <= edges[k + 1] or (k == 0 and s <= 0):
                    cnt[k] += 1
                    break
        w(f"Depth-{Dmax} nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; "
          f"the cut zone is 244):")
        w("")
        w("| " + " | ".join(names) + " |")
        w("|" + "---|" * len(names))
        w("| " + " | ".join(str(x) for x in cnt) + " |")
        w("")
    if phases:
        w("Cut phases:")
        w("")
        w("| node | depth | first lb | phase best | status | iterations | CPU s |")
        w("|---|---|---|---|---|---|---|")
        for r in phases:
            w(f"| {r['fixes']} | {r['depth']} | {r['lb']:.2f} | {r.get('cut_best', float('nan')):.2f} | "
              f"{r.get('cut_status')} | {r.get('cut_iters')} | {r.get('cut_cpu', 0):.0f} |")
    else:
        w("Cut phases: none ran (no first solve ended within the zone, 0.005 x 48,816 = 244 units of the cutoff).")
    w("")
    # branching lines
    w("Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each "
      "child's location, `col j` = location j):")
    w("")
    for depth in range(Dmax + 1):
        ds = sorted(by.get(depth, []), key=lambda d: [int(x) if x.isdigit() else -1 for x in d["id"].split(".")])
        cnt = collections.Counter(d["line"].split()[0] for d in ds if d["n_children"])
        w(f"- depth {depth}: {len(ds)} nodes, {sum(1 for d in ds if d['n_children'])} branched: "
          + ", ".join(f"{v} on a {k}" for k, v in sorted(cnt.items())))
        if depth < Dmax and ds:
            w("  " + "; ".join(f"`{d['spec'] or 'root'}` {d['line']} ({d['n_children']})" for d in ds))
    w("")
    if Dmax >= 2 and by.get(1):
        w(f"Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, "
          f"their lines, the frontier nodes they left):")
        w("")
        w("| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |")
        w("|---|---|---|---|---|---|---|---|")
        for d in sorted(by[1], key=lambda d: int(d["id"].split(".")[1])):
            kids = [x for x in by.get(2, []) if x["parent"] == d["id"]]
            t = d["trace"][0] if d["trace"] else None
            lc = collections.Counter(x["line"] for x in kids if x["n_children"])
            w(f"| `{d['spec']}` | {t['short']:.1f} | {d['line']} | {d['n_children']} | {len(kids)} | "
              f"{sum(1 for x in kids if x['n_children'] == 0)} | "
              + ", ".join(f"{k} x{v}" for k, v in lc.most_common(4)) + (" ..." if len(lc) > 4 else "")
              + f" | {sum(x['n_children'] for x in kids):,} |")
        w("")
    # children counts at the last depth
    last = by.get(Dmax, [])
    kc = collections.Counter(d["n_children"] for d in last)
    w(f"Children per depth-{Dmax} node: " + ", ".join(f"{k}: {v} nodes" for k, v in sorted(kc.items())))
    nfront = sum(d["n_children"] for d in last)
    pop = sum(d["n_children"] for d in by.get(Dmax - 1, [])) if Dmax >= 1 else 1
    if cfg.get("sample_last") or len(last) < pop:
        mk = nfront / len(last) if last else 0.0
        sd = statistics.pstdev([d["n_children"] for d in last]) if len(last) > 1 else 0.0
        fpc = ((pop - len(last)) / max(1, pop - 1)) ** 0.5
        w(f"**Depth-{Dmax} sample: {len(last)} of {pop} nodes, {nfront:,} children; expected frontier at depth "
          f"{Dmax + 1}: {pop * mk:,.0f} nodes (1 SE {pop * sd / max(1, len(last)) ** 0.5 * fpc:,.0f})** "
          f"[e: {pop} x the sample's mean children]")
    else:
        w(f"**Frontier at depth {Dmax + 1}: {nfront:,} nodes** (Burnside orbit count 29,580 for wil50 at depth 3, "
          f"roadmap/tree_size.md section 1, before any fathoming)")
    root = docs.get("r")
    if root and root.get("forced_root_facility") is not None:
        w(f"Root line forced: facility {root['forced_root_facility']} in the 15 J1 cells (the rule's line was "
          f"{root['rule_line']}; out.rule.state)")
    w("")
    # the frontier's inherited bounds
    fr = [f for d in last for f in d.get("frontier", [])]
    if fr and cutoff:
        sh = [cutoff - f["bound"] for f in fr]
        psh = [cutoff - f["pred"] for f in fr]
        w(f"Frontier nodes' inherited bound (the parent's best bound), shortfall: min {min(sh):.1f}, 10% {q(sh, .1):.1f}, "
          f"median {q(sh, .5):.1f}, 90% {q(sh, .9):.1f}, max {max(sh):.1f}. QPB predicted bound (parent's QPB + U_ij), "
          f"shortfall: min {min(psh):.1f}, median {q(psh, .5):.1f}, max {max(psh):.1f}.")
        kept = sum(1 for d in last if d.get("retained"))
        hs = sorted(f["hash"] for f in fr)
        h400 = hs[399] if len(hs) >= 400 else None
        need = sum(1 for d in last if not d.get("retained") and any(f["hash"] <= (h400 or 1) for f in d["frontier"]))
        w(f"Parents kept with their DNN state for stage B: {kept} of {len(last)} (tau history "
          f"{cfg.get('tau_history', [[cfg['tau'], 0]])}: [tau, tasks done when it was set]); the 400th smallest hash "
          f"is {h400}; parents of the 400 that stage B must rebuild: {need}")
    # sanity
    w("")
    w("Sanity checks:")
    root = docs.get("r")
    if root:
        n1 = root["n_children"]
        rg = None
        w(f"- the root: {n1} children (expected 15 orbits)")
    d1 = by.get(1, [])
    if d1:
        cand = sum(d["stats"]["l_level"][1] for d in d1)
        kept = sum(d["n_children"] for d in d1)
        elim = sum(d["stats"]["e_level"][1] for d in d1)
        w(f"- depth 2 before bounding: {cand} candidate children of the {len(d1)} depth-1 nodes ({kept} kept, {elim} "
          f"eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)")
    os.makedirs(a.out, exist_ok=True)
    txt = "\n".join(L) + "\n"
    open(os.path.join(a.out, "stageA_summary.md"), "w").write(txt)
    if fr:
        with open(os.path.join(a.out, "frontier_stats.tsv"), "w") as f:
            f.write("parent\tparent_spec\tparent_line\tparent_first_lb\tparent_status\tchildren\tmin_child_pred\t"
                    "max_child_pred\n")
            for d in sorted(last, key=lambda d: [int(x) if x.isdigit() else -1 for x in d["id"].split(".")]):
                t = d["trace"][0] if d["trace"] else None
                pr = [x["pred"] for x in d.get("frontier", [])]
                f.write(f"{d['id']}\t{d['spec']}\t{d['line']}\t{t['lb'] if t else ''}\t{t['status'] if t else 'QPB'}\t"
                        f"{d['n_children']}\t{min(pr) if pr else ''}\t{max(pr) if pr else ''}\n")
    print(txt)


if __name__ == "__main__":
    main()
