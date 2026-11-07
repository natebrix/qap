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
    w("| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations | DNN CPU s | CPU s (all) | "
      "CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
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
        it = sum(r["iters"] + r.get("cut_iters", 0) for r in tr)
        dcpu = sum(d["stats"]["dnn_time"][0] for d in ds)
        cpu = sum(d["stats"]["times"][0] for d in ds)
        warm = sum(r["warm"] for r in tr)
        ph = [r for r in tr if r["cut"]]
        phases += ph
        kept = sum(d["n_children"] for d in ds)
        elim = sum(d["stats"]["e_level"][depth] for d in ds)
        w(f"| {depth} | {nb} | {nf} | {dn} | {df} | {it:,} | {dcpu:,.0f} | {cpu:,.0f} | {cpu / max(1, nb):,.1f} | "
          f"{warm} | {len(ph)} | {kept:,} | {elim:,} |")
        tot.update(nodes=nb, fathomed=nf, dnn=dn, dnn_f=df, iters=it, dnn_cpu=dcpu, cpu=cpu)
    w(f"| all | {tot['nodes']} | {tot['fathomed']} | {tot['dnn']} | {tot['dnn_f']} | {tot['iters']:,} | "
      f"{tot['dnn_cpu']:,.0f} | {tot['cpu']:,.0f} | | | {len(phases)} | | |")
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
          + ", ".join(f"{k} {v}" for k, v in sorted(cnt.items())))
        if depth < Dmax and ds:
            w("  " + "; ".join(f"`{d['spec'] or 'root'}` {d['line']} ({d['n_children']})" for d in ds))
    w("")
    # children counts at the last depth
    last = by.get(Dmax, [])
    kc = collections.Counter(d["n_children"] for d in last)
    w(f"Children per depth-{Dmax} node: " + ", ".join(f"{k}: {v} nodes" for k, v in sorted(kc.items())))
    nfront = sum(d["n_children"] for d in last)
    w(f"**Frontier at depth {Dmax + 1}: {nfront:,} nodes** (Burnside orbit count 29,580 for wil50 at depth 3, "
      f"roadmap/tree_size.md section 1, before any fathoming)")
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
        w(f"Parents kept with their DNN state for stage B: {kept} of {len(last)}; frontier nodes below the hash "
          f"threshold: {sum(1 for f in fr if f['hash'] < cfg['tau'])}")
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
