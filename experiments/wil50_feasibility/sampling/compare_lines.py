#!/usr/bin/env python3
"""The root-line comparison (README.md, "Root lines"): stage A under the rule's
root line (facility 21) against the forced lines (run_samples.sh: 45, 25).

  compare_lines.py [--runs 21=runs/stageA 45=runs/stageA_f45 25=runs/stageA_f25] [--laptop 2d26817]

Per facility: the 15 depth-1 first solves next to the laptop's 1,500-iteration
cold bounds of the same children (tuning/lines/root/children/f<F>_c<cell>.log
at --laptop); for the depth-2 nodes, the seeded sample of 200 (the 200 task ids
with the smallest order hash, the same ids in every run; for the full run of 21
also the whole population): the fathomed fraction, the shortfall histogram, the
expected frontier (635 x the mean children), the lines cost
635 x mean(47^(max(gap, 0)/130)) over the depth-2 nodes with a bootstrap
interval, CPU per node, and at each depth-1 node the line it branched on and the
spread of that line's depth-2 children. Writes results/root_lines.md.
"""
import argparse
import json
import math
import os
import random
import re
import statistics
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import D, R  # noqa: E402
from stage_a import StageA  # noqa: E402

CELLS = [0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 20, 21, 22, 23, 24]
BINS = [-1e9, 0, 50, 100, 150, 200, 244, 300, 350, 400, 450, 500, 600, 1e9]
NAMES = ["fathomed", "0-50", "50-100", "100-150", "150-200", "200-244", "244-300", "300-350", "350-400",
         "400-450", "450-500", "500-600", ">= 600"]


def load(run):
    docs = {}
    for t in os.listdir(os.path.join(run, "tasks")):
        p = os.path.join(run, "tasks", t, "done.json")
        if os.path.exists(p):
            docs[t] = json.load(open(p))
    return docs, json.load(open(os.path.join(run, "config.json")))


def laptop_bounds(rev, F, sub="children"):
    out = {}
    for c in CELLS:
        path = f"experiments/wil50_feasibility/tuning/lines/root/{sub}/f{F}_c{c}.log"
        r = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=R, capture_output=True, text=True)
        m = re.search(r"lb=(\S+).*iters=(\d+)", r.stdout)
        out[c] = (float(m.group(1)), int(m.group(2))) if m else None
    return out


def short(d):
    t = d["trace"][0] if d["trace"] else None
    if t is None or t["fathomed"]:
        return -1.0
    return t["short"]


def hist(nodes):
    cnt = [0] * len(NAMES)
    for d in nodes:
        s = short(d)
        for k in range(len(NAMES)):
            if (k == 0 and s <= 0) or (k > 0 and BINS[k] < s <= BINS[k + 1]):
                cnt[k] += 1
                break
    return cnt


def wilson(k, n, z=1.0):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def lines_cost(nodes, floor=False):
    """47^(max(gap, 0)/130) per node (a fathomed node counts 1); floor: a miss
    counts at least itself and its children (1 + children), which the
    exponential model understates for misses a few units short"""
    out = []
    for d in nodes:
        v = 47 ** (max(short(d), 0.0) / 130)
        if floor and d["n_children"]:
            v = max(v, 1.0 + d["n_children"])
        out.append(v)
    return out


def boot(vals, pop, seed=1, B=4000):
    rng = random.Random(seed)
    n = len(vals)
    ms = sorted(pop * sum(vals[rng.randrange(n)] for _ in range(n)) / n for _ in range(B))
    return ms[int(0.025 * B)], ms[int(0.975 * B) - 1]


def cpu_of(d):
    t = d["trace"][0] if d["trace"] else None
    first = d["stats"]["dnn_time"][0]
    phase = t.get("cut_cpu", 0.0) if t else 0.0
    return first, phase, d["stats"]["times"][0]


def depth2_block(name, nodes, pop, lines):
    n = len(nodes)
    f = sum(1 for d in nodes if d["n_children"] == 0)
    lo, hi = wilson(f, n)
    kids = [d["n_children"] for d in nodes]
    mk = statistics.mean(kids)
    se = statistics.pstdev(kids) / math.sqrt(n) * math.sqrt(max(0, pop - n) / max(1, pop - 1)) if n > 1 else 0
    lc = lines_cost(nodes)
    blo, bhi = boot(lc, pop)
    lf = lines_cost(nodes, floor=True)
    flo, fhi = boot(lf, pop)
    c = [cpu_of(d) for d in nodes]
    ph = sum(1 for d in nodes if d["trace"] and d["trace"][0]["cut"])
    lines.append(f"| {name} | {n} of {pop} | {f} ({f / n:.3f}; 1-sigma {lo:.3f}-{hi:.3f}) | {ph} | "
                 f"{pop * mk:,.0f} ± {pop * se:,.0f} | {pop * statistics.mean(lc):.3g} ({blo:.3g}-{bhi:.3g}) | "
                 f"{pop * statistics.mean(lf):.3g} ({flo:.3g}-{fhi:.3g}) | "
                 f"{statistics.mean(x[0] for x in c):.0f} + {statistics.mean(x[1] for x in c):.0f} = "
                 f"{statistics.mean(x[2] for x in c):.0f} |")
    return hist(nodes)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=[f"21={D}/runs/stageA", f"45={D}/runs/stageA_f45",
                                                  f"25={D}/runs/stageA_f25"])
    ap.add_argument("--laptop", default="2d26817", help="the laptop's 1,500-iteration children")
    ap.add_argument("--laptop300", default="165780c", help="its 300-iteration calibration (children_i300)")
    ap.add_argument("--out", default=os.path.join(D, "results", "root_lines.md"))
    a = ap.parse_args()
    runs = {}
    for kv in a.runs:
        F, path = kv.split("=", 1)
        if os.path.exists(os.path.join(path, "config.json")):
            runs[int(F)] = load(path)
    L = ["# Root lines: facility 21 (the rule's line, full stage A) against forced lines", ""]
    # 1. depth-1 bounds
    L += ["## Depth-1 first solves against the laptop's children", "",
          "Ours: the production first solve (QPB warm start from the root, cold DNN, give-up at 300 iterations "
          "when hopeless; Xeon, OpenBLAS). Laptop: the same children cold with `-i 300` "
          f"(`tuning/lines/root/children_i300/`, rev {a.laptop300}) and `-i 1,500` (`children/`, rev {a.laptop}). "
          "Bound (shortfall to the cutoff 48,814.05).", ""]
    hdr = "| cell |" + "".join(f" f{F} ours, 300 it | f{F} laptop, 300 it | f{F} laptop, 1,500 it | f{F} line, children |"
                                for F in sorted(runs))
    L += [hdr, "|---|" + "---|---|---|---|" * len(runs)]
    lap = {F: laptop_bounds(a.laptop, F) for F in runs}
    lap3 = {F: laptop_bounds(a.laptop300, F, "children_i300") for F in runs}
    rows = {F: {} for F in runs}
    for F, (docs, cfg) in runs.items():
        for k, c in enumerate(CELLS):
            rows[F][c] = docs.get(f"r.{k}")
    fmt = lambda lb: f" {lb[0]:.1f} ({48814.05 - lb[0]:.1f}) |" if lb else " - |"
    for c in CELLS:
        line = f"| {c} |"
        for F in sorted(runs):
            d = rows[F][c]
            t = d["trace"][0] if d and d["trace"] else None
            line += (f" {t['lb']:.1f} ({t['short']:.1f}) |" if t else " - |")
            line += fmt(lap3[F].get(c)) + fmt(lap[F].get(c))
            line += (f" {d['line']}, {d['n_children']} |" if d else " - |")
        L.append(line)
    L += ["", "Lines cost over the 15 children, sum 47^(gap/130), and the ratio to facility 45's:", "",
          "| facility | ours, 300 it | laptop, 300 it | laptop, 1,500 it | largest difference ours - laptop at 300 it |",
          "|---|---|---|---|---|"]
    cost = lambda xs: sum(47 ** (x / 130) for x in xs)
    costs = {}
    for F in sorted(runs):
        ours = [rows[F][c]["trace"][0]["short"] for c in CELLS if rows[F][c] and rows[F][c]["trace"]]
        t3 = [48814.05 - lap3[F][c][0] for c in CELLS if lap3[F].get(c)]
        t15 = [48814.05 - lap[F][c][0] for c in CELLS if lap[F].get(c)]
        diff = max((abs(rows[F][c]["trace"][0]["lb"] - lap3[F][c][0]) for c in CELLS
                    if rows[F][c] and rows[F][c]["trace"] and lap3[F].get(c)), default=float("nan"))
        costs[F] = (cost(ours) if len(ours) == 15 else None, cost(t3) if len(t3) == 15 else None,
                    cost(t15) if len(t15) == 15 else None, diff)
    base = costs.get(45)
    for F in sorted(runs):
        cells = []
        for k in range(3):
            v = costs[F][k]
            r = f", {v / base[k]:.2f}x" if v and base and base[k] else ""
            cells.append(f"{v:.3g}{r}" if v else "-")
        L.append(f"| {F} | " + " | ".join(cells) + f" | {costs[F][3]:.2f} |")
    # 2-6. depth 2
    L += ["", "## Depth 2", "",
          "The sample: the 200 depth-2 task ids with the smallest seeded order hash (`stage_a.py --sample-last`), "
          "the same ids in every run. Lines cost: 635 x mean over the nodes of 47^(max(gap, 0)/130), gap = the "
          "node's shortfall after its first solve and cut phase (a fathomed node counts 1), with a 95% bootstrap "
          "interval (4,000 resamples); expected frontier: 635 x the mean children (± 1 SE).", "",
          "| line | depth-2 nodes | fathomed (fraction) | cut phases | expected frontier | lines cost (95%) | "
          "lines cost, a miss at least 1 + children (95%) | CPU s per node: first solve + phase = all |",
          "|---|---|---|---|---|---|---|---|"]
    hists = {}
    for F in sorted(runs):
        docs, cfg = runs[F]
        d2 = {t: d for t, d in docs.items() if d["depth"] == 2}
        pop = sum(d["n_children"] for t, d in docs.items() if d["depth"] == 1)
        ids = sorted((t for t in d2), key=lambda t: StageA.order_hash(cfg["seed"], t))
        if cfg.get("sample_last"):
            samp = [d2[t] for t in ids]
            hists[f"f{F} sample"] = depth2_block(f"f{F}, sample", samp, pop, L)
        else:
            if len(d2) == pop:
                hists[f"f{F} all"] = depth2_block(f"f{F}, all (population)", list(d2.values()), pop, L)
            allids = sorted((f"r.{k}.{j}" for k, c in enumerate(CELLS) for j in range(rows[F][c]["n_children"])),
                            key=lambda t: StageA.order_hash(cfg["seed"], t))[:200]
            samp = [d2[t] for t in allids if t in d2]
            if len(samp) == 200:
                hists[f"f{F} matched sample"] = depth2_block(f"f{F}, matched sample", samp, pop, L)
    L += ["", "Shortfall histogram of the depth-2 nodes (after the first solve and cut phase; the zone is 244):", "",
          "| nodes | " + " | ".join(NAMES) + " |", "|---|" + "---|" * len(NAMES)]
    for k, h in hists.items():
        L.append(f"| {k} | " + " | ".join(str(x) for x in h) + " |")
    # 7. the depth-1 question
    L += ["", "## The line at each depth-1 node and its children", "",
          "Rows-only nodes (the cell in the grid's middle row, 20-24: the J2 symmetry) branch by QPB row "
          "lookahead (`select_row_lookahead`, no DNN marginals); the others by rule 9's lookahead between the "
          "DNN-marginal row and column (`select_ind_lookahead` with `xbar_lines`). Children: the depth-2 nodes "
          "bounded (all for the full run, the sample's otherwise).", "",
          "| line | depth-1 node | how | line | children | bounded | fathomed | gap min / median / max |",
          "|---|---|---|---|---|---|---|---|"]
    for F in sorted(runs):
        docs, cfg = runs[F]
        for k, c in enumerate(CELLS):
            d = docs.get(f"r.{k}")
            if not d:
                continue
            kids = [x for t, x in docs.items() if x["depth"] == 2 and x["parent"] == f"r.{k}"]
            gs = sorted(short(x) for x in kids if short(x) > 0)
            how = "rows-only QPB lookahead" if c >= 20 else "rule 9 xbar lines"
            sp = f"{gs[0]:.0f} / {statistics.median(gs):.0f} / {gs[-1]:.0f}" if gs else "-"
            L.append(f"| f{F} | `{d['spec']}` | {how} | {d['line']} | {d['n_children']} | {len(kids)} | "
                     f"{sum(1 for x in kids if x['n_children'] == 0)} | {sp} |")
    txt = "\n".join(L) + "\n"
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
