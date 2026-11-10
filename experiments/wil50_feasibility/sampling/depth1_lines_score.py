#!/usr/bin/env python3
"""Scores for depth1_lines.py's runs (results/depth1_lines.md).

  depth1_lines_score.py [--run runs/depth1_lines] [--results results] [--out-dir results/depth1_lines]

Per depth-1 node and candidate line, over the line's 49 depth-2 children, each
bounded to its first DNN solve (no cut phase):
  - fathomed: first solve fathomed (or QPB did, with no DNN solve);
  - the shortfall histogram in stage A's bins (summarize_a.py, compare_lines.py);
  - beyond the zone: first-solve shortfall above cut.zone x 48,816 = 244.08,
    where production runs no cut phase and the children start without a pool;
  - the lines cost, sum over the children of 48^(max(gap, 0)/130);
  - the expected subtree E, if the model is fitted (below): sum over the
    children of E(s), the expected total of a depth-2 node's depth-3 children's
    anchor-model sizes (tree_estimate.size_anchor) given its first-solve
    shortfall s:
        E(s) = 0                     when the first solve fathoms,
        E(s) = q(s) x 48 x m(s)      otherwise,
    with q(s) the fraction of stage A's 635 depth-2 nodes in s's bin (by their
    first-solve shortfall) that the production cut phase did not fathom, and
    m(s) the mean anchor size of stage B's sampled depth-3 nodes whose parents'
    first-solve shortfall is in the bin (the converged shortfall where stage B
    converged the node, else the give-up one; 1 for a fathomed node). 95%
    bootstrap: stage A's 635 and stage B's 400 resampled (4,000 times); the
    line's own 49 children are the whole line, not a sample, so they are not
    resampled. The lines cost and the counts are exact over the 49 children
    and get no interval.
Also: the rewrite check (rule 9's rewritten children against its children
straight from the solver, child by child), HEAD's first solves of rule 9's
children against stage A's (a667727), and CPU.

Writes OUT_DIR/summary.json (every child), OUT_DIR/tables.md, OUT_DIR/model.md,
and per node and line the concatenated trace lines and solver logs (gzipped).
"""
import argparse
import gzip
import json
import math
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import D, parse_trace  # noqa: E402
from tree_estimate import size_anchor  # noqa: E402

ZONE = 0.005 * 48816                      # cut.zone x |incumbent|: the phase runs at shortfall <= 244.08
BINS = [-1e9, 0, 50, 100, 150, 200, 244, 300, 350, 400, 450, 500, 600, 1e9]
NAMES = ["fathomed", "0-50", "50-100", "100-150", "150-200", "200-244", "244-300", "300-350", "350-400",
         "400-450", "450-500", "500-600", ">= 600"]
# the model's bins of the first-solve shortfall (first-solve fathoms apart); the
# last bin is open: stage B has 14 nodes under parents 400-458 short and none beyond
MBINS = [(0, 100), (100, 150), (150, 200), (200, ZONE), (ZONE, 300), (300, 350), (350, 400), (400, 1e9)]
K3 = 48                                   # depth-3 children of a depth-2 node outside the rows-only cells


def short_of(r):
    """first-solve shortfall, -1 for a fathomed node (the DNN's or QPB's)"""
    t = r["trace"]
    if t is None:
        return -1.0 if r["children"] == 0 else None
    return -1.0 if t["fathomed"] else t["short"]


def hist(ss):
    cnt = [0] * len(NAMES)
    for s in ss:
        for k in range(len(NAMES)):
            if (k == 0 and s <= 0) or (k > 0 and BINS[k] < s <= BINS[k + 1]):
                cnt[k] += 1
                break
    return cnt


def lines_cost(ss):
    return sum(48 ** (max(s, 0.0) / 130) for s in ss)


def mbin(s):
    for k, (lo, hi) in enumerate(MBINS):
        if lo < s <= hi or (k == 0 and s <= hi):
            return k
    return len(MBINS) - 1


# ---------------------------------------------------------------- the model's data
def stage_data(results):
    """stage A's depth-2 nodes [(first-solve shortfall, fathomed in the end)] and
    stage B's sampled nodes [(parent's first-solve shortfall, anchor size)]"""
    tr = [r for r in parse_trace(os.path.join(results, "stageA_trace.txt")) if r["depth"] == 2]
    A = [(r["cutoff"] - r["lb"], r["fathomed"] == 1 or r["short"] <= 0) for r in tr]
    first = {r["fixes"]: r["cutoff"] - r["lb"] for r in tr}
    fs = {}
    for line in open(os.path.join(results, "frontier_stats.tsv")):
        t = line.rstrip("\n").split("\t")
        if t[0] != "parent":
            fs[t[0]] = t[1]
    B = []
    for line in open(os.path.join(results, "stageB", "summary.md")):
        if not line.startswith("| ") or line.startswith("| node"):
            continue
        t = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(t) < 15:
            continue
        try:
            g = float(t[4])
        except ValueError:
            continue
        s = -1.0 if t[5] == "fathomed" else (float(t[12]) if t[12] else g)
        B.append((first[fs[t[1]]], size_anchor(s), t[0]))
    return A, B


def fit(A, B):
    """per model bin: q (not fathomed after the phase, stage A) and m (mean anchor
    size of the children, stage B); None where a bin has no data"""
    q, m = [], []
    for k in range(len(MBINS)):
        a = [f for s, f in A if s > 0 and mbin(s) == k]
        b = [z for s, z, _ in B if mbin(s) == k]
        q.append(None if not a else sum(1 for f in a if not f) / len(a))
        m.append(None if not b else sum(b) / len(b))
    return q, m


def E_of(s, q, m):
    if s <= 0:
        return 0.0
    k = mbin(s)
    if q[k] == 0:
        return 0.0
    if q[k] is None or m[k] is None:
        # no data in the bin: the nearest bin below with both
        j = k
        while j > 0 and (q[j] is None or m[j] is None):
            j -= 1
        return (q[j] or 0.0) * K3 * (m[j] or 0.0)
    return q[k] * K3 * m[k]


def centers(A):
    """each model bin's mean first-solve shortfall over stage A's nodes in it"""
    c = []
    for k in range(len(MBINS)):
        v = [x for x, _ in A if x > 0 and mbin(x) == k]
        c.append(sum(v) / len(v) if v else None)
    return c


def E_loglin(s, q, m, cen):
    """E(s) log-linear in s between the bins' centres (their mean first-solve
    shortfall), constant beyond the last centre; 0 at or below 100 (stage A: every
    node there fathomed after the phase) and for a fathomed node; between 100 and
    the first nonzero centre, that bin's own value"""
    if s <= 0 or s <= MBINS[0][1]:
        return 0.0
    pts = [(cen[k], E_of(cen[k], q, m)) for k in range(1, len(MBINS)) if cen[k] is not None]
    pts = [(x, v) for x, v in pts if v > 0]
    if not pts:
        return 0.0
    if s <= pts[0][0]:
        return E_of(s, q, m)
    for (x0, v0), (x1, v1) in zip(pts, pts[1:]):
        if s <= x1:
            return math.exp(math.log(v0) + (math.log(v1) - math.log(v0)) * (s - x0) / (x1 - x0))
    return pts[-1][1]


def fmt(x):
    return f"{x:.3g}" if abs(x) < 1e4 else f"{x:.2e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=os.path.join(D, "runs", "depth1_lines"))
    ap.add_argument("--results", default=os.path.join(D, "results"))
    ap.add_argument("--out-dir", default=os.path.join(D, "results", "depth1_lines"))
    ap.add_argument("--boot", type=int, default=4000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    prep = json.load(open(os.path.join(a.run, "prep.json")))
    cfg = json.load(open(os.path.join(a.run, "config.json")))
    root = json.load(open(os.path.join(a.run, "root", "done.json")))
    # every child's result
    res = {}
    for node in prep["nodes"]:
        nd = os.path.join(a.run, "kids", node)
        if not os.path.isdir(nd):
            continue
        for tag in sorted(os.listdir(nd)):
            for k in sorted(os.listdir(os.path.join(nd, tag))):
                p = os.path.join(nd, tag, k, "result.json")
                if os.path.exists(p):
                    res.setdefault((node, tag), []).append(json.load(open(p)))
    A, B = stage_data(a.results)
    q, m = fit(A, B)
    stageA = {r["fixes"]: r for r in parse_trace(os.path.join(a.results, "stageA_trace.txt"))}

    # settings warnings: only cut.maxdepth and par.hash may differ from the node's record
    bad_warn = []
    for (node, tag), rows in res.items():
        for r in rows:
            w = r.get("settings_warning") or ""
            mm = re.search(r"written with settings \[(\S+) ([^\]]*)\], this run has \[(\S+) ([^\]]*)\]", w)
            if not mm:
                bad_warn.append((r["key"], "no warning parsed"))
                continue
            f1 = dict(x.split("=", 1) for x in mm.group(2).split())
            f2 = dict(x.split("=", 1) for x in mm.group(4).split())
            diff = {k for k in set(f1) | set(f2) if f1.get(k) != f2.get(k)}
            if diff - {"cut.maxdepth", "par.hash"}:
                bad_warn.append((r["key"], sorted(diff)))

    L = []
    # ------------------------------------------------ setup
    L += ["## Setup", "",
          f"Binary: `{cfg['binaries']['bin']['sha256'][:12]}` (commit {cfg['commit'][:7]}, `latest/` unmodified: "
          f"{not cfg['latest_modified']}); diagnostic binary `{cfg['binaries']['bin_diag']['sha256'][:12]}` (the same "
          f"sources plus the QAP_XBAR_DEBUG=2 print). Root: QPB bound {root['root_bound']:.3f}, line {root['line']}, "
          f"{len(root['children'])} children; CPU {root['stats']['times'][0]:.1f} s.", "",
          "| node | DNN bound (HEAD) | shortfall | iterations | CPU s | rule 9's line | lookahead score row / col | "
          "rule 9's line's rank by rule 8's score | diagnostic build: trace, children, warm start |",
          "|---|---|---|---|---|---|---|---|---|"]
    for node, nd in prep["nodes"].items():
        t = nd["trace"]
        x9 = nd["rule9_lookahead"]
        L.append(f"| {node} `{nd['spec']}` | {t['lb']:.4f} | {t['short']:.1f} | {t['iters']} | {nd['cpu']['bin']:.1f} | "
                 f"{nd['rule9_line']} | {x9[0]} {x9[1]}: {float(x9[2]):,.1f} / {float(x9[3]):,.1f} | "
                 f"{nd['rule9_rank'] + 1} | {'equal' if nd['diag_same_trace'] else 'DIFFER'}, "
                 f"{'equal' if nd['diag_same_children'] else 'DIFFER'}, {'equal' if nd['diag_same_warm'] else 'DIFFER'} |")
    L += ["", "Rule 8's score (sum of squared DNN marginals) of the ten best lines, rows and columns ranked together "
          "(the candidates in bold):", ""]
    for node, nd in prep["nodes"].items():
        cands = set(nd["candidates"])
        L.append(f"- {node}: " + ", ".join((f"**{r['line']}** {r['sumsq']:.4f}" if r["line"] in cands
                                           else f"{r['line']} {r['sumsq']:.4f}") for r in nd["rank"][:10]))
    if bad_warn:
        L += ["", f"**Settings warnings other than cut.maxdepth/par.hash: {bad_warn}**"]
    else:
        n_w = sum(len(v) for v in res.values())
        L += ["", f"Every child run ({n_w}) resumed with `--allow-settings-change`; its warning lists cut.maxdepth and "
              "par.hash as the only differences from the node's settings record."]

    # ------------------------------------------------ the rewrite check
    L += ["", "## The rewrite check: rule 9's children rewritten against straight from the solver", "",
          "| node | children | identical first lb, iterations, status | largest abs. lb difference | "
          "same depth-2 line | predicted bound: production (min-max) / rewrite |", "|---|---|---|---|---|---|"]
    check = {}
    for node in prep["nodes"]:
        rw = {r["spec"]: r for r in res.get((node, "rule9"), [])}
        pr = {r["spec"]: r for r in res.get((node, "prod"), [])}
        both = sorted(set(rw) & set(pr))
        same = sum(1 for s in both if rw[s]["trace"] and pr[s]["trace"] and
                   (rw[s]["trace"]["lb"], rw[s]["trace"]["iters"], rw[s]["trace"]["status"]) ==
                   (pr[s]["trace"]["lb"], pr[s]["trace"]["iters"], pr[s]["trace"]["status"]))
        dl = [abs(rw[s]["trace"]["lb"] - pr[s]["trace"]["lb"]) for s in both if rw[s]["trace"] and pr[s]["trace"]]
        sl = sum(1 for s in both if rw[s]["child_line"] == pr[s]["child_line"])
        preds = [pr[s]["predicted"] for s in both]
        check[node] = {"n": len(both), "identical": same, "max_dlb": max(dl) if dl else None,
                       "diffs": [(s, rw[s]["trace"]["lb"], pr[s]["trace"]["lb"], rw[s]["trace"]["iters"],
                                  pr[s]["trace"]["iters"]) for s in both
                                 if rw[s]["trace"] and pr[s]["trace"] and rw[s]["trace"]["lb"] != pr[s]["trace"]["lb"]]}
        if both:
            L.append(f"| {node} | {len(both)} | {same} | {max(dl) if dl else 0:.4f} | {sl} of {len(both)} | "
                     f"{min(preds):.1f}-{max(preds):.1f} / {rw[both[0]]['inherited']:.4f} (inherited) |")
    for node, c in check.items():
        for d in c["diffs"][:10]:
            L.append(f"- {node} `{d[0]}`: rewrite {d[1]:.4f} ({d[3]} it), production {d[2]:.4f} ({d[4]} it)")

    # ------------------------------------------------ HEAD against stage A on rule 9's children
    L += ["", "## HEAD against stage A (a667727) on rule 9's children (production predicted bounds)", "",
          "| node | children matched | identical first lb | median / largest abs. difference | HEAD higher | "
          "same status | stage A fathomed in the first solve / HEAD |", "|---|---|---|---|---|---|---|"]
    for node in prep["nodes"]:
        pr = res.get((node, "prod"), [])
        pairs = [(r, stageA.get(r["spec"])) for r in pr]
        pairs = [(r, s) for r, s in pairs if s is not None and r["trace"]]
        if not pairs:
            continue
        d = sorted(r["trace"]["lb"] - s["lb"] for r, s in pairs)
        L.append(f"| {node} | {len(pairs)} | {sum(1 for x in d if x == 0)} | {sorted(abs(x) for x in d)[len(d) // 2]:.3f} / "
                 f"{max(abs(x) for x in d):.3f} | {sum(1 for x in d if x > 0)} | "
                 f"{sum(1 for r, s in pairs if r['trace']['status'] == s['status'])} | "
                 f"{sum(1 for r, s in pairs if s['status'] == 'fathomed')} / "
                 f"{sum(1 for r, s in pairs if r['trace']['status'] == 'fathomed')} |")

    # ------------------------------------------------ the lines
    rng = random.Random(20261109)
    boots = []
    nA, nB = len(A), len(B)
    for _ in range(a.boot):
        Ab = [A[rng.randrange(nA)] for _ in range(nA)]
        Bb = [B[rng.randrange(nB)] for _ in range(nB)]
        boots.append(fit(Ab, Bb))
    cen = centers(A)
    pct = lambda xs, p: sorted(xs)[min(len(xs) - 1, int(p * len(xs)))]
    summary = {"prep": prep, "root": {"bound": root["root_bound"], "line": root["line"]},
               "model": {"q": q, "m": m, "bins": MBINS, "centers": cen}, "rewrite_check": check, "lines": {},
               "children": []}
    L += ["", "## Per node and line (49 depth-2 children each, first DNN solve, no cut phase)", "",
          "Predicted bound of every child: the inherited bound (the node's DNN bound), for every line, rule 9's "
          "rewritten included; \"prod\" is rule 9's children straight from the solver (production's predicted "
          "bounds), the rewrite check. Ratios are rule 9's score over the line's: above 1, the line is better. "
          "E: the expected-subtree score (model.md) with the step model (bins) and the log-linear one (between "
          "the bins' mean shortfalls), 95% bootstrap over stages A and B.", "",
          "| node | line | rule 8 score (rank) | children | fathomed | in the zone, not fathomed | beyond 244 | "
          "fathomed after a phase, modelled (95%) | "
          "shortfall min / median / max | lines cost 48^(gap/130) | rule 9 / line | top-3 share | E bins (95%) | "
          "rule 9 / line (95%) | E log-linear | rule 9 / line (95%) | CPU s (mean iterations) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]

    def fath_after(ss, qq):
        """the first-solve fathoms plus, for each child in the zone, the fraction of
        stage A's depth-2 nodes in its bin that the production cut phase fathomed"""
        return sum(1.0 if x <= 0 else (1.0 - qq[mbin(x)] if x <= ZONE and qq[mbin(x)] is not None else 0.0)
                   for x in ss)
    H = ["", "Shortfall histogram (first solve; stage A's bins; the zone ends at 244.08):", "",
         "| node | line | " + " | ".join(NAMES) + " |", "|---|---|" + "---|" * len(NAMES)]

    def ratio(num, den):
        return num / den if den > 0 else (1.0 if num == 0 else math.inf)

    for node, nd in prep["nodes"].items():
        rank = {r["line"]: i + 1 for i, r in enumerate(nd["rank"])}
        score = {r["line"]: r["sumsq"] for r in nd["rank"]}
        tags = ["rule9"] + [c.replace(" ", "") for c in nd["candidates"][1:]] + ["prod"]
        lines = {}
        for tag in tags:
            rows = res.get((node, tag), [])
            if rows:
                lines[tag] = {"rows": rows, "ss": [short_of(r) for r in rows]}
        base = lines.get("rule9")
        for tag in tags:
            if tag not in lines:
                continue
            rows, ss = lines[tag]["rows"], lines[tag]["ss"]
            line = nd["rule9_line"] if tag in ("rule9", "prod") else tag.replace("row", "row ").replace("col", "col ")
            fath = sum(1 for s in ss if s <= 0)
            inzone = sum(1 for s in ss if 0 < s <= ZONE)
            beyond = sum(1 for s in ss if s > ZONE)
            lc = lines_cost(ss)
            top3 = sum(sorted((48 ** (max(s, 0.0) / 130) for s in ss), reverse=True)[:3]) / lc
            out = {}
            cells = []
            for var in ("bins", "loglin"):
                if var == "bins":
                    f = lambda x, qq, mm: E_of(x, qq, mm)
                else:
                    f = lambda x, qq, mm: E_loglin(x, qq, mm, cen)
                E = sum(f(x, q, m) for x in ss)
                Eb = [sum(f(x, qq, mm) for x in ss) for qq, mm in boots]
                lo, hi = pct(Eb, 0.025), pct(Eb, 0.975)
                if base and tag != "rule9":
                    rb = [ratio(sum(f(x, qq, mm) for x in base["ss"]), sum(f(x, qq, mm) for x in ss)) for qq, mm in boots]
                    r0 = ratio(sum(f(x, q, m) for x in base["ss"]), E)
                    rr = [r0, pct(rb, 0.025), pct(rb, 0.975)]
                    rat = f"{rr[0]:.2f} ({rr[1]:.2f}-{rr[2]:.2f})"
                else:
                    rr, rat = None, "1"
                out[var] = {"E": E, "E_95": [lo, hi], "rule9_over_line": rr}
                cells.append((f"{fmt(E)} ({fmt(lo)}-{fmt(hi)})" if var == "bins" else fmt(E), rat))
            lcr = ratio(lines_cost(base["ss"]), lc) if base else None
            miss = sorted(x for x in ss if x > 0)
            cpu = sum(r["stats"]["times"][0] for r in rows)
            its = sum(r["trace"]["iters"] for r in rows if r["trace"]) / max(1, sum(1 for r in rows if r["trace"]))
            name = {"rule9": f"{line} (rule 9)", "prod": f"{line} (rule 9, prod)"}.get(tag, line)
            lcr_s = "-" if lcr is None else f"{lcr:.2f}"
            fa = fath_after(ss, q)
            fab = [fath_after(ss, qq) for qq, _ in boots]
            out["fathomed_after_phase_modelled"] = [fa, pct(fab, 0.025), pct(fab, 0.975)]
            L.append(f"| {node} | {name} | {score.get(line, 0):.4f} ({rank.get(line)}) | {len(rows)} | {fath} | {inzone} | "
                     f"{beyond} | {fa:.1f} ({pct(fab, 0.025):.1f}-{pct(fab, 0.975):.1f}) | "
                     + (f"{miss[0]:.1f} / {miss[len(miss) // 2]:.1f} / {miss[-1]:.1f}" if miss else "-")
                     + f" | {fmt(lc)} | {lcr_s} | {top3:.2f} | {cells[0][0]} | {cells[0][1]} | "
                     f"{cells[1][0]} | {cells[1][1]} | {cpu:,.0f} ({its:.0f}) |")
            H.append(f"| {node} | {name} | " + " | ".join(str(x) for x in hist(ss)) + " |")
            summary["lines"][f"{node}/{tag}"] = {
                "line": line, "children": len(rows), "fathomed": fath, "in_zone": inzone, "beyond_zone": beyond,
                "lines_cost": lc, "lines_cost_rule9_over_line": lcr, "top3_share": top3, "E": out, "cpu": cpu,
                "hist": hist(ss)}
            for r, x in zip(rows, ss):
                t = r["trace"] or {}
                summary["children"].append({"node": node, "line": tag, "spec": r["spec"], "short": x,
                                            "lb": t.get("lb"), "status": t.get("status", "qpb-fathomed"),
                                            "iters": t.get("iters"), "dnn_cpu": t.get("cpu"),
                                            "cpu": r["stats"]["times"][0], "wall": r["wall"],
                                            "children": r["children"], "depth2_line": r["child_line"],
                                            "inherited": r["inherited"], "predicted": r["predicted"]})
    L += H
    # ------------------------------------------------ CPU
    allr = [r for v in res.values() for r in v]
    L += ["", f"CPU: {len(allr)} child tasks, {sum(r['stats']['times'][0] for r in allr):,.0f} CPU s "
          f"(first solves {sum((r['trace'] or {}).get('cpu', 0) for r in allr):,.0f} s), wall per task "
          f"{sum(r['wall'] for r in allr) / max(1, len(allr)):.0f} s on average; prep (root and 4 depth-1 runs) "
          f"{root['stats']['times'][0] + sum(sum(nd['cpu'].values()) for nd in prep['nodes'].values()):,.0f} CPU s."]
    open(os.path.join(a.out_dir, "tables.md"), "w").write("\n".join(L) + "\n")

    # ------------------------------------------------ the model
    M = ["# The expected-subtree model (depth1_lines_score.py)", "",
         "E(s) for a depth-2 node of first-solve shortfall s (no cut phase): 0 when the first solve fathoms, else "
         "q x 48 x m, q = the fraction of stage A's depth-2 nodes in the bin (by first-solve shortfall) the production "
         "cut phase did not fathom, m = the mean anchor-model size (tree_estimate.size_anchor) of stage B's sampled "
         "depth-3 nodes whose parents are in the bin. 95% bootstrap intervals over stage A's 635 and stage B's 400.", "",
         "| first-solve shortfall | stage A nodes | their mean shortfall | q (not fathomed after the phase) | stage B nodes | "
         "m (mean child size) | largest child size | E at the mean shortfall (95%) |", "|---|---|---|---|---|---|---|---|"]
    cen = centers(A)
    for k, (lo_, hi_) in enumerate(MBINS):
        a_n = sum(1 for s, f in A if s > 0 and mbin(s) == k)
        b_v = [z for s, z, _ in B if mbin(s) == k]
        e = [(qq[k] or 0) * K3 * (mm[k] or 0) for qq, mm in boots if qq[k] is not None and mm[k] is not None]
        e.sort()
        ci = f" ({fmt(e[int(0.025 * len(e))])}-{fmt(e[int(0.975 * len(e)) - 1])})" if e else ""
        M.append(f"| {lo_:.0f}-{'' if hi_ > 1e8 else f'{hi_:.2f}'} | {a_n} | "
                 + (f"{cen[k]:.1f}" if cen[k] is not None else "-") + " | "
                 + (f"{q[k]:.3f}" if q[k] is not None else "-") + f" | {len(b_v)} | "
                 + (f"{m[k]:.2f}" if m[k] is not None else "-") + f" | {max(b_v) if b_v else 0:.1f} | "
                 + (f"{fmt(E_of(cen[k], q, m))}{ci}" if q[k] is not None else "-") + " |")
    M += ["", f"Stage A depth-2 nodes fathomed in the first solve: {sum(1 for s, f in A if s <= 0)} (E = 0)."]
    open(os.path.join(a.out_dir, "model.md"), "w").write("\n".join(M) + "\n")
    with open(os.path.join(a.out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    # ------------------------------------------------ traces and logs per node and line
    for (node, tag), rows in res.items():
        tdir = os.path.join(a.out_dir, "traces")
        ldir = os.path.join(a.out_dir, "logs")
        os.makedirs(tdir, exist_ok=True)
        os.makedirs(ldir, exist_ok=True)
        with open(os.path.join(tdir, f"{node}_{tag}.txt"), "w") as ft, \
                gzip.open(os.path.join(ldir, f"{node}_{tag}.log.gz"), "wt") as fl:
            for r in sorted(rows, key=lambda r: r["k"]):
                wd = os.path.join(a.run, "kids", r["key"])
                ft.write(open(os.path.join(wd, "trace.txt")).read())
                fl.write(f"===== {r['key']} [{r['spec']}] =====\n")
                fl.write(open(os.path.join(wd, "solver.log"), errors="replace").read())
    for node in prep["nodes"]:
        for b in ("bin", "bin_diag"):
            wd = os.path.join(a.run, "d1", node, b)
            ddir = os.path.join(a.out_dir, "depth1")
            os.makedirs(ddir, exist_ok=True)
            for fn in ("trace.txt", "solver.log"):
                open(os.path.join(ddir, f"{node}_{b}_{fn}"), "w").write(open(os.path.join(wd, fn), errors="replace").read())
    ddir = os.path.join(a.out_dir, "depth1")
    open(os.path.join(ddir, "root_solver.log"), "w").write(open(os.path.join(a.run, "root", "solver.log")).read())
    print("\n".join(L))
    print("\n".join(M))


if __name__ == "__main__":
    main()
