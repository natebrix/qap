#!/usr/bin/env python3
"""Does the marginal score of a branching line predict its subtree? (tho30, Nov 9 2026)

At a depth-1 node of tho30 under the record recipe (../cutd2_cp.par), the
candidate lines are rule 9's choice, the top line of the other side (row vs
column) and the next lines by marginal score (sum of squared DNN marginals,
QAP_XBAR_DEBUG=2). For each line, every child gets the same two treatments:
  (a) score: the child's first DNN solve from a cold start (dnntest, give-up 300,
      no cut phase); its shortfall is cutoff - lb (negative when it fathoms)
  (b) truth: the child's subtree solved to completion with -F (the node's fix
      plus the child's) under the recipe; the line's truth is the sum over its
      children plus 1, in nodes and in CPU seconds.
The report compares the ranking of the lines by (a) with the ranking by (b).

  python3 line_truth.py --node 5:2 [--node 5:1 ...] [--lines 8] [--jobs 8] [--out results]

Reads qaplib/tho30.dat and the par files of this directory; writes
results/<node>.json (every child) and results/summary.md.
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LATEST = os.path.join(ROOT, "latest")
BNBDRV = os.path.join(LATEST, "bnbdrv")
DNNTEST = os.path.join(LATEST, "dnntest")
DAT = os.path.join(ROOT, "qaplib", "tho30.dat")
PAR_R8 = os.path.join(HERE, "..", "cutd2_cp.par")        # the record recipe (rule 8, cut phase through depth 2)
PAR_R9 = os.path.join(HERE, "cutd2_cp_r9.par")           # the same with rule 9, for its line choice
UB = 149936
N = 30
ENV = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")


def run(cmd, env=None, cwd=LATEST):
    out = subprocess.run(cmd, capture_output=True, text=True, env=env or ENV, cwd=cwd)
    return out.returncode, out.stdout, out.stderr


def parse_fix(s):
    i, j = s.split(":")
    return int(i), int(j)


def ranking(node):
    """The node's lines by marginal score and rule 9's choice (the solver's own prints)."""
    env = dict(ENV, QAP_XBAR_DEBUG="2")
    rc, out, err = run([BNBDRV, "-p", PAR_R9, "-u", str(UB), "-F", node, "-d", str(node.count(":") + 1),
                        "-x", "1:7:1:200:1", DAT], env=env)
    if rc != 0:
        sys.exit("ranking run failed (rc=%d):\n%s\n%s" % (rc, out[-2000:], err[-2000:]))
    depth = node.count(",") + 1
    lines, choice = None, None
    for ln in err.splitlines():
        m = re.match(r"xbar-lines: depth (\d+) (.*)", ln)
        if m and int(m.group(1)) == depth and lines is None:
            lines = [(t.split(":")[0], float(t.split(":")[1])) for t in m.group(2).split()]
        m = re.match(r"xbar9: depth (\d+) (row|col) (\d+)", ln)
        if m and int(m.group(1)) == depth and choice is None:
            choice = m.group(2) + m.group(3)
    if lines is None or choice is None:
        sys.exit("no ranking for %s in the solver's output:\n%s" % (node, err[-3000:]))
    return lines, choice


def children(node, line):
    """(fixes string, child label) for every child of the line at the node."""
    fixes = [parse_fix(t) for t in node.split(",")]
    used_i = {i for i, _ in fixes}
    used_j = {j for _, j in fixes}
    side, idx = line[:3], int(line[3:])
    kids = []
    if side == "row":
        assert idx not in used_i, line
        for l in range(N):
            if l not in used_j:
                kids.append(("%s,%d:%d" % (node, idx, l), "%d:%d" % (idx, l)))
    else:
        assert idx not in used_j, line
        for k in range(N):
            if k not in used_i:
                kids.append(("%s,%d:%d" % (node, k, idx), "%d:%d" % (k, idx)))
    return kids


def score_child(fixes):
    """(a): the first DNN solve of the child, cold, give-up 300."""
    rc, out, err = run([DNNTEST, "-p", PAR_R8, "-x", fixes, "-u", str(UB), "-H", "300", "-i", "4000", DAT])
    if rc != 0:
        return {"error": "dnntest rc=%d: %s" % (rc, (out + err)[-500:])}
    for ln in out.splitlines():
        if ln.startswith("RESULT"):
            d = dict(kv.split("=", 1) for kv in ln.split()[1:])
            return {"lb": float(d["lb"]), "status": d["status"], "iters": int(d["iters"]),
                    "cpu": float(d["cpu"]), "fathomed": d.get("fathomed") == "1"}
    return {"error": "no RESULT line"}


def truth_child(fixes):
    """(b): the child's subtree to completion under the recipe."""
    t0 = time.time()
    rc, out, err = run([BNBDRV, "-p", PAR_R8, "-u", str(UB), "-F", fixes, DAT])
    if rc != 0:
        return {"error": "bnbdrv rc=%d: %s" % (rc, (out + err)[-500:])}
    r = {"wall": time.time() - t0}
    m = re.search(r"Total Time: ([0-9.eE+-]+) seconds", out)
    if m: r["cpu"] = float(m.group(1))
    m = re.search(r"Nodes: (\d+) / (\d+) / (\d+)", out)
    if m: r["nodes"], r["fathomed"], r["other"] = int(m.group(1)), int(m.group(2)), int(m.group(3))
    m = re.search(r"DNN bound: (\d+) nodes, (\d+) fathomed", out)
    if m: r["dnn_nodes"], r["dnn_fathomed"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"cut phase: (\d+) nodes, (\d+) fathomed", out)
    if m: r["cut_nodes"], r["cut_fathomed"] = int(m.group(1)), int(m.group(2))
    if "nodes" not in r:
        r["error"] = "no totals: " + out[-500:]
    return r


def candidates(lines, choice, k):
    """rule 9's line, the top line of the other side, then by score to k lines."""
    top_row = next(l for l, _ in lines if l.startswith("row"))
    top_col = next(l for l, _ in lines if l.startswith("col"))
    cand = [choice]
    for l in (top_row, top_col):
        if l not in cand: cand.append(l)
    for l, _ in lines:
        if len(cand) >= k: break
        if l not in cand: cand.append(l)
    return cand


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]: j += 1
            for t in range(i, j + 1): r[order[t]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den > 0 else float("nan")


def study(node, k, jobs, outdir, pool):
    path = os.path.join(outdir, node.replace(":", "-").replace(",", "_") + ".json")
    if os.path.exists(path):
        with open(path) as f: rec = json.load(f)
        if rec.get("complete"):
            print("%s: done already (%s)" % (node, path)); return rec
    else:
        rec = {"node": node, "ub": UB, "lines": {}, "complete": False}
    if "ranking" not in rec:
        print("%s: ranking ..." % node, flush=True)
        lines, choice = ranking(node)
        rec["ranking"] = lines; rec["rule9"] = choice
        rec["candidates"] = candidates(lines, choice, k)
        save(path, rec)
    print("%s: rule 9 %s; candidates %s" % (node, rec["rule9"], " ".join(rec["candidates"])), flush=True)
    for line in rec["candidates"]:
        L = rec["lines"].setdefault(line, {"score": dict(rec["ranking"])[line], "children": {}})
        kids = children(node, line)
        todo = [(f, lab) for f, lab in kids if "score" not in L["children"].get(lab, {})]
        if todo:
            print("  %s (a) %d children ..." % (line, len(todo)), flush=True)
            for (f, lab), res in zip(todo, pool.map(lambda t: score_child(t[0]), todo)):
                L["children"].setdefault(lab, {})["score"] = res
            save(path, rec)
        todo = [(f, lab) for f, lab in kids if "truth" not in L["children"].get(lab, {})]
        if todo:
            print("  %s (b) %d subtrees ..." % (line, len(todo)), flush=True)
            for (f, lab), res in zip(todo, pool.map(lambda t: truth_child(t[0]), todo)):
                L["children"][lab]["truth"] = res
            save(path, rec)
        print("  %s: %s" % (line, line_summary(L, rec)), flush=True)
    rec["complete"] = True
    save(path, rec)
    return rec


def line_totals(L, cutoff):
    sc = [c["score"] for c in L["children"].values() if "score" in c and "lb" in c["score"]]
    tr = [c["truth"] for c in L["children"].values() if "truth" in c and "nodes" in c["truth"]]
    short = [cutoff - s["lb"] for s in sc]
    return {"children": len(L["children"]),
            "shortfall_sum": sum(max(0.0, x) for x in short), "shortfall_mean": sum(short) / len(short) if short else float("nan"),
            "fathomed_a": sum(1 for s in sc if s["lb"] > cutoff), "score_cpu": sum(s["cpu"] for s in sc),
            "nodes": 1 + sum(t["nodes"] for t in tr), "cpu": sum(t["cpu"] for t in tr),
            "truth_done": len(tr)}


def cutoff_of(rec):
    return rec["ub"] - 1.95   # even objective values: the solver's fathoming tolerance for tho30


def line_summary(L, rec):
    t = line_totals(L, cutoff_of(rec))
    return ("score %.4f | (a) shortfall sum %.1f, %d of %d children fathom cold | (b) %s nodes, %s CPU s"
            % (L["score"], t["shortfall_sum"], t["fathomed_a"], t["children"],
               t["nodes"] if t["truth_done"] == t["children"] else "..", "%.0f" % t["cpu"] if t["truth_done"] == t["children"] else ".."))


def save(path, rec):
    tmp = path + ".tmp"
    with open(tmp, "w") as f: json.dump(rec, f, indent=1)
    os.replace(tmp, path)


def summary(recs, outdir):
    out = ["# Line score against subtree truth, tho30 depth-1 nodes (%s)\n" % time.strftime("%Y-%m-%d"),
           "Recipe `../cutd2_cp.par`; (a) cold first solve of every child, give-up 300, no cut phase, "
           "shortfall = cutoff - lb summed over the children where positive; (b) every child's subtree "
           "to completion with `-F`, summed plus 1. Ranks: 1 = best (smallest).\n"]
    for rec in recs:
        cutoff = cutoff_of(rec)
        rows = []
        for line in rec["candidates"]:
            L = rec["lines"][line]
            t = line_totals(L, cutoff)
            rows.append((line, L["score"], t))
        by_score = sorted(rows, key=lambda r: -r[1])
        by_a = sorted(rows, key=lambda r: r[2]["shortfall_sum"])
        by_nodes = sorted(rows, key=lambda r: r[2]["nodes"])
        by_cpu = sorted(rows, key=lambda r: r[2]["cpu"])
        rank = lambda lst, line: 1 + [r[0] for r in lst].index(line)
        out.append("## Node %s (rule 9 chose %s)\n" % (rec["node"], rec["rule9"]))
        out.append("| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |")
        out.append("|---|---|---|---|---|---|---|---|---|---|")
        for line, s, t in by_nodes:
            tag = " (rule 9)" if line == rec["rule9"] else ""
            out.append("| %s%s | %.4f | %d | %.1f | %d/%d | %d | %d | %d | %.0f | %d |"
                       % (line, tag, s, rank(by_score, line), t["shortfall_sum"], t["fathomed_a"], t["children"],
                          rank(by_a, line), t["nodes"], rank(by_nodes, line), t["cpu"], rank(by_cpu, line)))
        xs = [r[2]["shortfall_sum"] for r in rows]
        sc = [-r[1] for r in rows]
        ns = [r[2]["nodes"] for r in rows]
        cs = [r[2]["cpu"] for r in rows]
        best = by_nodes[0]; r9 = next(r for r in rows if r[0] == rec["rule9"])
        out.append("")
        out.append("Spearman rank correlation with (b) nodes: (a) shortfall sum %.2f, marginal score %.2f; "
                   "with (b) CPU: (a) %.2f, score %.2f." % (spearman(xs, ns), spearman(sc, ns), spearman(xs, cs), spearman(sc, cs)))
        out.append("Rule 9's line %s: rank %d of %d by (b) nodes, %d by (b) CPU; the best line %s has %.2fx its nodes and %.2fx its CPU."
                   % (rec["rule9"], rank(by_nodes, rec["rule9"]), len(rows), rank(by_cpu, rec["rule9"]), best[0],
                      best[2]["nodes"] / r9[2]["nodes"], best[2]["cpu"] / r9[2]["cpu"] if r9[2]["cpu"] > 0 else float("nan")))
        out.append("")
    with open(os.path.join(outdir, "summary.md"), "w") as f: f.write("\n".join(out))
    print("\n".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--node", action="append", required=True, help="depth-1 node as i:j (repeatable)")
    ap.add_argument("--lines", type=int, default=8, help="candidate lines per node (default 8)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--summary-only", action="store_true", help="rewrite summary.md from the saved results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    recs = []
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for node in a.node:
            if a.summary_only:
                with open(os.path.join(a.out, node.replace(":", "-") + ".json")) as f: recs.append(json.load(f))
            else:
                recs.append(study(node, a.lines, a.jobs, a.out, pool))
    summary(recs, a.out)


if __name__ == "__main__":
    main()
