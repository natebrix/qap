#!/usr/bin/env python3
"""Score the candidate lines of the line-truth nodes by the signals the branching rules could use, and
compare them with the truth (README.md):

  python3 experiments/line_signals/signals.py

Reads dumps/*.err (dump.sh: U, the bound, the cutoff, the DNN marginals and every line's lookahead score at
the node) and the truth: ../tho30_solve/line_truth/results/*.json and ../tho40_line_truth/results/*.json
(every child's subtree to the end, (b)), the tai60b line study's four lines (../tai60b_lines: branching
children and the rough subtree model). Signals, each oriented so that larger means a better line:
  marg   the marginal score, sum of squared DNN marginals (rule 8; rule 9's shortlist)
  r2     rule 2: children with bound + U_ij >= cutoff, ties by the U sum over the others
  usum   the U sum over the line (rule 1)
  la     rule 4's lookahead score: every child bounded by QPB (rule 9's arbiter)
  zero   children the DNN marginals rule out (x_ij < 1e-3)
  short  (a), the children's cold first DNN solve: minus the shortfall sum (truth-side, costly)
"""
import collections, glob, json, math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)


def parse_dump(path):
    d = {"U": [], "X": []}
    for ln in open(path):
        if ln.startswith("ls-node:"):
            f = ln.split()
            kv = dict(zip(f[1::2], f[2::2]))
            d.update(depth=int(kv["depth"]), m=int(kv["m"]), bound=float(kv["bound"]), inc=float(kv["inc"]),
                     delta=float(kv["delta"]))
        elif ln.startswith("ls-rows:"): d["rows"] = [int(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-cols:"): d["cols"] = [int(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-U "): d["U"].append([float(x) for x in ln.split(":", 1)[1].split()])
        elif ln.startswith("ls-X "): d["X"].append([float(x) for x in ln.split(":", 1)[1].split()])
        elif ln.startswith("ls-LArow:"): d["LArow"] = [float(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-LAcol:"): d["LAcol"] = [float(x) for x in ln.split()[1:]]
        elif ln.startswith("xbar9:") and "rule9" not in d:
            m = re.match(r"xbar9: depth \d+ (row|col) (\d+)", ln)
            d["rule9"] = m.group(1) + m.group(2)
    return d


def line_signals(d):
    """{line: {signal: value}}, lines named row<i> / col<j> by original index"""
    U, X, m = d["U"], d["X"], d["m"]
    cut = d["inc"] - d["delta"]
    out = {}
    for side in ("row", "col"):
        for k in range(m):
            cells = [(U[k][j], X[k][j] if X else None) for j in range(m)] if side == "row" else \
                    [(U[i][k], X[i][k] if X else None) for i in range(m)]
            fath = sum(1 for u, _ in cells if u + d["bound"] >= cut)
            usum_open = sum(u for u, _ in cells if u + d["bound"] < cut)
            name = side + str((d["rows"] if side == "row" else d["cols"])[k])
            s = {"r2": fath * 1e12 + usum_open, "r2_count": fath, "usum": sum(u for u, _ in cells),
                 "la": (d["LArow"] if side == "row" else d["LAcol"])[k]}
            if X:
                s["marg"] = sum(x * x for _, x in cells)
                s["zero"] = sum(1 for _, x in cells if x < 1e-3)
            out[name] = s
    return out


def ranks(vals, larger_better=True):
    """rank 1 = best; ties share the mean rank"""
    order = sorted(range(len(vals)), key=lambda i: -vals[i] if larger_better else vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num / den if den else float("nan")


def truth_tables():
    """[(name, dump file, {line: {"nodes", "cpu", "short"}}, rule 9's line in the study)]"""
    nodes = []
    for inst, pat in (("tho30", os.path.join(EXP, "tho30_solve/line_truth/results/*.json")),
                      ("tho40", os.path.join(EXP, "tho40_line_truth/results/*.json"))):
        for path in sorted(glob.glob(pat)):
            rec = json.load(open(path))
            cutoff = rec["ub"] - (1.95 if inst == "tho40" or rec["ub"] % 2 == 0 else 0.95)
            T = {}
            for line, L in rec["lines"].items():
                kids = L["children"].values()
                if not all("truth" in k and "nodes" in k["truth"] for k in kids):
                    continue
                T[line] = {"nodes": 1 + sum(k["truth"]["nodes"] for k in kids),
                           "cpu": sum(k["truth"]["cpu"] for k in kids),
                           "short": sum(max(0.0, cutoff - k["score"]["lb"]) for k in kids)}
            tag = inst + "_" + rec["node"].replace(":", "-").replace(",", "_")
            nodes.append((inst + " " + rec["node"], os.path.join(HERE, "dumps", tag.replace("_", "-", 0) + ".err"),
                          T, rec.get("rule9")))
    # tai60b 43:58: the four lines of ../tai60b_lines, by the rough subtree model there (README.md)
    T = {"row1": {"nodes": 19600}, "row45": {"nodes": 22900}, "col26": {"nodes": 38800}, "row8": {"nodes": 67100}}
    nodes.append(("tai60b 44:45,43:58 (arms' state)", os.path.join(HERE, "dumps", "tai60b_43-58.err"), T, "col26"))
    return nodes


def dump_name(node_name):
    inst, node = node_name.split()[:2]
    if inst == "tai60b":
        return os.path.join(HERE, "dumps", "tai60b_43-58.err")
    return os.path.join(HERE, "dumps", inst + "_" + node.replace(":", "-").replace(",", "-") + ".err")


SIGS = ["marg", "r2", "usum", "la", "zero", "short"]


def main():
    agg = {s: [] for s in SIGS}
    hits = {s: [0, 0, 0] for s in SIGS}   # truth-best line is the signal's top 1 / within top 2 / within top 3
    loss = {s: [] for s in SIGS}          # truth nodes of the signal's top line / the best line's
    count = 0
    for name, _, T, r9 in truth_tables():
        path = dump_name(name)
        if not os.path.exists(path):
            print(f"## {name}: no dump yet\n")
            continue
        d = parse_dump(path)
        if "LArow" not in d:
            print(f"## {name}: dump incomplete\n")
            continue
        S = line_signals(d)
        lines = [l for l in T if l in S]
        missing = [l for l in T if l not in S]
        tn = [T[l]["nodes"] for l in lines]
        print(f"## {name}  (bound {d['bound']:,.1f}, cutoff {d['inc'] - d['delta']:,.2f}, m {d['m']}; "
              f"rule 9 here: {d.get('rule9')}, in the study: {r9})\n")
        if missing:
            print("lines not free here:", missing)
        rank_t = ranks(tn, larger_better=False)
        cols = [s for s in SIGS if s != "short" or "short" in T[lines[0]]]
        vals = {}
        for s in cols:
            v = [(-T[l]["short"] if s == "short" else S[l].get(s, float("nan"))) for l in lines]
            vals[s] = v
        print("| line | truth nodes | rank | " + " | ".join(f"{s} rank" for s in cols) + " | rule-2 count |")
        print("|---|---|---|" + "---|" * len(cols) + "---|")
        rk = {s: ranks(vals[s]) for s in cols}
        for i in sorted(range(len(lines)), key=lambda i: tn[i]):
            print(f"| {lines[i]} | {tn[i]:,} | {rank_t[i]:g} | " + " | ".join(f"{rk[s][i]:g}" for s in cols) +
                  f" | {S[lines[i]]['r2_count']} |")
        best = min(tn)
        print("\nSpearman with the truth: " + ", ".join(f"{s} {spearman(vals[s], [-x for x in tn]):+.2f}" for s in cols))
        print()
        count += 1
        for s in cols:
            agg[s].append(spearman(vals[s], [-x for x in tn]))
            top = sorted(range(len(lines)), key=lambda i: rk[s][i])
            ib = tn.index(best)
            hits[s][0] += rk[s][ib] <= 1
            hits[s][1] += rk[s][ib] <= 2
            hits[s][2] += rk[s][ib] <= 3
            loss[s].append(tn[top[0]] / best)
    print(f"## Over the {count} nodes\n")
    print("| signal | mean Spearman | best line ranked 1st | in top 2 | in top 3 | its top line / best, mean (worst) |")
    print("|---|---|---|---|---|---|")
    for s in SIGS:
        if not agg[s]:
            continue
        n = len(agg[s])
        print(f"| {s} | {sum(agg[s]) / n:+.2f} | {hits[s][0]}/{n} | {hits[s][1]}/{n} | {hits[s][2]}/{n} | "
              f"{sum(loss[s]) / n:.2f}x ({max(loss[s]):.2f}x) |")


def arbiter():
    """rule 9's shortlist widened to the top K rows and K columns by marginal score, the line chosen by the
    children's first-solve shortfall (a): its truth nodes against the best of the study's lines"""
    print("\n## A DNN arbiter over the top K rows and K columns by marginal score\n")
    print("Truth nodes of the line chosen, over the best of the study's lines; (a) = the children's cold "
          "first solve.\n")
    print("| node | rule 9 (QPB arbiter, K = 1) | (a), K = 1 | K = 2 | K = 3 | K = 4 | all the study's lines |")
    print("|---|---|---|---|---|---|---|")
    tot = collections.defaultdict(list)
    for inst, pat in (("tho30", "tho30_solve/line_truth/results/*.json"), ("tho40", "tho40_line_truth/results/*.json")):
        for path in sorted(glob.glob(os.path.join(EXP, pat))):
            rec = json.load(open(path))
            cut = rec["ub"] - 1.95
            T = {}
            for line, L in rec["lines"].items():
                kids = L["children"].values()
                T[line] = (1 + sum(k["truth"]["nodes"] for k in kids),
                           sum(max(0.0, cut - k["score"]["lb"]) for k in kids), L["score"])
            best = min(v[0] for v in T.values())
            cells = [T[rec["rule9"]][0] / best]
            for K in (1, 2, 3, 4):
                cand = [l for side in ("row", "col") for l in
                        sorted((l for l in T if l.startswith(side)), key=lambda l: -T[l][2])[:K]]
                cells.append(T[min(cand, key=lambda l: T[l][1])][0] / best)
            cells.append(T[min(T, key=lambda l: T[l][1])][0] / best)
            for i, c in enumerate(cells):
                tot[i].append(c)
            print(f"| {inst} {rec['node']} | " + " | ".join(f"{c:.2f}x" for c in cells) + " |")
    print("| mean | " + " | ".join(f"{sum(v) / len(v):.2f}x" for _, v in sorted(tot.items())) + " |")


if __name__ == "__main__":
    main()
    arbiter()
