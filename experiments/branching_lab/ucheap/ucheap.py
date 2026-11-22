#!/usr/bin/env python3
"""ucheap: can the size score be read off QPB's dual matrix U, for free? (README.md)

  python3 ucheap.py [COGNAC]     COGNAC: a checkout of cognac with experiments/line_signals (default /home/user/wt-port)

A child (i, j) of a node has QPB's bound + U_ij as a lower bound, so U gives every child a shortfall
s_ij = max(0, cutoff - bound - U_ij) without a solve. The size score of DNN strong branching (sb.power) on those
shortfalls: a line's score is the sum over its children of (s_ij / s_max)^P, the least wins. P = 1 is rule 1's U
sum upside down; rule 2 is the count of children with s_ij = 0, ties by that sum.

A. The line-truth nodes (tho30 depth 1 x 3, tho40 depth 3 x 3: 8 lines each, every child's subtree solved to the
   end; tai60b 43:58: 4 lines by the line study's rough model): each score against the truth.
B. The grids' roots: U's shortfalls against the children's 150-iteration DNN shortfalls (../root, ../lsla/grids),
   cell by cell and row by row, and against the whole trees run from six nug30 rows and three tho30 rows.
"""
import glob, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
EXP = os.path.dirname(LAB)
COGNAC = sys.argv[1] if len(sys.argv) > 1 else "/home/user/wt-port"
sys.path.insert(0, os.path.join(COGNAC, "experiments", "line_signals"))
import signals as S   # parse_dump, ranks, spearman

POWERS = (1, 2, 4, 8)


def score(shortfalls, p, smax):
    return sum((s / smax) ** p for s in shortfalls if s > 0) if smax > 0 else 0.0


def pearson(a, b):
    n = len(a); ma, mb = sum(a) / n, sum(b) / n
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / den if den else float("nan")


def part_a():
    print("## A. The line-truth nodes\n")
    names = ["U, P = %d" % p for p in POWERS] + ["DNN first solve, P = %d" % p for p in POWERS]
    agg = {k: [] for k in names}; loss = {k: [] for k in names}; first = {k: 0 for k in names}
    cell_u, cell_d = [], []
    for inst, pat in (("tho30", "tho30_solve/line_truth/results/*.json"), ("tho40", "tho40_line_truth/results/*.json")):
        for path in sorted(glob.glob(os.path.join(EXP, pat))):
            rec = json.load(open(path))
            d = S.parse_dump(os.path.join(COGNAC, "experiments/line_signals/dumps",
                                          inst + "_" + rec["node"].replace(":", "-").replace(",", "-") + ".err"))
            cut = d["inc"] - d["delta"]
            ri = {r: k for k, r in enumerate(d["rows"])}; ci = {c: k for k, c in enumerate(d["cols"])}
            lines, truth, su, sd = [], [], {}, {}
            smax_u = smax_d = 0.0
            for line, L in rec["lines"].items():
                kids = L["children"]
                if not all("truth" in k and "nodes" in k["truth"] for k in kids.values()):
                    continue
                u, dn = [], []
                for cell, k in kids.items():
                    r, c = (int(x) for x in cell.split(":"))
                    u.append(max(0.0, cut - d["bound"] - d["U"][ri[r]][ci[c]]))
                    dn.append(0.0 if k["score"].get("fathomed") else max(0.0, cut - k["score"]["lb"]))
                    cell_u.append((u[-1], dn[-1], k["truth"]["nodes"]))
                lines.append(line); truth.append(1 + sum(k["truth"]["nodes"] for k in kids.values()))
                su[line], sd[line] = u, dn
                smax_u = max(smax_u, max(u)); smax_d = max(smax_d, max(dn))
            best = min(truth)
            row = []
            for nm in names:
                p = int(nm.split("= ")[1]); src, smax = (su, smax_u) if nm.startswith("U") else (sd, smax_d)
                v = [-score(src[l], p, smax) for l in lines]
                rho = S.spearman(v, [-t for t in truth])
                top = max(range(len(lines)), key=lambda i: v[i])
                agg[nm].append(rho); loss[nm].append(truth[top] / best); first[nm] += truth[top] == best
                row.append(f"{rho:+.2f} ({truth[top] / best:.2f}x)")
            fath = sum(1 for l in lines for s in su[l] if s == 0)
            print(f"- {inst} {rec['node']}: {len(lines)} lines, QPB {(cut - d['bound']) / cut:.1%} under the cutoff, "
                  f"children U fathoms: {fath}. Spearman with the truth (its first line / the best): "
                  + "; ".join(f"{n}: {r}" for n, r in zip(names, row)))
    print("\n| score | mean Spearman with the truth | best line first | its first line / the best, mean (worst) |\n|---|---|---|---|")
    for nm in names:
        n = len(agg[nm])
        print(f"| {nm} | {sum(agg[nm]) / n:+.2f} | {first[nm]}/{n} | {sum(loss[nm]) / n:.2f}x ({max(loss[nm]):.2f}x) |")
    open_cells = [(u, dn, t) for u, dn, t in cell_u]
    print(f"\nChild by child ({len(open_cells)} children of the studied lines): Spearman of the child's subtree nodes with "
          f"U's shortfall {S.spearman([c[0] for c in open_cells], [c[2] for c in open_cells]):+.2f}, with the DNN first "
          f"solve's {S.spearman([c[1] for c in open_cells], [c[2] for c in open_cells]):+.2f}; U's with the DNN's "
          f"{S.spearman([c[0] for c in open_cells], [c[1] for c in open_cells]):+.2f}.")


ROOTS = [   # name, DNN children of every row (150 iterations), cutoff of that run, whole trees (solver's nodes), production's row
    ("nug30", "root/runs/nug30_root150/children.txt", 6126 - 1.95,
     {23: 3623, 0: 4325, 1: 4566, 4: 4616, 6: 4854, 2: 10165}, 0),
    ("tho30", "root/runs/tho30_root150/children.txt", 149938 - 1.95, {4: 13812, 5: 14815, 25: 21095}, 5),
    ("tho40", "lsla/grids/runs/tho40/children.txt", 240518 - 1.95, {}, 26),
    ("sko49", "lsla/grids/runs/sko49/children.txt", 23386 - 1.95, {}, None),
    ("wil50", "lsla/grids/runs/wil50/children.txt", 48816 - 1.95, {}, None),
]


def part_b():
    print("\n## B. The roots\n")
    for name, kids_file, cut, trees, prod in ROOTS:
        d = S.parse_dump(os.path.join(HERE, "dumps", name + "_root.err"))
        dnn = {}
        for ln in open(os.path.join(LAB, kids_file)):
            f = ln.split()
            if f[0] != "child": continue
            r, c = (int(x) for x in f[f.index("cell") + 1].split(":"))
            dnn.setdefault(r, {})[c] = max(0.0, cut - float(f[f.index("lb") + 1]))
        rows = sorted(dnn)
        u = {r: {c: max(0.0, cut - d["bound"] - d["U"][r][c]) for c in dnn[r]} for r in rows}
        cu = [u[r][c] for r in rows for c in dnn[r]]; cd = [dnn[r][c] for r in rows for c in dnn[r]]
        smax_u, smax_d = max(cu), max(cd)
        print(f"### {name}: {len(rows)} rows x {len(dnn[rows[0]])} children; QPB {(cut - d['bound']) / cut:.1%} under the cutoff, "
              f"the DNN children {sum(cd) / len(cd) / cut:.2%} on average; children U fathoms: {sum(1 for x in cu if x == 0)}\n")
        print(f"- Cell by cell: U's shortfall against the DNN child's: Pearson {pearson(cu, cd):+.2f}, Spearman {S.spearman(cu, cd):+.2f}. "
              f"U's shortfalls span {min(cu) / smax_u:.2f} to 1 of their largest, the DNN's {min(cd) / smax_d:.2f} to 1.")
        ref = [-score(dnn[r].values(), 8, smax_d) for r in rows]        # the port's root score
        ref_rank = S.ranks(ref)
        out = []
        for p in POWERS:
            v = [-score(u[r].values(), p, smax_u) for r in rows]
            rk = S.ranks(v)
            top = rows[min(range(len(rows)), key=lambda i: rk[i])]
            s = f"P = {p}: first row r{top} (the DNN score's rank {ref_rank[rows.index(top)]:g} of {len(rows)}), Spearman with the DNN score {S.spearman(v, ref):+.2f}"
            if trees:
                tr = [r for r in rows if r in trees]
                s += f", with the trees {S.spearman([v[rows.index(r)] for r in tr], [-trees[r] for r in tr]):+.2f}"
                if top in trees: s += f", its tree {trees[top]:,}"
            out.append(s)
        print("- U's row score: " + "; ".join(out) + ".")
        dtop = rows[min(range(len(rows)), key=lambda i: ref_rank[i])]
        v8 = S.ranks([-score(u[r].values(), 8, smax_u) for r in rows])
        s = f"- The DNN score's first row r{dtop}: U's rank {v8[rows.index(dtop)]:g} (P = 8)"
        if prod is not None:
            s += f"; production's r{prod}: DNN rank {ref_rank[rows.index(prod)]:g}, U's rank {v8[rows.index(prod)]:g}"
        print(s + ".")
        if trees:
            tr = sorted(trees, key=trees.get)
            print("- Rows with a whole tree, best first: " + ", ".join(
                f"r{r} {trees[r]:,} (DNN rank {ref_rank[rows.index(r)]:g}, U rank {v8[rows.index(r)]:g})" for r in tr)
                + f". Spearman of the DNN score with the trees {S.spearman([ref[rows.index(r)] for r in tr], [-trees[r] for r in tr]):+.2f}.")
        print()


if __name__ == "__main__":
    part_a()
    part_b()
