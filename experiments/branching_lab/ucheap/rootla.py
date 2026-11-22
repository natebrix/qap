#!/usr/bin/env python3
"""ucheap/rootla.py: the size score on the children's own QPB bounds, at the grids' roots (README.md).

  python3 rootla.py

dumps/*_rootla.err: every row's children bounded by QPB as rule 4's row lookahead bounds them (ls-LAcell row col
bound; the experiment build with QAP_LINE_SCORES_LA at a root with symmetry left). A row's score from those
shortfalls s = cutoff - bound: the sum (P = 1), sum of (s / s_max)^P, sum of exp(L s / G) with G the root's QPB gap,
or the largest. Against the DNN strong branching's root score (P = 8 on the children's 150-iteration DNN bounds)
and the whole trees.
"""
import math, os, sys
sys.path.insert(0, "/home/user/wt-port/experiments/line_signals")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import signals as S, ucheap as UC

for name, kids_file, cut, trees, prod in UC.ROOTS:
    q, choice, bound = {}, None, None
    for ln in open(os.path.join(HERE, "dumps", name + "_rootla.err")):
        f = ln.split()
        if f[0] == "ls-LAcell": q.setdefault(int(f[1]), {})[int(f[2])] = cut - float(f[3])
        elif f[0] == "ls-LAchoice:": choice = int(f[2])
        elif f[0] == "ls-node:": bound = float(f[f.index("bound") + 1])
    dnn = {}
    for ln in open(os.path.join(UC.LAB, kids_file)):
        f = ln.split()
        if f[0] != "child": continue
        r, c = (int(x) for x in f[f.index("cell") + 1].split(":"))
        dnn.setdefault(r, {})[c] = max(0.0, cut - float(f[f.index("lb") + 1]))
    rows = sorted(dnn)
    assert all(set(q[r]) == set(dnn[r]) for r in rows), name
    G = cut - bound
    allq = [v for r in rows for v in q[r].values()]; alld = [v for r in rows for v in dnn[r].values()]
    smax, smax_d = max(allq), max(alld)
    ref = [-UC.score(dnn[r].values(), 8, smax_d) for r in rows]; rr = S.ranks(ref)
    print(f"### {name}: the QPB row lookahead's own choice r{choice}" + (f" (production's r{prod})" if prod is not None else "")
          + f"; children's QPB shortfalls {min(allq) / G:.2f} to {smax / G:.2f} of the root's QPB gap; cell by cell with the "
          f"DNN children's: Pearson {UC.pearson(allq, alld):+.2f}, Spearman {S.spearman(allq, alld):+.2f}\n")
    print("| score on the children's QPB bounds | first row | its rank by the DNN score | the DNN score's first row here | Spearman with the DNN score | with the trees | first row's tree |")
    print("|---|---|---|---|---|---|---|")
    forms = [("sum (P = 1)", lambda r: sum(q[r].values()))] + \
            [(f"power {p}", lambda r, p=p: sum((s / smax) ** p for s in q[r].values())) for p in (8, 32)] + \
            [(f"exp, L = {L}", lambda r, L=L: sum(math.exp(L * (s - smax) / G) for s in q[r].values())) for L in (16, 32, 64)] + \
            [("largest shortfall", lambda r: max(q[r].values()))]
    dtop = rows[min(range(len(rows)), key=lambda i: rr[i])]
    for nm, fn in forms:
        v = [-fn(r) for r in rows]; rk = S.ranks(v)
        top = rows[min(range(len(rows)), key=lambda i: rk[i])]
        tr = [r for r in rows if r in trees]
        print(f"| {nm} | r{top} | {rr[rows.index(top)]:g} of {len(rows)} | rank {rk[rows.index(dtop)]:g} | {S.spearman(v, ref):+.2f} | "
              + (f"{S.spearman([v[rows.index(r)] for r in tr], [-trees[r] for r in tr]):+.2f}" if tr else "") + " | "
              + (f"{trees[top]:,}" if top in trees else "") + " |")
    if trees:
        v = [-forms[4][1](r) for r in rows]; rk = S.ranks(v)
        print("\nRows with a whole tree, best first: " + ", ".join(
            f"r{r} {trees[r]:,} (exp L = 32 rank {rk[rows.index(r)]:g}, DNN rank {rr[rows.index(r)]:g})" for r in sorted(trees, key=trees.get)))
    print()
