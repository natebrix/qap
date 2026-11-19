#!/usr/bin/env python3
"""Step 3: does a model trained on an instance's traces rank the candidate lines at the line-truth nodes?

  OMP_NUM_THREADS=2 python3 truth.py > truth_out.md

The truth (read-only, /home/user/wt-cognac/experiments): tho30_solve/line_truth (tho30, depth-1 nodes 5:1,
5:2, 5:3) and tho40_line_truth (tho40t, three depth-3 nodes), 8 lines each, every child's cold first DNN
solve (a) and its whole subtree (b); and tai60b's 44:45,43:58 (four lines, the rough subtree model of
tai60b_lines). line_signals/dumps give each node's U (QPB's dual), the DNN marginals X and rule 4's
lookahead score per line.

Training sets, never containing the truth node's subtree (the node's own children included):
  tho30   tho30-dnn (the Nov 8 DNN trace of tho30, rule 9) less the subtree of each truth node
  tho40   tho40ab-32 (the brief's choice: another depth-2 subtree, 26:10,8:0); also tho40-dnn (the Nov 8
          trace the truth nodes were taken from) less their subtrees, and both together
  tai60b  tai60b-261117 less the subtree of 44:45,43:58
A child is a subproblem whatever line created it (fixes F + (r, c)), so its features do not name the line.
The model: gbt (common.m_gbt, no child QPB: the truth sets do not have the children's QPB bounds) on ylog;
a line's score is the sum over its children of the predicted shortfall (ylog transformed back). Also
pc (the pseudocost table) and the expected count of open children (a classifier would do; here the
share of predictions above the median of an open child's in training is not used: the sum is the score).
Child-level baselines at the truth nodes: bound + U_ij (QPB's predicted child bound) and the marginal x_ij."""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from common import *
from build_table import Inst, read_dnn_trace, DNN_TRACES

EXP = "/home/user/wt-cognac/experiments"
DUMPS = EXP + "/line_signals/dumps"


def parse_dump(path):
    d = {"U": [], "X": []}
    for ln in open(path):
        if ln.startswith("ls-node:"):
            f = ln.split(); kv = dict(zip(f[1::2], f[2::2]))
            d.update(depth=int(kv["depth"]), m=int(kv["m"]), bound=float(kv["bound"]), inc=float(kv["inc"]),
                     delta=float(kv["delta"]))
        elif ln.startswith("ls-rows:"): d["rows"] = [int(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-cols:"): d["cols"] = [int(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-U "): d["U"].append([float(x) for x in ln.split(":", 1)[1].split()])
        elif ln.startswith("ls-X "): d["X"].append([float(x) for x in ln.split(":", 1)[1].split()])
        elif ln.startswith("ls-LArow:"): d["LArow"] = [float(x) for x in ln.split()[1:]]
        elif ln.startswith("ls-LAcol:"): d["LAcol"] = [float(x) for x in ln.split()[1:]]
    return d


def ranks(v, larger_better=True):
    return list(rankdata([-x for x in v] if larger_better else v))


def sp(a, b):
    return spearmanr(a, b).statistic


def fixes_of(s):
    return [tuple(int(x) for x in p.split(":")) for p in s.split(",")]


def node_info(table, fx):
    """the truth node's own row in a child table (it is a child of its parent there)"""
    S = set(fx)
    for _, r in table[table.depth == len(fx)].iterrows():
        pk = set(fixes_of(r.pkey)) if r.pkey else set()
        if pk | {(r.row, r.col)} == S:
            return r
    return None


def in_subtree(table, fx):
    S = set(fx)
    return table.apply(lambda r: S <= (set(fixes_of(r.pkey)) if r.pkey else set()) | {(r.row, r.col)}, axis=1)


def truth_nodes():
    out = []
    for inst, pat, dat, sln, inv, train in (
            ("tho30", "tho30_solve/line_truth/results/*.json", "tho30", None, False, ["tho30-dnn"]),
            ("tho40", "tho40_line_truth/results/*.json", "tho40t", "tho40", True, ["tho40ab-32"])):
        for path in sorted(glob.glob(os.path.join(EXP, pat))):
            rec = json.load(open(path))
            cutoff = rec["ub"] - 1.95
            T = {}
            for line, L in rec["lines"].items():
                kids = L["children"]
                T[line] = dict(nodes=1 + sum(k["truth"]["nodes"] for k in kids.values()),
                               short=sum(max(0.0, cutoff - k["score"]["lb"]) for k in kids.values()),
                               kids={tuple(int(x) for x in f.split(":")): (max(0.0, cutoff - k["score"]["lb"]),
                                                                          k["score"]["fathomed"], k["truth"]["nodes"])
                                     for f, k in kids.items()})
            dump = os.path.join(DUMPS, f"{inst}_" + rec["node"].replace(":", "-").replace(",", "-") + ".err")
            out.append(dict(name=f"{inst} {rec['node']}", inst=dat, sln=sln, inv=inv, fixes=fixes_of(rec["node"]),
                            cutoff=cutoff, T=T, dump=dump, train=train, rule9=rec.get("rule9")))
    T = {l: dict(nodes=v, kids=None) for l, v in (("row1", 19600), ("row45", 22900), ("col26", 38800), ("row8", 67100))}
    out.append(dict(name="tai60b 44:45,43:58", inst="tai60b", sln="tai60b", inv=False, fixes=[(44, 45), (43, 58)],
                    cutoff=608215054.05, T=T, dump=os.path.join(DUMPS, "tai60b_43-58.err"), train=["tai60b-261117"],
                    rule9="col26"))
    return out


def candidate_rows(N, I, info, D):
    """one feature row per child of each candidate line"""
    fx = N["fixes"]; cutoff = N["cutoff"]; a = abs(cutoff)
    frow = {r for r, _ in fx}; fcol = {c for _, c in fx}
    free_r = [r for r in range(I.n) if r not in frow]; free_c = [c for c in range(I.n) if c not in fcol]
    pdnn = info.cdnn if info is not None else np.nan
    pb = np.nanmax([info.cqpb, info.cdnn, info.ccut]) if info is not None else D["bound"]
    rows = []
    for line in N["T"]:
        kind, k = line[:3], int(line[3:])
        pairs = [(k, c) for c in free_c] if kind == "row" else [(r, k) for r in free_r]
        for r, c in pairs:
            dfix, gl = I.fix_features(fx, r, c)
            u = x = np.nan
            if D and r in D["rows"] and c in D["cols"]:
                i, j = D["rows"].index(r), D["cols"].index(c)
                u = D["U"][i][j]; x = D["X"][i][j] if D["X"] else np.nan
            rows.append(dict(line=line, pkey=line, row=r, col=c, depth=len(fx) + 1, nfree=I.n - len(fx) - 1,
                             cutoff=cutoff, pshort=(cutoff - pdnn) / a, prelgap=(cutoff - pb) / a,
                             pcut_ran=float(info is not None and not np.isnan(info.ccut)),
                             rsum=I.rsum[r], rnnz=I.rnnz[r], csum=I.csum[c], cnnz=I.cnnz[c],
                             rclu=int(I.rclu[r]), cclu=int(I.cclu[c]), dfix=dfix / a, gl=gl / a,
                             agree=float(I.best[r] == c) if I.best is not None else np.nan,
                             par_agree=float(all(I.best[p] == q for p, q in fx)) if I.best is not None else np.nan,
                             U=u, X=x, Ubound=D["bound"] + u if D else np.nan))
    return pd.DataFrame(rows)


def main():
    out = sys.stdout
    print("# Step 3: ranking the candidate lines at the line-truth nodes (truth.py)\n", file=out)
    variants = {"tho40": [("tho40ab-32", ["tho40ab-32"]), ("tho40-dnn less the truth subtrees", ["tho40-dnn"]),
                          ("both", ["tho40ab-32", "tho40-dnn"])]}
    tables = {}
    agg = {}
    childagg = {}
    for N in truth_nodes():
        I = Inst(N["inst"], N["sln"], N["inv"])
        D = parse_dump(N["dump"]) if os.path.exists(N["dump"]) else None
        inst = N["name"].split()[0]
        print(f"## {N['name']}\n", file=out)
        for vname, runs in variants.get(inst, [(N["train"][0], N["train"])]):
            parts = []
            info = None
            for run in runs:
                if run not in tables: tables[run] = load(run)
                t = tables[run]
                if info is None: info = node_info(t, N["fixes"])
                parts.append(t[~in_subtree(t, N["fixes"])])
            tr = pd.concat(parts, ignore_index=True)
            if info is None:   # not in the training trace (tho40ab-32), or a depth-1 node: the Nov 8 DNN trace
                src = {"tho30": "tho30-dnn", "tho40": "tho40-dnn"}[inst]
                for t in read_dnn_trace(DNN_TRACES[src][3]):
                    if set(t["fixes"]) == set(N["fixes"]):
                        info = pd.Series(dict(cdnn=t["dnn"], ccut=t["cut"], cqpb=np.nan))
                        break
            C = candidate_rows(N, I, info, D)
            fs = [f for f in feats(tr, False) if f in C]
            g = m_gbt(tr, False, fs=fs)(C)
            p = m_pseudocost(tr)(C.assign(pair=C.row * 1000 + C.col))
            a = abs(N["cutoff"])
            C["pred_s"] = np.expm1(np.clip(g, 0, None)) / SCALE * a
            C["pred_pc"] = np.expm1(np.clip(p, 0, None)) / SCALE * a
            lines = list(N["T"])
            truth = [N["T"][l]["nodes"] for l in lines]
            score = {"gbt: predicted shortfall sum": [-C[C.line == l].pred_s.sum() for l in lines],
                     "pc: predicted shortfall sum": [-C[C.line == l].pred_pc.sum() for l in lines]}
            if D and "LArow" in D:
                la = {**{f"row{r}": v for r, v in zip(D["rows"], D["LArow"])},
                      **{f"col{c}": v for c, v in zip(D["cols"], D["LAcol"])}}
                score["QPB lookahead (rule 4)"] = [la[l] for l in lines]
            if "short" in N["T"][lines[0]]:
                score["(a) the children's first DNN solve"] = [-N["T"][l]["short"] for l in lines]
            print(f"training: {vname}, {len(tr):,} children; node first solve "
                  f"{fmt(info.cdnn if info is not None else np.nan, 1)}, cutoff {N['cutoff']:,.2f}\n", file=out)
            print("| line | truth nodes | rank | " + " | ".join(f"{s.split(':')[0]} rank" for s in score) +
                  " | gbt shortfall sum |", file=out)
            print("|---|---|---|" + "---|" * len(score) + "---|", file=out)
            rk = {s: ranks(v) for s, v in score.items()}
            rt = ranks(truth, larger_better=False)
            for i in sorted(range(len(lines)), key=lambda i: truth[i]):
                print(f"| {lines[i]} | {truth[i]:,} | {rt[i]:g} | " + " | ".join(f"{rk[s][i]:g}" for s in score) +
                      f" | {-score['gbt: predicted shortfall sum'][i]:,.0f} |", file=out)
            best = min(truth)
            print("\nSpearman with the truth: " + ", ".join(
                f"{s.split(':')[0]} {sp(v, [-x for x in truth]):+.2f}" for s, v in score.items()), file=out)
            for s, v in score.items():
                top = int(np.argmax(v))
                key = (inst if inst != "tho40" else f"tho40 [{vname}]", s)
                agg.setdefault(key, []).append((sp(v, [-x for x in truth]), truth[top] / best,
                                                rt[top] if False else None, ranks(v)[truth.index(best)]))
            # child level: the (a) solves of every candidate child
            if N["T"][lines[0]].get("kids"):
                C["a_short"] = [N["T"][l]["kids"][(r, c)][0] for l, r, c in zip(C.line, C.row, C.col)]
                C["a_open"] = [1 - int(N["T"][l]["kids"][(r, c)][1]) for l, r, c in zip(C.line, C.row, C.col)]
                ch = {"gbt": C.pred_s.values, "pc": C.pred_pc.values, "-(bound + U_ij)": -C.Ubound.values,
                      "x_ij (marginal; larger = more likely open)": C.X.values}
                cells = []
                for s, v in ch.items():
                    # within each line (the line's children), then over the node's children
                    spl = np.nanmean([sp(v[C.line.values == l], C.a_short.values[C.line.values == l]) for l in lines])
                    k = C.a_open.values.astype(bool)
                    r = rankdata(v); auc = (r[k].sum() - k.sum() * (k.sum() + 1) / 2) / (k.sum() * (~k).sum())
                    cells.append(f"{s} {spl:+.2f} / {auc:.2f}")
                    childagg.setdefault((inst if inst != "tho40" else f"tho40 [{vname}]", s), []).append((spl, auc))
                print(f"\nChild level, the node's {len(C)} candidate children against their (a) solves "
                      f"({int(C.a_open.sum())} open): Spearman within a line / AUC open vs fathomed over the node: "
                      + "; ".join(cells), file=out)
            print(file=out)
            out.flush()
    print("## Over the nodes\n", file=out)
    print("| training | score | nodes | mean Spearman | best line 1st | its top line / best, mean (worst) |", file=out)
    print("|---|---|---|---|---|---|", file=out)
    for (tr, s), v in agg.items():
        sps = [x[0] for x in v]; loss = [x[1] for x in v]; first = sum(1 for x in v if x[3] <= 1)
        print(f"| {tr} | {s} | {len(v)} | {np.nanmean(sps):+.2f} | {first}/{len(v)} | {np.mean(loss):.2f}x ({max(loss):.2f}x) |", file=out)
    print("\n## Child level over the nodes (Spearman within a line / AUC over the node)\n", file=out)
    print("| training | score | nodes | Spearman within a line | AUC open |", file=out)
    print("|---|---|---|---|---|", file=out)
    for (tr, s), v in childagg.items():
        print(f"| {tr} | {s} | {len(v)} | {np.nanmean([x[0] for x in v]):+.2f} | {np.mean([x[1] for x in v]):.2f} |", file=out)


if __name__ == "__main__":
    main()
