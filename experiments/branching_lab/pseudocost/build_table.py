#!/usr/bin/env python3
"""Step 1: one row per (parent, child) from the node traces (QAP_NODE_TRACE) and the older DNN traces
(QAP_DNN_TRACE, tho30 and tho40 of Nov 8), with the instance's structural features of the child's fix.

  python3 build_table.py            # writes data/<run>.pkl (not committed: see .gitignore)

Conventions, checked against QAPSOLN (check_orientation() below, output in build_table_out.txt):
a fix r:c assigns row r of the instance's first matrix A to row/column c of the second, B,
objective sum_ij A[i,j] B[p(i),p(j)]. For tai50b/tai60b/nug30 the .sln permutation p agrees with the
traces' fixes (the best known solution's chain is in the trace); for tho40t.dat it is the inverse of
QAPSOLN/tho40.sln (tho40t is the transposed copy). tai*b: A is the symmetric clustered distance-like
matrix, B the asymmetric flows. tho/nug: A the flows, B the grid distances.

Targets and censoring (README, "Censoring"):
  fath1   1 when the first DNN solve fathomed the child (outcome D) or QPB did (Q); 0 when it did not (C, B)
  splus   max(0, cutoff - first-solve DNN bound): exactly 0 for D and Q (the solve stops when it passes
          the cutoff, so the recorded bound is a floor, but its shortfall's positive part is exact); for
          C and B, the shortfall where the solve stopped (iteration limit or the give-up's extrapolation)
  rel     splus / (cutoff - parent's bound): the share of the parent's gap the child leaves open
  gaincap min(dnn, cutoff) - parent bound = gap - splus: the gain, capped at the cutoff (uncensored for D)
P (discarded when popped) and E (enumerated) children are dropped, and children whose first solve
did not run (no DNN bound and not Q).
"""
import collections, os, pickle, sys
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from qapio import COG, read_dat, read_sln, obj, read_node_trace, open_any

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

RUNS = {  # run: (instance file, sln, inverse?)
    "tai60b-261117": ("tai60b", "tai60b", False),
    "tai50b-261116": ("tai50b", "tai50b", False),
    "tai50b-r9": ("tai50b", "tai50b", False),
    "tho40ab-32": ("tho40t", "tho40", True),
    "tho40ab-0": ("tho40t", "tho40", True),
    "nug30-trace-a": ("nug30", "nug30", False),
    "nug30-c7a": ("nug30", "nug30", False),
}
DNN_TRACES = {  # older QAP_DNN_TRACE files: no QPB, no line
    "tho30-dnn": ("tho30", None, False, COG + "/experiments/tho40_estimate/traces/tho30.trace.gz"),
    "tho40-dnn": ("tho40t", "tho40", True, COG + "/experiments/tho40_estimate/traces/tho40.trace.gz"),
}


class Inst:
    def __init__(self, name, sln, inverse):
        self.n, self.A, self.B = read_dat(f"{COG}/qaplib/{name}.dat")
        n, A, B = self.n, self.A, self.B
        self.best = None
        if sln:
            o, p = read_sln(f"{COG}/QAPSOLN/{sln}.sln")
            if len(p) == n:
                if inverse:
                    q = [0] * n
                    for i, j in enumerate(p): q[j] = i
                    p = q
                assert abs(obj(A, B, p) - o) < 0.5, (name, obj(A, B, p), o)
                self.best = p
        As, Bs = A + A.T, B + B.T
        self.rsum = As.sum(1); self.rnnz = (As > 0).sum(1) / n
        self.csum = Bs.sum(1); self.cnnz = (Bs > 0).sum(1) / n
        # clusters on the distance-like side (dense, symmetric, zero diagonal): A for tai, B for nug/tho
        self.rclu = self._clusters(A); self.cclu = self._clusters(B)

    @staticmethod
    def _clusters(M):
        n = len(M)
        if not (M == M.T).all() or (M > 0).sum() < 0.9 * n * (n - 1):
            return np.zeros(n, dtype=int)
        Z = linkage(squareform(M, checks=False), "average")
        h = Z[:, 2]; gaps = np.diff(h)
        # cut at the largest jump among the last 12 merges (k = 2..13 clusters), if it is a clear one
        k0 = max(0, len(h) - 13)
        j = k0 + int(np.argmax(gaps[k0:]))
        if gaps[j] < 0.25 * h[-1]:
            return np.zeros(n, dtype=int)
        return fcluster(Z, t=(h[j] + h[j + 1]) / 2, criterion="distance") - 1

    def fix_features(self, fixes, r, c):
        """features of adding r:c to the fixes (list of (row, col)) of the parent"""
        A, B, n = self.A, self.B, self.n
        fr = [a for a, _ in fixes]; fc = [l for _, l in fixes]
        dfix = float((A[r, fr] * B[c, fc]).sum() + (A[fr, r] * B[fc, c]).sum()) if fixes else 0.0
        used_r = set(fr) | {r}; used_c = set(fc) | {c}
        R = [i for i in range(n) if i not in used_r]; C = [j for j in range(n) if j not in used_c]
        a = np.sort(A[r, R] + A[R, r]); b = np.sort(B[c, C] + B[C, c])[::-1]
        gl = float((a * b).sum()) / 2   # Gilmore-Lawler term of the pair on the free part (symmetrised)
        return dfix, gl


def check_orientation(out):
    for run, (name, sln, inv) in RUNS.items():
        I = Inst(name, sln, inv)
        tr = read_node_trace(f"{COG}/infra/phase1/runs/{run}/nodes.txt.gz")
        if I.best is None: continue
        q = I.best; qi = [0] * I.n
        for i, j in enumerate(q): qi[j] = i
        dp = max((t["depth"] for t in tr if all(q[a] == l for a, l in t["fixes"])), default=0)
        di = max((t["depth"] for t in tr if all(qi[a] == l for a, l in t["fixes"])), default=0)
        print(f"{run}: deepest node agreeing with the best known solution as row->col: {dp}; "
              f"as its inverse: {di} (n {I.n})", file=out)


def child_rows(run, I, nodes, has_qpb=True):
    """nodes: list of dicts (idx, depth, fixes in history order, qpb, dnn, cut, cutoff, outcome, rule, line)"""
    by = {}
    for t in nodes:   # keep the last decided copy of a node (a killed unit's node is traced again)
        k = tuple(t["fixes"])
        if k in by and t["outcome"] == "P" and by[k]["outcome"] != "P": continue
        by[k] = t
    rows = []
    nan = float("nan")
    for k, t in by.items():
        if not k: continue
        par = by.get(k[:-1])
        if par is None or par["outcome"] != "B": continue
        o = t["outcome"]
        if o in ("P", "E"): continue
        if np.isnan(t["dnn"]) and o != "Q": continue
        cutoff = t["cutoff"] if not np.isnan(t["cutoff"]) else par["cutoff"]
        pb = np.nanmax([par["qpb"], par["dnn"], par["cut"]])
        gap = cutoff - pb
        if not gap > 0: continue
        r, c = k[-1]
        sp = 0.0 if o == "Q" else max(0.0, cutoff - t["dnn"])
        if o == "D": sp = 0.0
        line = par["line"]
        dfix, gl = I.fix_features(list(k[:-1]), r, c)
        rows.append(dict(
            run=run, pkey=",".join(f"{a}:{b}" for a, b in k[:-1]), pidx=par["idx"], cidx=t["idx"],
            depth=t["depth"], prule=par["rule"], pline=line,
            line_is_row=(line[0] == "r") if isinstance(line, str) and line not in ("-", "") else np.nan,
            row=r, col=c, cutoff=cutoff, pbound=pb, pgap=gap, prelgap=gap / abs(cutoff),
            pqpb=par["qpb"], pdnn=par["dnn"], pcut=par["cut"],
            pshort_rel=(cutoff - par["dnn"]) / abs(cutoff) if not np.isnan(par["dnn"]) else nan,
            pcut_ran=float(not np.isnan(par["cut"])),
            pdnn_minus_qpb=(par["dnn"] - par["qpb"]) / abs(cutoff) if has_qpb and not np.isnan(par["dnn"]) else nan,
            cqpb=t["qpb"], cqpb_short=(cutoff - t["qpb"]) / gap if has_qpb and not np.isnan(t["qpb"]) else nan,
            cdnn=t["dnn"], ccut=t["cut"], outcome=o, cstatus=t.get("status", ""), ccpu=t["cpu"],
            fath1=float(o in ("D", "Q")), fath_any=float(o in ("D", "Q", "C")),
            splus=sp, rel=sp / gap, gaincap=gap - sp,
            agree=float(I.best[r] == c) if I.best is not None else nan,
            par_agree=float(all(I.best[a] == b for a, b in k[:-1])) if I.best is not None else nan,
            rsum=I.rsum[r], rnnz=I.rnnz[r], csum=I.csum[c], cnnz=I.cnnz[c],
            rclu=int(I.rclu[r]), cclu=int(I.cclu[c]),
            dfix=dfix / abs(cutoff), gl=gl / abs(cutoff), nfree=I.n - t["depth"],
        ))
    return rows


def read_dnn_trace(path):
    """QAP_DNN_TRACE lines as node dicts; fixes sorted by facility, so the history is rebuilt from the
    parent that exists (the set less one pair); no QPB, no line: a parent is any node with children."""
    nodes = []
    for k, ln in enumerate(open_any(path)):
        f = ln.split("|")
        h = f[0].split()
        depth = int(h[3]); fx = [tuple(int(x) for x in p.split(":")) for p in h[5:]]
        first = f[1].split()[1:]; kv = dict(zip(first[::2], first[1::2]))
        cut_lb = float("nan")
        for seg in f[2:]:
            s = seg.split()
            if s and s[0] == "cut" and len(s) > 1 and s[1] == "1":
                cut_lb = float(s[s.index("best") + 1])
            if s and s[0] == "cutoff":
                cutoff = float(s[1]); fath = int(s[3])
        nodes.append(dict(idx=k, depth=depth, fixset=frozenset(fx), dnn=float(kv["lb"]), status=kv["status"],
                          cut=cut_lb, cutoff=cutoff, fathomed=fath, qpb=float("nan")))
    byset = {t["fixset"]: t for t in nodes}
    if frozenset() not in byset:   # the root is not DNN-bounded in these traces: a bare parent for depth 1
        root = dict(idx=-1, depth=0, fixset=frozenset(), dnn=float("nan"), status="none", cut=float("nan"),
                    cutoff=nodes[0]["cutoff"], fathomed=0, qpb=float("nan"))
        byset[frozenset()] = root; nodes = [root] + nodes
    # history order: walk up through existing subsets
    out = []
    for t in nodes:
        fs = t["fixset"]
        cands = [p for p in fs if (fs - {p}) in byset]
        if len(fs) > 0 and len(cands) != 1:
            continue   # root-level, orphan or ambiguous
        t["last"] = cands[0] if fs else None
        out.append(t)
    have_kids = collections.Counter()
    for t in out:
        if t["last"] is not None: have_kids[t["fixset"] - {t["last"]}] += 1
    hist = {}
    def history(fs):
        if fs in hist: return hist[fs]
        t = byset[fs]
        if not fs: h = []
        elif t.get("last") is None: return None
        else:
            hp = history(fs - {t["last"]})
            h = None if hp is None else hp + [t["last"]]
        hist[fs] = h
        return h
    res = []
    for t in out:
        h = history(t["fixset"]) if t["fixset"] in byset else None
        if h is None and t["fixset"]: continue
        fathomed_first = t["status"] == "fathomed"
        o = "D" if fathomed_first else ("C" if t["fathomed"] else ("B" if have_kids[t["fixset"]] else "X"))
        res.append(dict(idx=t["idx"], depth=t["depth"], fixes=h or [], qpb=float("nan"), dnn=t["dnn"],
                        cut=t["cut"], cutoff=t["cutoff"], outcome=o, rule="-", line="-", cpu=float("nan"),
                        status=t["status"]))
    return res


def main():
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(HERE, "build_table_out.txt"), "w") as out:
        check_orientation(out)
        for run, (name, sln, inv) in RUNS.items():
            I = Inst(name, sln, inv)
            nodes = read_node_trace(f"{COG}/infra/phase1/runs/{run}/nodes.txt.gz")
            rows = child_rows(run, I, nodes)
            pickle.dump(rows, open(os.path.join(DATA, run + ".pkl"), "wb"))
            oc = collections.Counter(r["outcome"] for r in rows)
            print(f"{run}: {len(nodes)} trace lines, {len(rows)} child rows, outcomes {dict(oc)}, "
                  f"row clusters {len(set(I.rclu))}, col clusters {len(set(I.cclu))}", file=out)
        for run, (name, sln, inv, path) in DNN_TRACES.items():
            I = Inst(name, sln, inv)
            nodes = read_dnn_trace(path)
            # an "X" node: unfathomed, no children in the trace (its subtree below the DNN depths, or unfinished)
            rows = child_rows(run, I, nodes, has_qpb=False)
            pickle.dump(rows, open(os.path.join(DATA, run + ".pkl"), "wb"))
            oc = collections.Counter(r["outcome"] for r in rows)
            print(f"{run}: {len(nodes)} DNN-trace nodes placed, {len(rows)} child rows, outcomes {dict(oc)}", file=out)
    print(open(os.path.join(HERE, "build_table_out.txt")).read())


if __name__ == "__main__":
    main()
