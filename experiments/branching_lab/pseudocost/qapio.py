"""Shared readers: QAPLIB instances and solutions, the solver's node traces (QAP_NODE_TRACE) and DNN traces."""
import gzip, os
import numpy as np

COG = "/home/user/wt-cognac"


def read_dat(path):
    tok = open(path).read().split()
    n = int(tok[0]); v = np.array(tok[1:1 + 2 * n * n], dtype=float)
    return n, v[:n * n].reshape(n, n), v[n * n:].reshape(n, n)


def read_sln(path):
    tok = open(path).read().split()
    n = int(tok[0]); obj = int(tok[1]); p = [int(x) - 1 for x in tok[2:2 + n]]
    return obj, p


def obj(A, B, p):
    p = np.asarray(p)
    return float((A * B[np.ix_(p, p)]).sum())


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_node_trace(path):
    """list of dicts in trace order; fixes as a list of (row, col) in history order"""
    out = []
    num = lambda s: float("nan") if s == "-" else float(s)
    for k, ln in enumerate(open_any(path)):
        f = ln.split()
        if len(f) < 12 or f[0] != "node":
            continue
        fixes = [] if f[2] == "-" else [tuple(int(x) for x in p.split(":")) for p in f[2].split(",")]
        out.append(dict(idx=k, depth=int(f[1]), fixes=fixes, qpb=num(f[3]), dnn=num(f[4]), cut=num(f[5]),
                        cutoff=num(f[6]), outcome=f[7], rule=f[8], line=f[9], kids=f[10], cpu=float(f[11])))
    return out
