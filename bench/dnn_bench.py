#!/usr/bin/env python3
"""Benchmark the DNN bound (latest/dnntest) against the Python prototype.

  python3 bench/dnn_bench.py iter  [--out F]   # CPU per ADMM iteration, m = 18, 22, 25, 28
  python3 bench/dnn_bench.py nodes [--out F] [--count N] [--proto]
                                               # CPU per node on QPB-open depth-2 nodes

Everything runs on one core (BLAS threads = 1), one job at a time. Results
are appended to a JSONL file as they are produced and finished items are
skipped on a rerun, so an interrupted run can simply be restarted.

iter:  dnntest in double (dsyevd) and single precision (-f), and dnn.py,
       each on the same depth-2 node of nug20 / nug24 / nug27 / nug30.
nodes: depth-2 frontiers come from bnbdrv (PARAM/rule2.par, incumbent =
       optimum, -d 1 -o state: the children of the depth-1 nodes), as in
       the nug20/nug30 studies. A seeded random sample of them is bounded
       with the solver's cutoff (optimum - 1.95). With --proto each node is
       also run in prototype mode (dnntest -P), which reproduces dnn.py's
       iterates and stopping rule, so its iteration count times dnn.py's
       measured CPU per iteration estimates the prototype's CPU per node.
"""

import argparse
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LATEST = os.path.join(ROOT, "latest")
DNNTEST = os.path.join(LATEST, "dnntest")
PROTO = os.path.join(ROOT, "experiments", "dnn_bound", "nug20")

# instance, optimum, m of a depth-2 node, node used for the per-iteration timing
INSTANCES = [("nug20", 2570, 18, "0:19,1:0"), ("nug24", 3488, 22, "0:0,1:1"),
             ("nug27", 5234, 25, "0:0,1:1"), ("nug30", 6124, 28, "4:29,6:7")]


def dat(inst):
    return os.path.join(ROOT, "qaplib", inst + ".dat")


def load_done(out, key):
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                done.add(key(json.loads(line)))
            except (ValueError, KeyError):
                pass
    return done


def append(out, rec):
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


def dnntest(args):
    out = subprocess.run([DNNTEST] + args, capture_output=True, text=True)
    res = [dict(kv.split("=", 1) for kv in l.split()[1:])
           for l in out.stdout.splitlines() if l.startswith("RESULT")]
    if out.returncode != 0 or not res:
        raise RuntimeError(out.stdout + out.stderr)
    return res


def python_iter_cost(inst, fixes, iters):
    sys.path.insert(0, PROTO)
    import qaputil
    import dnn as proto
    n, a, b, c = qaputil.read_qaplib(dat(inst))
    p = [-1] * n
    for t in fixes.split(","):
        i, l = map(int, t.split(":"))
        p[i] = l
    F, L, AFF, BLL, Cp, const = qaputil.reduce_subproblem(a, b, c, p)
    r = proto.dnn_admm(AFF, BLL, Cp, const, maxit=iters, check_every=iters)
    # one certificate (an extra eigendecomposition) is included, as in C
    return r["cpu"] / r["it"]


def cmd_iter(args):
    out = args.out or os.path.join(ROOT, "bench", "results", "dnn_iter.jsonl")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = load_done(out, lambda r: (r["inst"], r["impl"]))
    for inst, opt, m, fixes in INSTANCES:
        for impl, flags in (("C-dsyevd", ["-e", "0"]), ("C-ssyevd", ["-e", "1"]),
                            ("C-default", []), ("python", None)):
            if (inst, impl) in done:
                continue
            iters = args.iters if flags is not None else max(5, args.iters // 3)
            if flags is None:
                spi = python_iter_cost(inst, fixes, iters)
            else:
                r = dnntest(["-i", str(iters), "-k", str(iters), "-x", fixes] + flags +
                            [dat(inst)])[0]
                spi = float(r["cpu"]) / int(r["iters"])
            append(out, dict(inst=inst, m=m, node=fixes, impl=impl, iters=iters,
                             ms_per_iter=round(1e3 * spi, 2),
                             when=time.strftime("%Y-%m-%d %H:%M")))


def read_state(fn):
    """Nodes of a solver state file: (p, history, depth, bound)."""
    n = None
    nodes = []
    for line in open(fn):
        t = line.split()
        if t and t[0] == "n":
            n = int(t[1])
        elif t and t[0] == "node":
            nodes.append(([int(x) for x in t[5:5 + n]], [int(x) for x in t[5 + n:5 + 2 * n]],
                          int(t[1]), float(t[2])))
    return nodes


def frontier(inst, opt):
    """QPB-open depth-2 nodes of inst as (fixes, own QPB). The depth-2 nodes
    are the children written by -d 1; the QPB-open ones are those that have
    children in the -d 2 output, whose bound field is the parent's QPB."""
    cache = os.path.join(ROOT, "bench", "results", "dnn_frontier_%s.txt" % inst)
    if not os.path.exists(cache):
        with tempfile.TemporaryDirectory() as tmp:
            st = {}
            for d in (1, 2):
                st[d] = os.path.join(tmp, "d%d.state" % d)
                subprocess.run(["./bnbdrv", "-p", "PARAM/rule2.par", "-u", str(opt), "-d",
                                str(d), "-o", st[d], dat(inst)], cwd=LATEST,
                               capture_output=True, check=True)
            qpb = {}
            for p, hist, depth, bound in read_state(st[2]):
                last = hist[depth - 1]
                key = tuple((i, l) for i, l in enumerate(p) if l >= 0 and i != last)
                qpb[key] = bound
            rows = []
            for p, hist, depth, bound in read_state(st[1]):
                key = tuple((i, l) for i, l in enumerate(p) if l >= 0)
                if key in qpb:
                    rows.append("%s %.4f" % (",".join("%d:%d" % kv for kv in key), qpb[key]))
        with open(cache, "w") as f:
            f.write("\n".join(rows) + "\n")
    return [(l.split()[0], float(l.split()[1])) for l in open(cache) if l.strip()]


def cmd_nodes(args):
    out = args.out or os.path.join(ROOT, "bench", "results", "dnn_nodes.jsonl")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = load_done(out, lambda r: (r["inst"], r["node"], r["impl"]))
    insts = [x for x in INSTANCES if not args.inst or x[0] in args.inst]
    for inst, opt, m, _ in insts:
        nodes = frontier(inst, opt)
        rng = random.Random(args.seed)
        sample = rng.sample(nodes, min(args.count, len(nodes)))
        print("%s: %d QPB-open depth-2 nodes, sampling %d" % (inst, len(nodes), len(sample)),
              flush=True)
        impls = [("C-default", [])]
        if args.proto:
            impls.append(("proto", ["-P"]))
        for fx, qpb in sample:
            for impl, flags in impls:
                if (inst, fx, impl) in done:
                    continue
                r = dnntest(["-u", str(opt), "-i", str(args.max_iter), "-x", fx] + flags +
                            [dat(inst)])[0]
                append(out, dict(inst=inst, m=int(r["m"]), node=fx, impl=impl, qpb=qpb,
                                 lb=float(r["lb"]), status=r["status"], iters=int(r["iters"]),
                                 cpu=float(r["cpu"]), cutoff=opt - 1.95,
                                 when=time.strftime("%Y-%m-%d %H:%M")))


def cmd_study(args):
    """Rerun the nodes of a Python study (dnn.py results JSONL: key, it, cpu, lb) with
    the default C bound and the study's cutoff (optimum - 1.95)."""
    recs = [json.loads(l) for l in open(args.study)]
    inst = "nug30" if recs[0]["m"] >= 25 else "nug20"
    opt = dict((x[0], x[1]) for x in INSTANCES)[inst]
    out = args.out or os.path.join(ROOT, "bench", "results", "dnn_study.jsonl")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tag = inst + ":" + os.path.basename(args.study)
    done = load_done(out, lambda r: (r["study"], r["node"]))
    for rec in recs:
        fx = ",".join("%d:%d" % tuple(t) for t in rec["key"])
        if (tag, fx) in done:
            continue
        r = dnntest(["-u", str(opt), "-i", str(args.max_iter), "-x", fx, dat(inst)])[0]
        append(out, dict(study=tag, inst=inst, m=int(r["m"]), node=fx, qpb=rec.get("qpb"),
                         lb=float(r["lb"]), status=r["status"], iters=int(r["iters"]),
                         cpu=float(r["cpu"]), py_lb=rec["lb"], py_iters=rec["it"],
                         py_cpu=rec["cpu"], py_fathom=rec["lb"] > opt - 1.95,
                         when=time.strftime("%Y-%m-%d %H:%M")))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=["iter", "nodes", "study"])
    p.add_argument("--out")
    p.add_argument("--iters", type=int, default=30, help="iter: ADMM iterations per timing")
    p.add_argument("--count", type=int, default=8, help="nodes: nodes per instance")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-iter", type=int, default=4000)
    p.add_argument("--inst", nargs="*", help="nodes: restrict to these instances")
    p.add_argument("--proto", action="store_true", help="nodes: also run prototype mode")
    p.add_argument("--study", help="study: JSONL of a nug20/nug30 study (results_depth*.jsonl)")
    args = p.parse_args()
    if not os.path.exists(DNNTEST):
        sys.exit("build latest/dnntest first: make -C latest dnntest")
    {"iter": cmd_iter, "nodes": cmd_nodes, "study": cmd_study}[args.cmd](args)


if __name__ == "__main__":
    main()
