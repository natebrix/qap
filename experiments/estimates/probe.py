#!/usr/bin/env python3
"""A sampled estimate of subtrees too large to finish, under an instance's settings file.

  SAMPLING_SETTINGS=<settings.sh> SAMPLING_PARAMS=<file.par> \\
  probe.py --run DIR [--k 24] [--beta 0.02] [--slots 4] [--cpus 0,1,2,3] [--seed 1] NAME=SPEC ...

Each NAME=SPEC is one node (a -F spec, i:j,...; its depth is the number of fixes).
The probe goes down the node's subtree level by level until nothing is open:

  the node       -F SPEC -d L0, cold, as subtree.sh starts a subtree
  each level     the open children of the nodes bounded at the level above. When
                 there are at most --k of them, all are bounded. Otherwise --k draws
                 by systematic sampling, the children ordered by their parent's
                 shortfall (implicit strata), the parent and the child's predicted
                 bound, with probability proportional to weight x a: weight the
                 number of nodes the parent stands for, a = exp(beta x the parent's
                 shortfall at its best bound), since a subtree's size grows about
                 exponentially in its shortfall and the draws should go where the
                 tree is. A drawn child stands for hits x S / (k a) nodes, S the
                 level's sum of weight x a: the Hansen-Hurwitz estimator, unbiased
                 for each level's node count and CPU whatever beta is (on synthetic
                 trees, 60 probes each: mean 0.99-1.02 of the truth; one probe's sd
                 0.67 of the truth at beta 0 and 0.36 at a beta matched to the tree).

Every bounded node is a one-node task (-r, its parent's warm start with the
shipped DNN state and the cut pool), as a fleet run hands nodes out, pinned to one
CPU. The probes given share the slots; each goes level by level.

Output per probe, in DIR/NAME/:
  summary.txt  first line in subtree_summary.py's form (estimated: "N DNN nodes ...
               DNN CPU X s"), so estimate.py reads it as an anchor; then the levels
  nodes.tsv    every bounded node: level, weight, shortfalls, phase, CPU, children
The SE is the ultimate-cluster one (the draws of the first sampled level with
everything below them, as independent draws): approximate, and on the synthetic
trees a 95% interval of +-1.96 SE held the truth 75-90% of the time. With no sampled level
the probe is the complete subtree, run node by node, and the SE is 0.

Resumable: a finished task (result.json) is kept and each level's draws are saved;
run the same command again. SIGTERM or SIGINT drains: no new task starts, running
ones finish. Disk: a node with children keeps its out.state (its DNN state, about
25 MB at m = 46) until its drawn children are done."""
import argparse
import bisect
import json
import math
import os
import queue
import random
import signal
import sys
import threading
import time

H = os.path.dirname(os.path.abspath(__file__))
R = os.path.abspath(os.path.join(H, "../.."))
sys.path.insert(0, os.path.join(R, "experiments", "wil50_feasibility", "sampling"))
sys.path.insert(0, os.path.join(R, "bench"))
from common import (log_line, fullrun_env, params_with, dnn_flags, solver_env, private_binary,  # noqa: E402
                    node_fields, spec as spec_of, write_one, run_bnbdrv, parse_trace, stats_of, write_json,
                    disk_free_gb, PARAMS, SETTINGS)
from parallel import State  # noqa: E402


class Stop(Exception):
    pass


class Runner:
    """the shared slots: one CPU per running task"""

    def __init__(self, a):
        self.a = a
        self.run = os.path.abspath(a.run)
        os.makedirs(self.run, exist_ok=True)
        self.log_path = os.path.join(self.run, "probe.log")
        self.binary, self.commit = private_binary(a.binary, os.path.join(self.run, "bin"))
        settings, info = fullrun_env(self.binary)
        params = os.path.join(self.run, "params.par")
        params_with(PARAMS, {}, params, self.binary)
        cfg = {"instance": info["FULLRUN_INSTANCE"], "sym": info["FULLRUN_SYM"], "ub": float(info["FULLRUN_UB"]),
               "params": params, "params_base": PARAMS, "settings_file": SETTINGS, "settings_env": settings,
               "k": a.k, "beta": a.beta, "seed": a.seed, "commit": self.commit}
        cp = os.path.join(self.run, "config.json")
        if os.path.exists(cp):
            old = json.load(open(cp))
            for key in ("instance", "sym", "ub", "params_base", "settings_env", "k", "beta", "seed"):
                if old.get(key) != cfg[key]:
                    raise SystemExit(f"{cp}: {key} was {old.get(key)!r}, now {cfg[key]!r}; use another --run")
            cfg = old
        else:
            write_json(cp, cfg)
        self.cfg = cfg
        self.env = solver_env(cfg["settings_env"])
        self.cpus = queue.Queue()
        for c in [int(x) for x in a.cpus.split(",")][:a.slots]:
            self.cpus.put(c)
        self.stop = threading.Event()
        self.procs = {}

    def log(self, msg):
        log_line(self.log_path, msg)

    def task(self, wd, level, F=None, src=None, index=None):
        """bound one node to `level` (its children open in out.state); result.json"""
        res_p = os.path.join(wd, "result.json")
        if os.path.exists(res_p):
            return json.load(open(res_p))
        if self.stop.is_set():
            raise Stop()
        os.makedirs(wd, exist_ok=True)
        while disk_free_gb(wd) < self.a.min_free_gb:
            if self.stop.is_set():
                raise Stop()
            time.sleep(30)
        cpu = self.cpus.get()
        try:
            if self.stop.is_set():
                raise Stop()
            c = self.cfg
            out, inp = os.path.join(wd, "out.state"), os.path.join(wd, "in.state")
            args = ["-p", c["params"]] + dnn_flags("") + ["-s", c["sym"], "-u", repr(c["ub"]), "-d", str(level),
                                                          "-o", out]
            if F:
                args += ["-F", F]
            else:
                write_one(src, index, inp, warm="keep")
                args += ["-r", inp]
            args.append(c["instance"])
            trace, logf = os.path.join(wd, "trace.txt"), os.path.join(wd, "solver.log")
            rc, wall = run_bnbdrv(self.binary, args, self.env, logf, trace, cpu=cpu, procs=self.procs, tag=wd)
        finally:
            self.cpus.put(cpu)
        if os.path.exists(inp):
            os.remove(inp)
        if rc != 0 or not os.path.exists(out):
            raise RuntimeError(f"{wd}: bnbdrv exit {rc} (see {logf})")
        st = State(out)
        kids = []
        for text, _ in st.nodes:
            d, b, pb, f = node_fields(text, st.n)
            if d != level + 1:
                raise RuntimeError(f"{wd}: a child at depth {d}, expected {level + 1}")
            kids.append({"spec": spec_of(f), "bound": b, "pred": pb})
        rows = parse_trace(trace)
        dnn = sum(r["cpu"] + r.get("cut_cpu", 0.0) for r in rows)
        doc = {"level": level, "wall": wall, "cpu_slot": cpu, "stats": stats_of(st), "trace": rows,
               "dnn_cpu": dnn, "task_cpu": st.times[0], "children": kids, "obj": st.obj,
               "finished": time.strftime("%F %T")}
        if st.obj < c["ub"] - 0.5:
            doc["new_incumbent"] = st.obj
        write_json(res_p, doc)
        if not kids:
            os.remove(out)          # fathomed: nothing below needs its state
        return doc


def systematic(weights, k, rng):
    """k systematic draws with probability proportional to weight: hits per unit"""
    cum, s = [], 0.0
    for w in weights:
        s += w
        cum.append(s)
    step = s / k
    u = rng.random() * step
    hits = [0] * len(weights)
    for m in range(k):
        hits[min(len(weights) - 1, bisect.bisect_right(cum, u + m * step))] += 1
    return s, hits


class Probe:
    def __init__(self, runner, name, spec):
        self.r, self.name, self.spec = runner, name, spec
        self.dir = os.path.join(runner.run, name)
        self.L0 = len(spec.split(","))
        self.nodes = []        # one dict per bounded node: level, id, parent, weight, hits, sampled, doc

    def tdir(self, level, j):
        return os.path.join(self.dir, f"L{level}", f"t{j}")

    def go(self):
        r = self.r
        os.makedirs(self.dir, exist_ok=True)
        L = self.L0
        doc = r.task(self.tdir(L, 0), L, F=self.spec)
        level = [{"level": L, "id": (L, 0), "parent": None, "weight": 1.0, "hits": 1, "sampled": False, "doc": doc}]
        self.nodes += level
        r.log(f"{self.name}: level {L}, the node: {self.describe(level)}")
        while True:
            units = []          # (parent's index in the level, child index, the child's predicted bound)
            for pi, p in enumerate(level):
                for i, kid in enumerate(p["doc"]["children"]):
                    units.append((pi, i, kid["pred"]))
            short = [self.short(p) for p in level]
            if not units:
                break
            L += 1
            if L > self.L0 + r.a.max_levels:
                raise RuntimeError(f"{self.name}: still open at level {L}")
            dpath = os.path.join(self.dir, f"L{L}", "draws.json")
            if os.path.exists(dpath):
                draws = json.load(open(dpath))
            else:
                k = r.a.k
                # by the parent's shortfall (implicit strata), parent, the child's predicted bound
                units.sort(key=lambda u: (short[u[0]], u[0], u[2]))
                if len(units) <= k:
                    draws = [{"p": u[0], "i": u[1], "hits": 1, "weight": level[u[0]]["weight"]} for u in units]
                    sampled = False
                else:
                    rng = random.Random(f"{r.cfg['seed']}|{self.name}|{L}")
                    # importance: a child's draw probability is its weight times exp(beta x its
                    # parent's shortfall), since a subtree's size grows about exponentially in it;
                    # a draw then stands for (sum of weight x a) / (k a) nodes (Hansen-Hurwitz)
                    b = r.cfg["beta"]
                    a = [math.exp(b * (short[u[0]] - max(short))) for u in units]
                    S, hits = systematic([level[u[0]]["weight"] * x for u, x in zip(units, a)], k, rng)
                    draws = [{"p": u[0], "i": u[1], "hits": h, "weight": h * S / (k * x)}
                             for u, x, h in zip(units, a, hits) if h]
                    sampled = True
                os.makedirs(os.path.dirname(dpath), exist_ok=True)
                write_json(dpath, {"sampled": sampled, "units": len(units), "draws": draws})
                draws = json.load(open(dpath))
            sampled = draws["sampled"]
            # parents no draw points at: their DNN state is not needed
            used = {d["p"] for d in draws["draws"]}
            for pi, p in enumerate(level):
                if pi not in used:
                    self.drop_state(p)
            nxt = [None] * len(draws["draws"])
            errs = []

            def one(j, d):
                p = level[d["p"]]
                src = os.path.join(self.tdir(*p["id"]), "out.state")
                doc = r.task(self.tdir(L, j), L, src=src, index=d["i"])
                nxt[j] = {"level": L, "id": (L, j), "parent": p["id"], "weight": d["weight"], "hits": d["hits"],
                          "sampled": sampled, "doc": doc}
            th = [threading.Thread(target=lambda j=j, d=d: self.guard(one, j, d, errs)) for j, d in
                  enumerate(draws["draws"])]
            for t in th:
                t.start()
            for t in th:
                t.join()
            if errs:
                raise errs[0]
            for p in level:
                self.drop_state(p)
            level = nxt
            self.nodes += level
            r.log(f"{self.name}: level {L}, {draws['units']:,} open, {len(level)} bounded"
                  f"{' (sampled)' if sampled else ''}: {self.describe(level)}")
        self.summarize(final=True)

    @staticmethod
    def guard(fn, j, d, errs):
        try:
            fn(j, d)
        except Exception as e:   # noqa: BLE001  (a drain's Stop included: the level is retried on resume)
            errs.append(e)

    @staticmethod
    def short(p):
        """a bounded node's shortfall at its best bound (after its phase)"""
        t = p["doc"]["trace"]
        return max(0.0, t[0]["short"]) if t else 0.0

    def drop_state(self, p):
        f = os.path.join(self.tdir(*p["id"]), "out.state")
        if os.path.exists(f):
            os.remove(f)

    @staticmethod
    def describe(level):
        fath = sum(1 for n in level if not n["doc"]["children"])
        kids = sum(len(n["doc"]["children"]) for n in level)
        cpu = sum(n["doc"]["dnn_cpu"] for n in level)
        ph = sum(1 for n in level for t in n["doc"]["trace"] if t["cut"])
        return f"{fath} fathomed, {ph} phases, {kids} children, DNN CPU {cpu:,.0f} s"

    def summarize(self, final=False):
        r = self.r
        cutoff = None
        by = {}
        for n in self.nodes:
            by.setdefault(n["level"], []).append(n)
        tot = {"nodes": 0.0, "dnn": 0.0, "cpu": 0.0, "task": 0.0}
        lines = []
        rows = ["level\tweight\thits\tparent\tfirst_short\tbest_short\tphase\tfathomed\tdnn_cpu\ttask_cpu\tchildren"]
        for L in sorted(by):
            lv = by[L]
            W = sum(n["weight"] for n in lv)
            dn = sum(n["weight"] * len(n["doc"]["trace"]) for n in lv)
            cpu = sum(n["weight"] * n["doc"]["dnn_cpu"] for n in lv)
            task = sum(n["weight"] * n["doc"]["task_cpu"] for n in lv)
            tot["nodes"] += W
            tot["dnn"] += dn
            tot["cpu"] += cpu
            tot["task"] += task
            fath = sum(1 for n in lv if not n["doc"]["children"])
            ph = sum(1 for n in lv for t in n["doc"]["trace"] if t["cut"])
            sh = sorted(t["short"] for n in lv for t in n["doc"]["trace"])
            med = sh[len(sh) // 2] if sh else float("nan")
            lines.append(f"| {L} | {W:,.1f} | {len(lv)} | {'yes' if lv[0]['sampled'] else 'no'} | {fath} | {ph} | "
                         f"{med:.1f} | {cpu / max(W, 1e-9):,.1f} | {cpu:,.0f} |")
            for n in lv:
                t = n["doc"]["trace"][0] if n["doc"]["trace"] else None
                if t:
                    cutoff = t["cutoff"]
                rows.append("\t".join(map(str, [L, round(n["weight"], 4), n["hits"], n["parent"],
                                                round(t["cutoff"] - t["lb"], 2) if t else "",
                                                round(t["short"], 2) if t else "", t["cut"] if t else "",
                                                int(not n["doc"]["children"]), round(n["doc"]["dnn_cpu"], 1),
                                                round(n["doc"]["task_cpu"], 1), len(n["doc"]["children"])])))
        se_n, se_c, Ls = self.se()
        k = r.cfg["k"]
        root = self.nodes[0]["doc"]["trace"]
        rt = root[0] if root else None
        exact = Ls is None
        head = (f"{self.name}: {tot['dnn']:,.0f} DNN nodes{'' if exact else ' (estimated)'}, "
                f"{tot['nodes']:,.0f} nodes in all, DNN CPU {tot['cpu']:,.0f} s"
                + ("" if exact else f" (SE {se_c:,.0f} s; nodes SE {se_n:,.0f})")
                + f", task CPU {tot['task']:,.0f} s; {len(self.nodes)} nodes bounded"
                + ("" if final else " (unfinished)"))
        txt = [head,
               f"  root (depth {self.L0}) {self.spec}: first solve "
               + (f"{rt['cutoff'] - rt['lb']:.1f} short after {rt['iters']} iterations"
                  + (f", phase to {rt['short']:.1f} short" if rt["cut"] else ", no phase") if rt else "no DNN line"),
               "  " + ("complete: every node bounded" if exact else f"sampled from level {Ls}, k = {k} draws a level")
               + f"; settings {os.path.relpath(r.cfg['settings_file'], R)}, binary {r.cfg['commit']['commit'][:10]}",
               "", "| level | nodes (est.) | bounded | sampled | fathomed | phases | median short | DNN CPU s per node | "
               "DNN CPU s at the level (est.) |", "|---|---|---|---|---|---|---|---|---|"] + lines
        if any("new_incumbent" in n["doc"] for n in self.nodes):
            txt.append("\nNEW INCUMBENT found: " + ", ".join(str(n["doc"]["new_incumbent"]) for n in self.nodes
                                                          if "new_incumbent" in n["doc"]))
        open(os.path.join(self.dir, "summary.txt"), "w").write("\n".join(txt) + "\n")
        open(os.path.join(self.dir, "nodes.tsv"), "w").write("\n".join(rows) + "\n")
        return "\n".join(txt)

    def se(self):
        """ultimate clusters: the nodes of the first sampled level and everything
        below them; (SE nodes, SE DNN CPU, that level) or (0, 0, None)"""
        Ls = min((n["level"] for n in self.nodes if n["sampled"]), default=None)
        if Ls is None:
            return 0.0, 0.0, None
        k = self.r.cfg["k"]
        anc = {}
        for n in self.nodes:
            if n["level"] == Ls:
                anc[n["id"]] = n["id"]
            elif n["level"] > Ls:
                anc[n["id"]] = anc[tuple(n["parent"])]
        zn, zc = {}, {}
        for n in self.nodes:
            if n["level"] >= Ls:
                a = anc[n["id"]]
                zn[a] = zn.get(a, 0.0) + n["weight"]
                zc[a] = zc.get(a, 0.0) + n["weight"] * n["doc"]["dnn_cpu"]
        hits = {n["id"]: n["hits"] for n in self.nodes if n["level"] == Ls}

        def var(z):
            t = [k * z[a] / hits[a] for a in z for _ in range(hits[a])]
            m = sum(t) / len(t)
            return sum((x - m) ** 2 for x in t) / (len(t) - 1) / len(t) if len(t) > 1 else 0.0
        return math.sqrt(var(zn)), math.sqrt(var(zc)), Ls


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--k", type=int, default=24, help="draws a level once more than k children are open")
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--cpus", default="0,1,2,3")
    ap.add_argument("--beta", type=float, default=0.02,
                    help="importance: draw probability x exp(beta x the parent's shortfall); 0 = by weight alone")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-levels", type=int, default=30)
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--binary", default=os.path.join(R, "latest", "bnbdrv"))
    ap.add_argument("nodes", nargs="+", help="NAME=SPEC")
    a = ap.parse_args()
    r = Runner(a)

    def drain(sig, frm):
        if r.stop.is_set():
            for p in list(r.procs.values()):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except OSError:
                    pass
        r.stop.set()
        r.log("draining: no new task starts")
    signal.signal(signal.SIGTERM, drain)
    signal.signal(signal.SIGINT, drain)
    probes = []
    for x in a.nodes:
        name, _, sp = x.partition("=")
        probes.append(Probe(r, name, sp))
    r.log(f"probe: {len(probes)} subtrees, k = {r.cfg['k']}, {a.slots} slots; settings {SETTINGS}")
    errs = []

    def run(p):
        try:
            p.go()
            r.log(f"{p.name}: done\n" + p.summarize(final=True))
        except Stop:
            p.summarize()
            r.log(f"{p.name}: stopped (drain); run again to resume")
        except Exception as e:   # noqa: BLE001
            errs.append((p.name, e))
            r.log(f"{p.name}: FAILED {e!r}")
    th = [threading.Thread(target=run, args=(p,)) for p in probes]
    for t in th:
        t.start()
    for t in th:
        t.join()
    sys.exit(1 if errs else (3 if r.stop.is_set() else 0))


if __name__ == "__main__":
    main()
