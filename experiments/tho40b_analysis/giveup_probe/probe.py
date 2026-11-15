"""Two questions for a faster run, on nodes tho40-b actually bounded (Nov 15).

1. The DNN give-up (dnn.giveup = 200): would a give-up node fathom if its DNN
   bound ran on, and would that save more than it costs? In tho40-b's part-way
   export, children of give-ups whose children all fathomed at once hold about
   40% of the CPU (giveups.py).
2. Per-core speed: the time of a DNN iteration with Intel AMX (dnn.amx = 1) and
   without it on this box (a Sapphire Rapids Xeon), against tho40-b's c7a times
   for the same nodes.

Each sampled node is rebuilt as the run bounded it: a chain of bnbdrv processes
from the root down its branching path, each run to its own depth (-d) and each
resumed (-r) from its parent's state file cut to the next node on the path, so
every node starts from its parent's shipped DNN state as a unit does
(path_check/dive.py --warm). Then three arms on the node itself:
  A0  the run's settings: must reproduce the run's DNN line (bound, iterations, status)
  A1  dnn.amx = 1
  G   dnn.giveup = 0 (never; the DNN bound stops at its fathom or its 4,000 iterations)

  sample EXPORT OUT.json [--per-depth 16] [--depths 4,5,6] [--seed 1]
      give-up nodes (DNN status hopeless, branched) drawn uniformly at each depth
      from a qapctl export (nodes.txt.gz, dnn.txt.gz), with what the run did
      below them
  run SAMPLE.json WORKDIR [--workers 4]
      the chains (shared prefixes once), then the arms; resumable: finished
      work is kept and skipped
  report SAMPLE.json WORKDIR
"""
import argparse
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../.."))
LATEST = os.path.join(ROOT, "latest")
sys.path.insert(0, os.path.dirname(HERE))
PAR = os.path.join(ROOT, "infra/phase1/recipes/tho40.par")
SYM = "SYM/tho40t.sym"
INSTANCE = os.path.join(ROOT, "qaplib/tho40t.dat")
UB = 240518
N = 40
ARMS = {"A0": {}, "A1": {"dnn.amx": "1"}, "G": {"dnn.giveup": "0"}}
DNN_RE = re.compile(r"depth (\d+) fixes ([0-9: ]+?) \| first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+)")


def key_of(pairs):
    return tuple(sorted(pairs, key=lambda p: int(p.split(":")[0])))


def sample(a):
    from ntrace import load
    depths = [int(x) for x in a.depths.split(",")]
    dnn = {}
    with gzip.open(os.path.join(a.export, "dnn.txt.gz"), "rt") as f:
        for ln in f:
            m = DNN_RE.search(ln)
            if m and int(m.group(1)) in depths:
                dnn[key_of(m.group(2).split())] = (float(m.group(3)), m.group(4), int(m.group(5)), float(m.group(6)))
    t = load(os.path.join(a.export, "nodes.txt.gz"))
    rng = random.Random(a.seed)
    out = []
    for d in depths:
        cand = [i for i in range(len(t)) if t.depth[i] == d and t.outcome[i] == "B"
                and dnn.get(key_of(t.fixes[i]), (0, ""))[1] == "hopeless"]
        for i in rng.sample(cand, min(a.per_depth, len(cand))):
            lb, status, iters, dcpu = dnn[key_of(t.fixes[i])]
            ch = t.children[i]
            out.append({"depth": d, "fixes": list(t.fixes[i]), "givesups_at_depth": len(cand),
                        "run": {"lb": lb, "status": status, "iters": iters, "dnn_cpu": dcpu, "cpu": t.cpu[i],
                                "cutoff": t.cutoff[i], "kids": t.kids[i], "kids_present": len(ch),
                                "kids_fathomed": sum(1 for c in ch if t.outcome[c] in "QDCP"),
                                "subtree_nodes": t.size[i], "subtree_cpu": t.sub_cpu[i], "complete": t.complete[i]}})
        print(f"depth {d}: {len(cand):,} give-ups, {min(a.per_depth, len(cand))} sampled")
    json.dump(out, open(a.out, "w"), indent=1)


def par_for(arm, work):
    if not ARMS[arm]:
        return PAR
    p = os.path.join(work, f"tho40_{arm}.par")
    if not os.path.exists(p):
        txt = open(PAR).read()
        for k, v in ARMS[arm].items():
            txt = re.sub(rf"^{re.escape(k)} = .*$", f"{k} = {v}", txt, flags=re.M) if re.search(rf"^{re.escape(k)} = ", txt, re.M) \
                else txt.rstrip("\n") + f"\n{k} = {v}\n"
        open(p, "w").write(txt)
    return p


def solve(depth, in_state, out_dir, par, allow=False):
    """bnbdrv on the node of in_state (the root when None), to its own depth; its traces and -o state in out_dir"""
    os.makedirs(out_dir, exist_ok=True)
    nt, dt, st = (os.path.join(out_dir, x) for x in ("nodes.txt", "dnn.txt", "out.state"))
    for f in (nt, dt):
        if os.path.exists(f):
            os.remove(f)
    cmd = ["./bnbdrv", "-p", par, "-s", SYM, "-u", str(UB), "-d", str(depth), "-o", st]
    if in_state:
        cmd += ["-r", in_state]
    if allow:
        cmd.append("--allow-settings-change")
    cmd.append(INSTANCE)
    env = dict(os.environ, QAP_NODE_TRACE=nt, QAP_DNN_TRACE=dt, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with open(os.path.join(out_dir, "log"), "w") as log:
        rc = subprocess.run(cmd, cwd=LATEST, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"bnbdrv exit {rc}: {out_dir}/log")
    node = [l.split() for l in open(nt) if l.startswith("node ") and int(l.split()[1]) == depth]
    dl = [DNN_RE.search(l) for l in open(dt)] if os.path.exists(dt) else []
    dl = [m for m in dl if m and int(m.group(1)) == depth]
    res = {"outcome": node[0][7] if node else None, "cpu": float(node[0][11]) if node else None}
    if dl:
        m = dl[-1]
        res.update(lb=float(m.group(3)), status=m.group(4), iters=int(m.group(5)), dnn_cpu=float(m.group(6)))
    return res


def cut_to(out_state, child, in_state):
    """out_state with the child (facility, location) of the next depth as its only node"""
    f_, l_ = (int(x) for x in child.split(":"))
    kept = 0
    with open(out_state) as fh, open(in_state + ".tmp", "w") as w:
        for line in fh:
            if line.startswith("nodes "):
                continue
            if line.startswith("node "):
                t = line.split()
                if int(t[5 + f_]) == l_ and kept == 0:
                    w.write("nodes 1\n" + line)
                    kept = 1
                continue
            w.write(line)
    if not kept:
        os.remove(in_state + ".tmp")
        raise RuntimeError(f"{out_state}: no child {child}")
    os.replace(in_state + ".tmp", in_state)


def chain_dir(work, prefix):
    return os.path.join(work, "chain", "root" if not prefix else "_".join(p.replace(":", "-") for p in prefix))


def run(a):
    smp = json.load(open(a.sample))
    work = os.path.abspath(a.work)
    os.makedirs(work, exist_ok=True)
    # chains: the in-state of every prefix of every sampled path, depth by depth (shared prefixes once)
    maxd = max(s["depth"] for s in smp)
    for k in range(0, maxd):
        prefixes = sorted({tuple(s["fixes"][:k]) for s in smp if s["depth"] > k})

        def step(prefix):
            d = chain_dir(work, prefix)
            done = os.path.join(d, "done.json")
            if os.path.exists(done):
                return
            in_state = os.path.join(chain_dir(work, prefix[:-1]), "next_" + prefix[-1].replace(":", "-") + ".state") if prefix else None
            r = solve(k, in_state, d, PAR)
            for s in smp:   # cut this node's children to the ones the sampled paths go through
                if tuple(s["fixes"][:k]) == prefix and s["depth"] > k:
                    nxt = s["fixes"][k]
                    cut_to(os.path.join(d, "out.state"), nxt, os.path.join(d, "next_" + nxt.replace(":", "-") + ".state"))
            json.dump(r, open(done, "w"))
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(step, prefixes))
        print(f"chain depth {k}: {len(prefixes)} nodes", flush=True)

    jobs = [(i, arm) for arm in ("A0", "A1", "G") for i in range(len(smp))]

    def arm_job(job):
        i, arm = job
        s = smp[i]
        d = os.path.join(work, "arms", f"{i:03d}_{arm}")
        done = os.path.join(d, "done.json")
        if os.path.exists(done):
            return
        prefix = tuple(s["fixes"][:-1])
        in_state = os.path.join(chain_dir(work, prefix), "next_" + s["fixes"][-1].replace(":", "-") + ".state")
        r = solve(s["depth"], in_state, d, par_for(arm, work), allow=bool(ARMS[arm]))
        json.dump(r, open(done, "w"))
        print(f"{i:03d} {arm} d{s['depth']}: {r}", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        list(ex.map(arm_job, jobs))
    # free the disk the chains used
    shutil.rmtree(os.path.join(work, "chain"), ignore_errors=True) if a.clean else None


def report(a):
    smp = json.load(open(a.sample))
    work = os.path.abspath(a.work)
    rows = []
    for i, s in enumerate(smp):
        r = {arm: json.load(open(p)) for arm in ARMS if os.path.exists(p := os.path.join(work, "arms", f"{i:03d}_{arm}", "done.json"))}
        rows.append((s, r))
    done = [(s, r) for s, r in rows if "A0" in r]
    print(f"{len(done)} of {len(rows)} nodes with A0")
    # 1. reproduction
    same = sum(1 for s, r in done if r["A0"].get("status") == s["run"]["status"] and r["A0"].get("iters") == s["run"]["iters"]
               and abs(r["A0"].get("lb", 0) - s["run"]["lb"]) < 0.01)
    print(f"A0 reproduces the run's DNN line (status, iterations, bound to 0.01): {same} of {len(done)}")
    # 2. per-core speed: ms per DNN iteration
    def ms(x):
        return 1000 * x["dnn_cpu"] / x["iters"] if x.get("iters") else None
    print("\nper-core speed, ms per DNN iteration (median over nodes)")
    for d in sorted({s["depth"] for s, _ in done}):
        sel = [(s, r) for s, r in done if s["depth"] == d and "A1" in r]
        if not sel:
            continue
        run_ms = sorted(1000 * s["run"]["dnn_cpu"] / s["run"]["iters"] for s, _ in sel)
        a0 = sorted(ms(r["A0"]) for _, r in sel)
        a1 = sorted(ms(r["A1"]) for _, r in sel)
        same1 = sum(1 for _, r in sel if r["A1"].get("status") == r["A0"].get("status"))
        med = lambda v: v[len(v) // 2]
        print(f"  depth {d}: c7a in tho40-b {med(run_ms):.1f}, here amx 0 {med(a0):.1f}, amx 1 {med(a1):.1f} "
              f"({len(sel)} nodes; amx 1 the same status on {same1})")
    # 3. the give-up: a limit L in place of the run's give-up. A node that the
    # G arm fathoms by iteration L fathoms under L too (its subtree saved, its own
    # extra iterations paid); one that does not pays for running to L (G's ms per
    # iteration, L less the run's iterations) and keeps its subtree. This box's CPU
    # is put in c7a seconds by each node's A0 against the run.
    print("\nthe give-up: a node that gave up, run on (dnn.giveup = 0) and judged at a limit L (c7a s per give-up)")
    print(f"{'depth':>5} {'nodes':>5} {'L':>5} {'fathom':>6} {'saved':>9} {'extra':>9} {'net':>9} {'x give-ups at depth':>20} {'net, c7a s':>12}")
    totals = {}
    for d in sorted({s["depth"] for s, _ in done}):
        sel = [(s, r) for s, r in done if s["depth"] == d and "G" in r and r["A0"].get("dnn_cpu")]
        if not sel:
            continue
        for L in (400, 1000, 4000):
            saved = extra = 0.0
            nf = 0
            for s, r in sel:
                conv = s["run"]["dnn_cpu"] / r["A0"]["dnn_cpu"]
                g = r["G"]
                per_it = g["dnn_cpu"] / max(1, g["iters"])
                if g.get("outcome") in ("D", "Q", "C") and g["iters"] <= L:
                    nf += 1
                    saved += s["run"]["subtree_cpu"] - s["run"]["cpu"]
                    extra += conv * max(0.0, g["dnn_cpu"] - r["A0"]["dnn_cpu"])
                else:
                    extra += conv * per_it * max(0, min(L, g["iters"]) - s["run"]["iters"])
            k = sel[0][0]["givesups_at_depth"]
            net = (saved - extra) / len(sel)
            totals[L] = totals.get(L, 0.0) + net * k
            print(f"{d:>5} {len(sel):>5} {L:>5} {nf:>6} {saved / len(sel):>9,.0f} {extra / len(sel):>9,.0f} {net:>9,.0f} {k:>20,} {net * k:>12,.0f}")
    for L, v in totals.items():
        print(f"  L = {L}: net {v:,.0f} c7a s over the sampled depths (positive: saved)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("export")
    s.add_argument("out")
    s.add_argument("--per-depth", type=int, default=16)
    s.add_argument("--depths", default="4,5,6")
    s.add_argument("--seed", type=int, default=1)
    r = sub.add_parser("run")
    r.add_argument("sample")
    r.add_argument("work")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--clean", action="store_true")
    p = sub.add_parser("report")
    p.add_argument("sample")
    p.add_argument("work")
    a = ap.parse_args()
    {"sample": sample, "run": run, "report": report}[a.cmd](a)


if __name__ == "__main__":
    main()
