#!/usr/bin/env python3
"""The node trace (QAP_NODE_TRACE, format in latest/bnb.C): a line per node
the search takes off the stack, written by the process that decides it.

  python3 tests/node_trace_test.py [--only serial,killed,parallel,coord,coordkill]

serial     nug12 QPB with the optimum as upper bound and without one (the
           second has nodes discarded when popped and enumerated nodes),
           nug16a with the DNN bound at depths 2-4 (QPB fathoms above and
           below), nug17 with the DNN bound, the cut phase and the pool
           (tests/par/ship.par). Each:
           * the run with the trace makes the statistics of the run without
             (state files compared, CPU times aside)
           * lines = node counter + popped: the counter (stats.C anotherNode)
             counts every node bounded or enumerated; a node discarded when
             popped (anotherDiscard) is counted in fathom and pop_f_level only
           * by depth, the lines other than P = n_level, and Q, D, C, P = the
             fathoms by cause (f_level less the rest, dnn_f_level, cut_f_level,
             pop_f_level); E = the enumerate counter
           * every node once; the parent of each non-root node (its fixes less
             the last pair) is a B line, whose kids are its children in the
             trace and whose line (r<facility>, c<location>) each child fixes
           * the values agree with the outcome (a fathomed node's bound above
             the cutoff, a branched node's not)
killed     a solver SIGKILLed after its -k checkpoints (nug15): its trace cut to
           the lines the last checkpoint accounts for (coord/worker.py
           traced_nodes and keep_lines, as for a killed unit), then the trace of
           a run resumed from that checkpoint, are the serial run's
parallel   bench/parallel.py (nug15, --jobs 2, the optimum as upper bound)
           passes QAP_NODE_TRACE to the frontier pass and every worker: their
           lines, appended to one file, are the serial run's
coord      the coordinator on the local backend with --trace, nug17 with the
           DNN (ship.bits = 64, so the tree is the serial one) in 0.3-second
           units: frontier/nodes.txt plus the nodes.txt of each unit's
           recorded attempt are the serial run's lines (fixes and outcome),
           each node once
coordkill  the same with nug17 QPB in 4-second units, the worker SIGKILLed once
           a unit has uploaded a checkpoint: the new worker resumes it from the
           checkpoint's trace (ckpt.nodes.txt)
           Both runs are public and publish every second, so the live tree
           (coord/tree.py, public/runs/<run>/tree/) grows as units merge. Read
           from 0.json through the boundary rows as the viewer reads it, the
           final tree is the serial run's: every node once as a child row, its
           parent, outcome, subtree size and fathoms. qapctl publish then
           rebuilds it from the traces (as a restarted coordinator does) into
           the same files: the same ids.

At most two solver processes run at once.
"""
import argparse
import collections
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
os.environ.setdefault("QAP_MIN_FREE_GB", "0.5")   # the worker's disk guard: a test machine's disk is small

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "latest", "bnbdrv")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench"))
from parallel import State  # noqa: E402  (latest/state.h)
import envpar  # noqa: E402

ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ENV.pop("QAP_NODE_TRACE", None)
DAT = lambda name: os.path.join(ROOT, "qaplib", name + ".dat")
SHIP_PAR = os.path.join(ROOT, "tests", "par", "ship.par")
# (name, instance, solver flags, environment)
SERIAL = [
    ("nug12", "nug12", ["-p", "PARAM/default.par", "-u", "578"], {}),
    ("nug12-noub", "nug12", ["-p", "PARAM/default.par"], {}),
    ("nug16a-dnn", "nug16a", ["-p", "../tests/par/default_ship64.par", "-u", "1610", "-x", "2:4:1:0:1"], {}),
    ("nug17-cut", "nug17", ["-p", SHIP_PAR, "-u", "1732", "-x", "1:3:1:20:1"], {"QAP_DNN_POOL_CHECK": "20"}),
]
TIMES = ("t_level", "t_strat", "dnn_time", "dnn_t_level", "cut_t_level")
TOL = 1e-3   # the trace prints 4 decimals


def solve(flags, dat, state, trace=None, env=()):
    e = dict(ENV, **dict(env))
    if trace:
        e["QAP_NODE_TRACE"] = trace
    r = subprocess.run([SOLVER] + flags + ["-o", state, dat], cwd=os.path.dirname(SOLVER),
                       capture_output=True, text=True, env=e, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"bnbdrv {' '.join(flags)} exited {r.returncode}:\n{r.stdout[-1500:]}{r.stderr}")
    return r.stdout


# ---------------------------------------------------------------- the trace

def num(x):
    return None if x == "-" else float(x)


def read_trace(*paths):
    rows = []
    for p in paths:
        with open(p) as f:
            for line in f:
                t = line.split()
                if len(t) != 12 or t[0] != "node":
                    raise ValueError(f"{p}: not a node line: {line!r}")
                rows.append({"depth": int(t[1]), "fixes": t[2], "qpb": num(t[3]), "dnn": num(t[4]),
                             "cut": num(t[5]), "cutoff": num(t[6]), "outcome": t[7], "rule": t[8],
                             "line": t[9], "kids": t[10], "cpu": float(t[11])})
    return rows


def pairs(fixes):
    return [] if fixes == "-" else [tuple(map(int, p.split(":"))) for p in fixes.split(",")]


def parent(fixes):
    return ",".join(fixes.split(",")[:-1]) or "-"


def tree_problems(rows):
    """structure and values, from the trace alone"""
    bad = []
    by_fixes = {}
    for r in rows:
        if r["fixes"] in by_fixes:
            bad.append(f"node {r['fixes']} twice")
        by_fixes[r["fixes"]] = r
    children = collections.defaultdict(list)
    for r in rows:
        p = pairs(r["fixes"])
        if r["depth"] != len(p):
            bad.append(f"node {r['fixes']}: depth {r['depth']}")
        if r["fixes"] != "-":
            children[parent(r["fixes"])].append(r)
        o, vals = r["outcome"], [v for v in (r["qpb"], r["dnn"], r["cut"]) if v is not None]
        if o in "PE" and vals:
            bad.append(f"node {r['fixes']}: {o} with a bound")
        if o in "QDC" and (r["cutoff"] is None or r["qpb"] is None):
            bad.append(f"node {r['fixes']}: {o} without a cutoff or a QPB value")
        elif o == "Q" and not max(vals) > r["cutoff"] - TOL:
            bad.append(f"node {r['fixes']}: Q with every bound below the cutoff")
        elif o == "D" and not (r["dnn"] is not None and r["dnn"] > r["cutoff"] - TOL and r["cut"] is None):
            bad.append(f"node {r['fixes']}: D, but the DNN value {r['dnn']} and cut {r['cut']}")
        elif o == "C" and not (r["cut"] is not None and r["cut"] > r["cutoff"] - TOL):
            bad.append(f"node {r['fixes']}: C, but the cut value {r['cut']}")
        if o == "B":
            if r["cutoff"] is not None and vals and max(vals) > r["cutoff"] + TOL:
                bad.append(f"node {r['fixes']}: B with a bound above the cutoff")
            if r["line"][0] not in "rc" or not r["rule"].isdigit() or not r["kids"].isdigit():
                bad.append(f"node {r['fixes']}: B with rule {r['rule']} line {r['line']} kids {r['kids']}")
        elif (r["rule"], r["line"], r["kids"]) != ("-", "-", "-"):
            bad.append(f"node {r['fixes']}: {o} with a rule, line or kids")
    roots = [r for r in rows if r["fixes"] == "-"]
    if len(roots) != 1:
        bad.append(f"{len(roots)} root lines")
    for p, kids in children.items():
        if p not in by_fixes or by_fixes[p]["outcome"] != "B":
            bad.append(f"parent {p} of {len(kids)} node(s) is {'missing' if p not in by_fixes else by_fixes[p]['outcome']}")
    for r in rows:
        if r["outcome"] != "B":
            continue
        kids = children.get(r["fixes"], [])
        if int(r["kids"]) != len(kids):
            bad.append(f"node {r['fixes']}: kids {r['kids']}, {len(kids)} in the trace")
        side, idx = r["line"][0], int(r["line"][1:])
        for k in kids:
            f, l = pairs(k["fixes"])[-1]
            if (f if side == "r" else l) != idx:
                bad.append(f"node {k['fixes']}: does not fix line {r['line']} of its parent")
    return bad


def stats_problems(rows, st):
    """the trace against the statistics of the state file the run wrote"""
    bad = []
    node, fathom, enum, _ = st.counters
    a = st.arrays
    popped = sum(a["pop_f_level"])
    if len(rows) != node + popped:
        bad.append(f"{len(rows)} lines, node counter {node} + popped {popped}")
    n = len(a["n_level"])
    count = collections.Counter((min(r["depth"], n - 1), r["outcome"]) for r in rows)
    for d in range(n):
        want = {"P": a["pop_f_level"][d], "D": a["dnn_f_level"][d], "C": a["cut_f_level"][d],
                "Q": a["f_level"][d] - a["pop_f_level"][d] - a["dnn_f_level"][d] - a["cut_f_level"][d]}
        got = {o: count[(d, o)] for o in want}
        evaluated = sum(count[(d, o)] for o in "QDCEB")
        if got != want or evaluated != a["n_level"][d]:
            bad.append(f"depth {d}: trace {got} and {evaluated} evaluated, state {want} and n_level {a['n_level'][d]}")
    e = sum(1 for r in rows if r["outcome"] == "E")
    if e != enum:
        bad.append(f"{e} E lines, enumerate counter {enum}")
    return bad


def same_stats(s1, s2):
    """the differences between two state files' statistics, CPU times aside"""
    d = []
    if s1.counters != s2.counters:
        d.append(f"counters {s1.counters} / {s2.counters}")
    if s1.obj != s2.obj:
        d.append(f"incumbent {s1.obj} / {s2.obj}")
    for k in s1.arrays:
        if k not in TIMES and s1.arrays[k] != s2.arrays.get(k):
            d.append(k)
    return d


def outcome_counts(rows):
    c = collections.Counter(r["outcome"] for r in rows)
    return " ".join(f"{o} {c[o]}" for o in "QDCPEB" if c[o])


# ---------------------------------------------------------------- serial

def serial(tmp):
    bad = 0
    for tag, inst, flags, env in SERIAL:
        trace, s_on, s_off = (os.path.join(tmp, f"{tag}.{x}") for x in ("nodes", "on.state", "off.state"))
        out = {}
        # with and without the trace at once (two processes)
        th = threading.Thread(target=lambda: out.setdefault("off", solve(flags, DAT(inst), s_off, None, env)))
        th.start()
        out["on"] = solve(flags, DAT(inst), s_on, trace, env)
        th.join()
        rows, on, off = read_trace(trace), State(s_on), State(s_off)
        problems = same_stats(on, off)
        problems = ([f"the trace changed the statistics: {', '.join(problems)}"] if problems else []) \
            + stats_problems(rows, on) + tree_problems(rows)
        bad += bool(problems)
        node, fathom, enum, fw = on.counters
        print(f"{'ok  ' if not problems else 'FAIL'} {tag}: {len(rows)} lines = node counter {node} + popped "
              f"{sum(on.arrays['pop_f_level'])}; {outcome_counts(rows)}; fathom {fathom}, enumerate {enum}, "
              f"fw {fw} with and without the trace")
        for p in problems[:10]:
            print(f"     {p}")
    return bad


# ---------------------------------------------------------------- a killed solver

def killed(tmp):
    from coord.worker import keep_lines, traced_nodes
    flags = ["-p", "PARAM/default.par", "-u", "1150"]
    ref, first, second, ck = (os.path.join(tmp, f"k15.{x}") for x in ("serial.nodes", "1.nodes", "2.nodes", "ckpt.state"))
    solve(flags, DAT("nug15"), os.path.join(tmp, "k15.state"), ref)
    p = subprocess.Popen([SOLVER] + flags + ["-k", "0.05", "-o", ck, DAT("nug15")], cwd=os.path.dirname(SOLVER),
                         env=dict(ENV, QAP_NODE_TRACE=first), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t_end = time.time() + 30
    while time.time() < t_end and p.poll() is None and not os.path.exists(ck):
        time.sleep(0.02)
    time.sleep(0.2)                     # past that checkpoint
    p.kill()
    p.wait()
    total = sum(1 for _ in open(first))
    kept = keep_lines(first, first, traced_nodes(ck))
    solve(["-p", "PARAM/default.par", "-r", ck], DAT("nug15"), os.path.join(tmp, "k15.end.state"), second)
    bad = compare("killed nug15", read_trace(ref), read_trace(first, second),
                  f"{kept} of the killed run's {total} lines, then the resumed run's")
    if kept == total:
        print("     (the killed run wrote no line after its last checkpoint: the cut was not exercised)")
    return bad


# ---------------------------------------------------------------- bench/parallel.py

def parallel(tmp):
    flags = ["-p", "PARAM/default.par", "-u", "1150"]
    ref = os.path.join(tmp, "p15.serial.nodes")
    solve(flags, DAT("nug15"), os.path.join(tmp, "p15.state"), ref)
    shared = os.path.join(tmp, "p15.parallel.nodes")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bench", "parallel.py"), DAT("nug15"), "--jobs", "2",
                        "--ub", "1150"], capture_output=True, text=True, env=dict(ENV, QAP_NODE_TRACE=shared))
    if r.returncode != 0:
        print(f"FAIL parallel.py exited {r.returncode}: {r.stderr[-1000:]}")
        return 1
    return compare("parallel.py nug15 --jobs 2", read_trace(ref), read_trace(shared), "one file, every process")


def compare(what, ref, got, how):
    key = lambda rows: collections.Counter((r["fixes"], r["outcome"]) for r in rows)
    k_ref, k_got = key(ref), key(got)
    twice = sum(1 for c in k_got.values() if c > 1)
    ok = k_ref == k_got and twice == 0
    print(f"{'ok  ' if ok else 'FAIL'} {what}: {len(got)} lines ({how}), serial {len(ref)}; "
          f"{outcome_counts(got)}; {twice} nodes twice"
          + ("" if ok else f"; {len(k_ref - k_got)} serial lines missing, {len(k_got - k_ref)} extra"))
    if not ok:
        for (f, o), c in list((k_ref - k_got).items())[:5]:
            print(f"     missing {f} {o}")
        for (f, o), c in list((k_got - k_ref).items())[:5]:
            print(f"     extra   {f} {o} x{c}")
    return 0 if ok else 1


# ---------------------------------------------------------------- the coordinator

def coord(tmp, kill, dnn=None):
    from coord.backend import open_backend
    dnn = (not kill) if dnn is None else dnn
    if dnn:
        # ship.bits = 64: a node resumed from a state file starts its DNN from its
        # parent's exact state, which goes through the parent-state store (coord/warm.py)
        params = os.path.join(tmp, "ship64.par")
        envpar.write_par(SHIP_PAR, {"ship.bits": "64"}, params, solver=SOLVER)
        inst, ub, extra, env = "nug17", 1732, "-x 1:3:1:20:1", {"QAP_DNN_POOL_CHECK": "20"}
    else:
        inst, ub, params, extra, env = "nug17", 1732, "PARAM/default.par", "", {}
    if kill:
        tag = "coordkilldnn" if dnn else "coordkill"
        timing = ["--budget", "4", "--target-unit-cpu", "4", "--steal-after", "3"]
    else:
        tag = "coord"
        timing = ["--budget", "0.3", "--target-unit-cpu", "0.5", "--steal-after", "0.5"]
    ref = os.path.join(tmp, f"{tag}.serial.nodes")
    solve(["-p", params, "-u", str(ub)] + extra.split(), DAT(inst), os.path.join(tmp, f"{tag}.state"), ref, env)
    be = f"local:{os.path.join(tmp, tag + '.store')}"
    ctl = lambda *a, log: subprocess.Popen([sys.executable, "-m", "coord.qapctl"] + list(a), cwd=ROOT, env=ENV,
                                           stdout=open(os.path.join(tmp, f"{tag}.{log}.log"), "a"),
                                           stderr=subprocess.STDOUT, start_new_session=True)
    init = ["init", tag, "--backend", be, "--instance", DAT(inst), "--params", params, "--ub", str(ub),
            "--extra", extra, "--trace", "--publish-every", "1", "--slots", "2", "--lease", "4",
            "--ckpt-upload", "0.5", "--heartbeat", "0.5", "--checkpoint-every", "0.2", "--poll", "0.3",
            "--status-every", "5"] + timing + [f"--env={k}={v}" for k, v in env.items()]
    if ctl(*init, log="init").wait() != 0:
        print(f"FAIL {tag}: init failed (logs in {tmp})")
        return 1
    worker = lambda w: ctl("worker", tag, "--backend", be, "--slots", "2", "--workdir", os.path.join(tmp, f"{tag}.{w}"),
                           log=w)
    c, w = ctl("coordinator", tag, "--backend", be, "--poll", "0.3", log="coordinator"), worker("w1")
    b = open_backend(be)
    # the parent states in the store (coord/warm.py) while the run goes on
    seen = {"max": 0}
    def watch():
        while c.poll() is None:
            seen["max"] = max(seen["max"], len(b.list(f"runs/{tag}/warm/")))
            time.sleep(0.2)
    threading.Thread(target=watch, daemon=True).start()
    if kill:
        t_end = time.time() + 60            # at once, while that unit still runs
        while time.time() < t_end and not any(b.exists(f"runs/{tag}/units/{u['uid']}/ckpt.state")
                                              for u in b.units_by_status(tag, "leased")):
            time.sleep(0.05)
        os.killpg(w.pid, signal.SIGKILL)    # the supervisor; its solvers (own sessions) run out their budgets
        w.wait()
        t_end = time.time() + 30            # before the new worker's two start
        while time.time() < t_end and subprocess.run(["pgrep", "-c", "-x", "bnbdrv"], capture_output=True,
                                                     text=True).stdout.strip() not in ("", "0"):
            time.sleep(0.2)
        w = worker("w2")
    try:
        c.wait(timeout=600)
        w.wait(timeout=120)
    except subprocess.TimeoutExpired:
        for p in (c, w):
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGKILL)
        print(f"FAIL {tag}: the run did not finish (logs in {tmp})")
        return 1
    units = b.all_units(tag)
    paths = [os.path.join(tmp, f"{tag}.frontier.nodes")]
    if not b.get_file(f"runs/{tag}/frontier/nodes.txt", paths[0]):
        print(f"FAIL {tag}: no frontier/nodes.txt")
        return 1
    missing = 0
    for u in units:
        p = os.path.join(tmp, f"{tag}.{u['uid']}.nodes")
        if u.get("status") != "merged" or not b.get_file(f"runs/{tag}/units/{u['uid']}/a{u['attempt']}/nodes.txt", p):
            missing += 1
            continue
        paths.append(p)
    resumed = sum(bool(u.get("resumed_from_ckpt")) for u in units)
    if dnn:
        # the DNN states went through the store (ship.bits = 64: the tree is the
        # serial one only if each came back exact), and nothing is left of them
        # or of the units' state files once the run is over
        left = len(b.list(f"runs/{tag}/warm/")) + len([k for k in b.list(f"runs/{tag}/units/") if k.endswith(".state")])
        ok = seen["max"] > 0 and left == 0
        print(f"{'ok  ' if ok else 'FAIL'} {tag} parent-state store: up to {seen['max']} parent states during the run, "
              f"{left} files left after it")
        if not ok:
            return 1
    how = f"frontier + {len(paths) - 1} of {len(units)} units, {sum(u['attempt'] > 1 for u in units)} retried, " \
          f"{resumed} from a checkpoint"
    bad = compare(f"{tag} {inst}{' ' + extra if extra else ''}", read_trace(ref), read_trace(*paths), how)
    if missing:
        print(f"     {missing} unit(s) without a merged attempt's nodes.txt")
        bad = 1
    if kill and not resumed:
        print("     (no unit resumed from a checkpoint: the carry of ckpt.nodes.txt was not exercised)")
    # qapctl nodes assembles the same tree (what collect.sh and the tree viewer use)
    whole = os.path.join(tmp, f"{tag}.whole.nodes.gz")
    r = subprocess.run([sys.executable, "-m", "coord.qapctl", "nodes", tag, "--backend", be, "--out", whole],
                       cwd=ROOT, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL {tag}: qapctl nodes exited {r.returncode}: {r.stdout}{r.stderr}")
        return 1
    import gzip
    plain = os.path.join(tmp, f"{tag}.whole.nodes")
    with gzip.open(whole, "rb") as src, open(plain, "wb") as dst:
        dst.write(src.read())
    bad |= compare(f"{tag} qapctl nodes", read_trace(ref), read_trace(plain), "qapctl nodes --out FILE.gz")
    # the live tree, as published while the run went on, then rebuilt
    live = published_tree(b, tag)
    bad |= published_tree_problems(f"{tag} live tree", read_trace(ref), live)
    r = subprocess.run([sys.executable, "-m", "coord.qapctl", "publish", tag, "--backend", be],
                       cwd=ROOT, env=ENV, capture_output=True, text=True)
    rebuilt = published_tree(b, tag)
    strip = lambda t: (dict(t["meta"], generatedAt=None), t["files"])
    same = r.returncode == 0 and strip(live) == strip(rebuilt)
    print(f"{'ok  ' if same else 'FAIL'} {tag} tree rebuilt by qapctl publish: "
          f"{'the same meta and ' + str(len(live['files'])) + ' shards' if same else 'differs (or exited ' + str(r.returncode) + ')'}")
    return bad | (0 if same else 1)


def published_tree(b, tag):
    """a run's published tree, read as the viewer reads it: 0.json, then the
    shard of every boundary row"""
    key = lambda name: f"public/runs/{tag}/tree/{name}"
    meta = json.loads(b.get_bytes(key("meta.json")) or b"null")
    files, todo = {}, [0]
    while todo:
        r = todo.pop()
        data = b.get_bytes(key(f"{r}.json"))
        if data is None:
            files[r] = None
            continue
        files[r] = json.loads(data)
        todo += [row[0] for row in files[r]["rows"] if row[20] == 1 and row[0] not in files]
    return {"meta": meta, "files": files}


def published_tree_problems(what, ref, tree):
    meta, files = tree["meta"], tree["files"]
    problems = []
    if meta is None or files.get(0) is None:
        problems.append("no meta.json or 0.json")
    else:
        rows, count = {0: files[0]["rows"][0]}, collections.Counter()
        for r, shard in files.items():
            if shard is None:
                problems.append(f"no {r}.json for a boundary row")
                continue
            for row in shard["rows"][1:]:
                rows[row[0]] = row
                count[row[0]] += 1
        if any(c > 1 for c in count.values()):
            problems.append(f"{sum(c > 1 for c in count.values())} nodes in two shards")
        path = {}
        def fixes(i):
            if i not in path:
                row = rows[i]
                path[i] = "-" if row[1] < 0 else ",".join(x for x in (fixes(row[1]), row[3]) if x != "-")
            return path[i]
        got = {fixes(i): row for i, row in rows.items()}
        want = {r["fixes"]: r for r in ref}
        size, fat = collections.Counter(), collections.Counter()
        for r in ref:
            p = pairs(r["fixes"])
            for k in range(len(p) + 1):
                anc = ",".join(f"{a}:{b}" for a, b in p[:k]) or "-"
                size[anc] += 1
                fat[anc] += r["outcome"] in "QDCP"
        if "-" not in want:
            size["-"] += 1                   # the root, never traced: synthesized
        if set(got) != set(size):
            problems.append(f"{len(set(size) - set(got))} serial nodes missing, {len(set(got) - set(size))} extra")
        wrong = [f for f, row in got.items() if f in size and
                 ((f in want and row[4] != want[f]["outcome"]) or row[16] != size[f] or row[17] != fat[f])]
        if wrong:
            problems.append(f"{len(wrong)} nodes with another outcome, subtree size or fathom count, e.g. {wrong[0]}")
        if meta["nodes"] != len(size) or not meta["final"] or meta["orphans"] or meta["duplicates"]:
            problems.append(f"meta: {meta['nodes']} nodes, final {meta['final']}, {meta['orphans']} orphans, "
                            f"{meta['duplicates']} duplicates")
    ok = not problems
    print(f"{'ok  ' if ok else 'FAIL'} {what}: " +
          (f"{meta['nodes']} nodes in {len(files)} shards, read from 0.json through the boundaries" if ok
           else "; ".join(problems)))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="serial,killed,parallel,coord,coordkill,coordkilldnn")
    ap.add_argument("--keep", action="store_true", help="keep the temporary directory")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="qap-nodetrace-")
    bad = 0
    for part in a.only.split(","):
        t0 = time.time()
        bad += {"serial": serial, "killed": killed, "parallel": parallel, "coord": lambda t: coord(t, False),
                "coordkill": lambda t: coord(t, True),
                "coordkilldnn": lambda t: coord(t, True, dnn=True)}[part](tmp)
        print(f"     ({part}: {time.time() - t0:.0f} s)", flush=True)
    if bad or a.keep:
        print(f"files in {tmp}")
    else:
        shutil.rmtree(tmp)
    print("node trace: all checks passed" if not bad else f"node trace: {bad} check(s) FAILED")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
