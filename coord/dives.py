"""Dives runs: an estimate of the rest of a paused solve (docs/estimates_plan.md, step 2).

An estimate of a tree like tai60b's is about as costly as solving a fifth of it
(experiments/estimator_replay), so an estimate run is the start of a solve: the solve runs
for a budget of wall time and is paused (`qapctl pause`: its units return their open nodes),
then a dives run estimates the subtrees of those open nodes, and the solve can be resumed
from its store if the estimate says to go on.

    qapctl dives SRC --as NEW [--dive-b 0.4] [--budget S]   # one dives unit per open unit of SRC
    qapctl coordinator NEW / qapctl worker NEW               # as for any run
    qapctl estimate NEW [--out DIR]                          # SRC's work so far + the dives' estimate

A dives unit is one of SRC's queued units: its state file (read where SRC keeps it, never
copied or deleted: SRC can resume) and its parents' states (SRC's warm/). The worker runs the
estimator on it (bnbdrv -e -r, QAP_EST_SURV: survivor dives, latest/estimate.C) for the
unit's budget of wall seconds and keeps the estimator's JSON lines (dives.jsonl): one per
stratum (an open node, bounded once with its children: the exact part) and one per dive (an
unbiased estimate of the stratum's subtree beyond it). `estimate` adds up, for nodes and CPU,
by depth:

    SRC's merged totals  +  sum over strata of (exact part + mean of its dives)

with a bootstrap interval (dives resampled within their strata), the dives' deepest level
against the deepest level with nodes, and the strata no dive reached (they count their exact
part only, a known low bias where they have survivors).
"""
import json
import os
import random
import statistics
import tempfile
import time
import zlib

from .common import State, run_key, unit_key


def init_dives(b, src, new, dive_b=0.4, budget=None, log=print):
    """create run NEW: a dives unit for every queued unit of SRC (paused, nothing leased)"""
    s = b.get_run(src)
    if s is None:
        raise SystemExit(f"no run {src}")
    if b.get_run(new):
        raise SystemExit(f"run {new} exists")
    if s.get("config", {}).get("kind") == "dives":
        raise SystemExit(f"{src} is a dives run")
    counts = b.count_by_status(src)
    busy = {k: v for k, v in counts.items() if k in ("leased", "done", "returned", "failed") and v}
    if s.get("status") not in ("paused", "stopped") or busy:
        raise SystemExit(f"{src} is {s.get('status')} with {busy or 'nothing'} in flight: pause it "
                         f"(qapctl pause) and wait until its units have returned and merged")
    units = b.units_by_status(src, "queued")
    if not units:
        raise SystemExit(f"{src} has no queued units: nothing open to estimate")
    cfg = dict(s["config"])
    cfg.update({"kind": "dives", "source_run": src, "trace": False, "keep_unit_states": True,
                "dive": {"b": dive_b}})
    if budget is not None:
        cfg["budget"] = budget
    for name in ("instance.dat", "params.par", "incumbent.json"):
        data = b.get_bytes(run_key(src, name))
        if data is not None:
            b.put_bytes(run_key(new, name), data)
    meta = json.loads(b.get_bytes(run_key(src, "totals.json")) or b"{}")
    tk = meta.get("totals_key", run_key(src, "totals.state"))
    b.put_bytes(run_key(new, "totals.state"), b.get_bytes(tk))     # SRC's work so far, for the status
    b.put_bytes(run_key(new, "totals.json"), json.dumps({"last_seq": 0}).encode())
    b.put_bytes(run_key(new, "source_totals.json"), json.dumps({"src": src, "totals_key": tk,
                                                                "last_seq": meta.get("last_seq", 0)}).encode())
    inc = json.loads(b.get_bytes(run_key(src, "incumbent.json")) or b"{}")
    b.put_run(new, {"run_id": new, "status": "running", "created": time.time(), "incumbent": inc,
                    "visibility": "private", "config": cfg, "solver": s["solver"]})
    now = time.time()
    for u in units:
        doc = {k: u[k] for k in ("uid", "parent", "nodes", "depth_min", "depth_max", "by_depth", "min_bound",
                                 "lane", "in_key") if k in u}
        doc.update({"run": new, "status": "queued", "attempt": 0, "budget": cfg["budget"], "created": now,
                    "warm": u.get("warm", [])})
        b.put_unit(new, doc)
        b.send(doc["lane"], {"run": new, "uid": doc["uid"]})
        b.update_unit(new, doc["uid"], {"sent_at": time.time()})
    log(f"dives: run {new} from {src}: {len(units)} units, {sum(u['nodes'] for u in units)} open nodes, "
        f"{cfg['budget']:g} s each, weights exp({dive_b} x)")
    return new


def dive_seed(run_id, uid, attempt):
    return zlib.crc32(f"{run_id}/{uid}/{attempt}".encode())


def read_lines(text):
    strata, dives = {}, {}
    for line in text.splitlines():
        if not line.strip():
            continue
        x = json.loads(line)
        if x["type"] == "stratum":
            strata[x["key"]] = x
        elif x["type"] == "dive":
            dives.setdefault(x["key"], []).append(x)
    return strata, dives


def _by_depth(d):
    return {int(k): float(v) for k, v in (d or {}).items()}


def combine(src_state, strata, dives, reps=400, seed=1):
    """SRC's totals plus every stratum's exact part and the mean of its dives; a bootstrap over dives"""
    nl, tl = src_state.arrays.get("n_level", []), src_state.arrays.get("t_level", [])
    done_n = {d: float(v) for d, v in enumerate(nl) if v}
    done_c = {d: float(v) for d, v in enumerate(tl) if v}
    est_n, est_c = dict(), dict()
    for s in strata.values():
        for d, v in _by_depth(s.get("exact_by_depth")).items():
            est_n[d] = est_n.get(d, 0.0) + v
        for d, v in _by_depth(s.get("exact_cpu_by_depth")).items():
            est_c[d] = est_c.get(d, 0.0) + v
    for key, ds in dives.items():
        for x in ds:
            for d, v in _by_depth(x.get("by_depth")).items():
                est_n[d] = est_n.get(d, 0.0) + v / len(ds)
            for d, v in _by_depth(x.get("cpu_by_depth")).items():
                est_c[d] = est_c.get(d, 0.0) + v / len(ds)
    exact_n = sum(s["exact_nodes"] for s in strata.values())
    exact_c = sum(s["exact_cpu"] for s in strata.values())
    rest_n = exact_n + sum(statistics.fmean(x["nodes"] for x in ds) for ds in dives.values())
    rest_c = exact_c + sum(statistics.fmean(x["cpu"] for x in ds) for ds in dives.values())
    rng = random.Random(seed)
    boots_n, boots_c = [], []
    for _ in range(reps):
        bn, bc = exact_n, exact_c
        for ds in dives.values():
            pick = [ds[rng.randrange(len(ds))] for _ in ds]
            bn += statistics.fmean(x["nodes"] for x in pick)
            bc += statistics.fmean(x["cpu"] for x in pick)
        boots_n.append(bn)
        boots_c.append(bc)
    boots_n.sort()
    boots_c.sort()
    q = lambda xs, f: xs[min(int(f * len(xs)), len(xs) - 1)] if xs else 0.0
    no_dive = [k for k, s in strata.items() if s.get("survivors", 0) and k not in dives]
    deepest = max([x.get("deepest", 0) for ds in dives.values() for x in ds] + [0])
    spent = sum(s.get("setup_cpu", 0.0) for s in strata.values()) + \
        sum(x.get("cost", 0.0) for ds in dives.values() for x in ds)
    done_nodes, done_cpu = float(src_state.counters[0]), float(src_state.times[0])
    return {
        "done": {"nodes": done_nodes, "cpu": done_cpu},
        "rest": {"nodes": rest_n, "cpu": rest_c,
                 "nodes_5_95": [q(boots_n, .05), q(boots_n, .95)], "cpu_5_95": [q(boots_c, .05), q(boots_c, .95)]},
        "total": {"nodes": done_nodes + rest_n, "cpu": done_cpu + rest_c,
                  "nodes_5_95": [done_nodes + q(boots_n, .05), done_nodes + q(boots_n, .95)],
                  "cpu_5_95": [done_cpu + q(boots_c, .05), done_cpu + q(boots_c, .95)]},
        "strata": len(strata), "with_survivors": sum(1 for s in strata.values() if s.get("survivors", 0)),
        "dived": len(dives), "dives": sum(len(ds) for ds in dives.values()),
        "strata_with_survivors_not_dived": len(no_dive),
        "deepest_dive": deepest,
        "estimator_cpu": spent,
        "by_depth": {str(d): {"done_nodes": done_n.get(d, 0.0), "rest_nodes": est_n.get(d, 0.0),
                              "done_cpu": done_c.get(d, 0.0), "rest_cpu": est_c.get(d, 0.0)}
                     for d in sorted(set(done_n) | set(est_n) | set(done_c) | set(est_c))},
    }


def _cpu(s):
    return f"{s / 3600.0:,.1f} CPU h" if s >= 360 else f"{s:,.1f} CPU s"


def report(e):
    h = lambda s: s / 3600.0
    t, r, d = e["total"], e["rest"], e["done"]
    lines = [
        f"done so far: {d['nodes']:,.0f} nodes, {_cpu(d['cpu'])}",
        f"estimated rest: {r['nodes']:,.0f} nodes (5-95%: {r['nodes_5_95'][0]:,.0f} to {r['nodes_5_95'][1]:,.0f}), "
        f"{_cpu(r['cpu'])} ({_cpu(r['cpu_5_95'][0])} to {_cpu(r['cpu_5_95'][1])})",
        f"estimated tree: {t['nodes']:,.0f} nodes (5-95%: {t['nodes_5_95'][0]:,.0f} to {t['nodes_5_95'][1]:,.0f}), "
        f"{_cpu(t['cpu'])} ({_cpu(t['cpu_5_95'][0])} to {_cpu(t['cpu_5_95'][1])})",
        f"{e['strata']} open nodes ({e['with_survivors']} with children that branch), {e['dives']} dives in "
        f"{e['dived']} of them, deepest dive node at depth {e['deepest_dive']}; estimator {_cpu(e['estimator_cpu'])}",
    ]
    if e["strata_with_survivors_not_dived"]:
        lines.append(f"LOW: {e['strata_with_survivors_not_dived']} open nodes with branching children got no dive: "
                     f"their exact part only is counted")
    lines.append("")
    lines.append("| depth | nodes done | nodes estimated, rest | CPU h done | CPU h estimated, rest |")
    lines.append("|---|---|---|---|---|")
    for k, v in e["by_depth"].items():
        lines.append(f"| {k} | {v['done_nodes']:,.0f} | {v['rest_nodes']:,.0f} | {h(v['done_cpu']):,.2f} | "
                     f"{h(v['rest_cpu']):,.2f} |")
    deep_nodes = [int(k) for k, v in e["by_depth"].items() if v["done_nodes"] + v["rest_nodes"] > 0]
    if deep_nodes and e["deepest_dive"] and max(deep_nodes) > e["deepest_dive"]:
        lines.append("")
        lines.append(f"(nodes at depth {max(deep_nodes)} are in the solve's own part; the dives reached depth "
                     f"{e['deepest_dive']})")
    return "\n".join(lines)


def estimate_run(b, run_id, out=None, log=print):
    """`qapctl estimate RUN`: a dives run's estimate (its merged units' dives.jsonl) with its source's totals"""
    r = b.get_run(run_id)
    if r is None or r.get("config", {}).get("kind") != "dives":
        raise SystemExit(f"{run_id} is not a dives run")
    src = json.loads(b.get_bytes(run_key(run_id, "source_totals.json")))
    tmp = tempfile.mkdtemp(prefix=f"qapest-{run_id}-")
    p = os.path.join(tmp, "src.state")
    with open(p, "wb") as f:
        f.write(b.get_bytes(src["totals_key"]) or b.get_bytes(run_key(run_id, "totals.state")))
    st = State(p)
    strata, dives = {}, {}
    pending = 0
    for u in b.all_units(run_id):
        if u.get("status") != "merged":
            pending += 1
            continue
        data = b.get_bytes(unit_key(run_id, u["uid"], f"a{u['attempt']}", "dives.jsonl"))
        if data is None:
            continue
        s1, d1 = read_lines(data.decode())
        strata.update(s1)
        for k, v in d1.items():
            dives.setdefault(k, []).extend(v)
    e = combine(st, strata, dives)
    e["source_run"] = src["src"]
    e["units_pending"] = pending
    text = report(e)
    if pending:
        text = f"PART-WAY: {pending} dives units not merged yet\n" + text
    b.put_bytes(run_key(run_id, "estimate.json"), json.dumps(e, indent=1).encode())
    b.put_bytes(run_key(run_id, "estimate.txt"), (text + "\n").encode())
    if out:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "estimate.json"), "w") as f:
            json.dump(e, f, indent=1)
        with open(os.path.join(out, "estimate.txt"), "w") as f:
            f.write(text + "\n")
    log(text)
    return e


def estimate_flow(b, run_id, solve_wall, dive_budget, dive_b=0.4, poll=2.0, log=print):
    """`qapctl estimate-flow RUN`: the estimate run, in place of `qapctl coordinator RUN` on the
    coordinating host. RUN is solved until solve_wall seconds after its start, paused (dives_run set
    on it, for the workers that follow it: `qapctl worker --follow`), its units left to return and
    merge, then RUN-dives made from its open units and coordinated to the end (its results.txt is the
    estimate). Restartable at every step: each looks at what the store says is done. RUN stays
    paused and can be resumed (`qapctl resume RUN`) to go on with the solve."""
    from .coordinator import Coordinator, PENDING
    new = f"{run_id}-dives"
    r = b.get_run(run_id)
    if r is None:
        raise SystemExit(f"no run {run_id}")
    if r.get("status") == "finished":
        log(f"estimate-flow: {run_id} is finished: the tree closed, the solve is the answer")
        return b.get_bytes(run_key(run_id, "results.txt")).decode()
    if not r.get("dives_run"):
        if r.get("status") == "running":
            out = Coordinator(b, run_id, log=log).loop(poll=poll, max_wall=solve_wall)
            if out is not None:
                log(f"estimate-flow: {run_id} finished inside its {solve_wall:g} s: no estimate needed")
                return out
            if b.get_run(run_id).get("status") == "stopped":
                raise SystemExit(f"estimate-flow: {run_id} was stopped")
        b.update_run(run_id, {"status": "paused", "dives_run": new})
        log(f"estimate-flow: {run_id} paused after {solve_wall:g} s; its units return, then {new}")
    c = Coordinator(b, run_id, log=log)
    while True:
        counts = c.tick()
        if not any(counts.get(s, 0) for s in PENDING if s != "queued"):
            break
        time.sleep(poll)
    c._save_totals()
    if b.get_run(new) is None:
        init_dives(b, run_id, new, dive_b=dive_b, budget=dive_budget, log=log)
    out = Coordinator(b, new, log=log).loop(poll=poll)
    return out
