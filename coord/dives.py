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
import math
import os
import random
import re
import statistics
import tempfile
import time
import zlib

from .backend import check_queue_name
from .common import State, run_key, unit_key


DIVE_ENV = re.compile(r"^QAP_EST_[A-Z_]+$")

# The estimate's decision point (experiments/estimates_nov19, "The decision point"): an estimate is
# worth its CPU only while it is a small share of the tree it estimates. tai50b-est-261119 spent 37.5
# CPU h of estimator on a rest of 84 (truth): its setup alone bounded the top of a tree whose CPU is
# there. So init_dives marks a pilot, a random pilot_size(n) of the n dives units (the first of its
# shuffled order, so they are sent first), and once every pilot unit has merged the estimate-flow
# compares the estimator CPU the whole dives run will take (the pilot's mean per unit, times n) with
# the rest of the tree as the pilot estimates it (its open nodes standing, depth by depth, for every
# other unit's); above DECIDE (a share of that rest), the dives stop and the solve goes on, on the same
# workers (estimate_flow). The pilot is fixed in advance and waited for whole: the first units to merge
# are the small ones (a unit runs long when its open nodes branch), and on tho40-est2-261119 the first
# 55 read a rest of 42 CPU h against the 5,610 its dives found, which would have stopped them. Replayed
# on the collected runs (decide_sim.py, 400 random pilots each): tai50b-est decides to solve in every
# draw, after 3.1 estimator CPU h (median) of the 37.5 its dives took; tho40-est2 in none (its share
# 4.8-8.4%, 10-90%; the whole dives run 6.7%).
DECIDE = 0.15
DECIDE_MIN = 0.10


def pilot_size(n):
    """the decision point's pilot: a tenth of the n dives units, and 10 at least (half of a smaller run's)"""
    return max(math.ceil(DECIDE_MIN * n), min(10, math.ceil(n / 2)))


def decision(e, frac):
    """the decision point on estimate e (estimate_run's, with its "pilot" summary): verdict 'solve' (the
    dives run would cost more than frac of the rest the pilot estimates), 'dive' (they go on), 'wait' (a
    pilot unit not merged yet, or no pilot: a dives run from before Nov 19) or 'off' (frac 0 or None)"""
    p = e.get("pilot") or {}
    units = e.get("units", 0)
    projected = rest = share = None
    if p.get("units") and p["merged"] == p["units"]:
        projected = p["cpu"] / p["units"] * units
        rest = p["rest_cpu"]
        share = projected / rest if rest > 0 else None
    if not frac:
        verdict = "off"
    elif projected is None:
        verdict = "wait"
    elif share is None or share > frac:
        verdict = "solve"
    else:
        verdict = "dive"
    return {"frac": frac, "verdict": verdict, "pilot": p.get("units", 0), "pilot_merged": p.get("merged", 0),
            "units": units, "merged": (e.get("estimator") or {}).get("merged", 0),
            "estimator_cpu": (e.get("estimator") or {}).get("merged_cpu", 0.0),
            "projected_cpu": projected, "rest_cpu": rest, "share": share, "at": time.time()}


def final_cost(e, frac):
    """a dives run that ran to the end: what its estimator CPU was against the rest it found"""
    spent = (e.get("estimator") or {}).get("merged_cpu", 0.0)
    rest = e["rest"]["cpu"]
    return {"frac": frac, "verdict": "done", "projected_cpu": spent, "rest_cpu": rest,
            "share": spent / rest if rest > 0 else None, "units": e.get("units", 0), "at": time.time()}


def decision_text(d):
    if not d or d["verdict"] == "off":
        return ""
    sh = f"{d['share']:.0%}" if d.get("share") is not None else "all"
    if d["verdict"] == "wait":
        if not d.get("pilot"):
            return ""
        return (f"decision check: {d['pilot_merged']} of the {d['pilot']} pilot dives units merged; when all have, "
                f"the estimate says whether its dives are worth their CPU")
    if d["verdict"] == "done":
        return (f"the estimate's cost: {_cpu(d['projected_cpu'])} of estimator CPU, {sh} of the estimated rest of "
                f"the tree ({_cpu(d['rest_cpu'])})" +
                (f"; the estimate-flow allows {d['frac']:.0%}" if d.get("frac") else ""))
    cost = (f"the dives run would cost {_cpu(d['projected_cpu'])} of estimator CPU, {sh} of the rest of the tree "
            f"as its {d['pilot']} pilot units estimate it ({_cpu(d['rest_cpu'] or 0)})")
    if d["verdict"] == "solve":
        return (f"DECIDED TO SOLVE: {cost}, over the {d['frac']:.0%} an estimate may cost: the dives stopped "
                f"after {d['merged']} of {d['units']} units and the solve went on")
    return f"decision check: {cost}, within the {d['frac']:.0%} an estimate may cost: the dives go on"


def init_dives(b, src, new, dive_b="0.4", budget=None, log=print, sample=1.0, weights="unif", env=None,
               decide=None):
    """create run NEW: a dives unit for every queued unit of SRC (paused, nothing leased). env: more
    estimator settings for every unit (QAP_EST_* only), kept in the run's config; decide: the decision
    point's share (estimate_flow), kept there too, for the readouts"""
    env = dict(env or {})
    bad = [k for k in env if not DIVE_ENV.match(k)]
    if bad:
        raise SystemExit(f"dives: --dive-env takes QAP_EST_* settings only, not {', '.join(bad)}")
    s = b.get_run(src)
    if s is None:
        raise SystemExit(f"no run {src}")
    if b.get_run(new):
        raise SystemExit(f"run {new} exists")
    check_queue_name(new)
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
                "dive": {"b": str(dive_b), "sample": sample, "weights": weights, "env": env, "decide": decide}})
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
    b.create_queues(new)
    if not s.get("dives_run"):      # (estimate-flow sets it at the pause; a hand-made dives run, here)
        b.update_run(src, {"dives_run": new})
    b.put_run(new, {"run_id": new, "status": "running", "created": time.time(), "incumbent": inc,
                    "visibility": "private", "queues": "run", "config": cfg, "solver": s["solver"]})
    now = time.time()
    # in a random order, in one lane: every dives unit takes about the same wall time (its budget), and a
    # dives run read part-way then has a random part of SRC's units done, whose open nodes stand for
    # the others' depth by depth (estimate_run). In SRC's order (large lane first) the units left at the
    # end were all small and deep (tho40-est-dives, Nov 18)
    random.Random(zlib.crc32(new.encode())).shuffle(units)
    npilot = pilot_size(len(units))      # the decision point's pilot (decision()): the first of the shuffle
    for i, u in enumerate(units):
        doc = {k: u[k] for k in ("uid", "parent", "nodes", "depth_min", "depth_max", "by_depth", "min_bound",
                                 "in_key") if k in u}
        doc.update({"run": new, "status": "queued", "attempt": 0, "budget": cfg["budget"], "created": now,
                    "lane": "large", "warm": u.get("warm", []), "pilot": i < npilot})
        b.put_unit(new, doc)
        b.send(new, doc["lane"], {"run": new, "uid": doc["uid"]})
        b.update_unit(new, doc["uid"], {"sent_at": time.time()})
    log(f"dives: run {new} from {src}: {len(units)} units, {sum(u['nodes'] for u in units)} open nodes, "
        f"{cfg['budget']:g} s each, " + (f"Knuth's dives, g^{str(dive_b)[1:]}" if str(dive_b).startswith("k")
                                         else f"survivor dives, weights exp({dive_b} x)") +
        (f", {sample:g} of the open nodes sampled ({weights})" if sample < 1 else "") +
        "".join(f", {k}={v}" for k, v in sorted(env.items())) +
        (f"; the dives stop for the solve if they would cost over {decide:.0%} of the estimated rest" if decide else ""))
    return new


def dive_seed(run_id, uid, attempt):
    return zlib.crc32(f"{run_id}/{uid}/{attempt}".encode())


def read_lines(text, samples=None):
    strata, dives = {}, {}
    for line in text.splitlines():
        if not line.strip():
            continue
        x = json.loads(line)
        if x["type"] == "stratum":
            strata[x["key"]] = x
        elif x["type"] == "dive":
            dives.setdefault(x["key"], []).append(x)
        elif x["type"] == "sample" and samples is not None:
            samples.append(x)
    return strata, dives


def _by_depth(d):
    return {int(k): float(v) for k, v in (d or {}).items()}


def combine(src_state, strata, dives, reps=400, seed=1):
    """SRC's totals plus every stratum's exact part and the mean of its dives, each weighted 1/pi (the
    stratum's probability of being sampled: 1 without QAP_EST_SAMPLE); a bootstrap over the dives within
    each stratum and, where strata were sampled (pi < 1), over those strata"""
    nl, tl = src_state.arrays.get("n_level", []), src_state.arrays.get("t_level", [])
    done_n = {d: float(v) for d, v in enumerate(nl) if v}
    done_c = {d: float(v) for d, v in enumerate(tl) if v}
    pi = {k: float(s.get("pi", 1.0)) or 1.0 for k, s in strata.items()}
    est_n, est_c = dict(), dict()
    for k, s in strata.items():
        for d, v in _by_depth(s.get("exact_by_depth")).items():
            est_n[d] = est_n.get(d, 0.0) + v / pi[k]
        for d, v in _by_depth(s.get("exact_cpu_by_depth")).items():
            est_c[d] = est_c.get(d, 0.0) + v / pi[k]
    for key, ds in dives.items():
        w = pi.get(key, 1.0)
        for x in ds:
            for d, v in _by_depth(x.get("by_depth")).items():
                est_n[d] = est_n.get(d, 0.0) + v / len(ds) / w
            for d, v in _by_depth(x.get("cpu_by_depth")).items():
                est_c[d] = est_c.get(d, 0.0) + v / len(ds) / w
    def part(k, pick=None):
        s = strata[k]
        ds = dives.get(k, [])
        if pick is not None and ds:
            ds = pick(ds)
        n = s["exact_nodes"] + (statistics.fmean(x["nodes"] for x in ds) if ds else 0.0)
        c = s["exact_cpu"] + (statistics.fmean(x["cpu"] for x in ds) if ds else 0.0)
        return n / pi[k], c / pi[k]
    rest_n = rest_c = 0.0
    for k in strata:
        n, c = part(k)
        rest_n += n
        rest_c += c
    rng = random.Random(seed)
    fixed = [k for k in strata if pi[k] >= 1.0]
    sampled = [k for k in strata if pi[k] < 1.0]
    boots_n, boots_c = [], []
    redraw = lambda ds: [ds[rng.randrange(len(ds))] for _ in ds]
    for _ in range(reps):
        bn = bc = 0.0
        for k in fixed + ([sampled[rng.randrange(len(sampled))] for _ in sampled] if sampled else []):
            n, c = part(k, redraw)
            bn += n
            bc += c
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
    if e.get("open_nodes") and e.get("open_sampled", 0) < e["open_nodes"]:
        lines.append(f"{e['open_sampled']:,} of the {e['open_nodes']:,} open nodes sampled (each weighted 1/pi)")
    if e.get("jackknife_cpu") is not None:
        lines.append(f"jackknife over units: standard error {e['jackknife_nodes']:,.0f} nodes, {_cpu(e['jackknife_cpu'])}")
    g = e.get("diag")
    if g:
        lines.append(f"the largest dive makes {g['top1_cpu']:.0%} of the rest's CPU, the 10 largest {g['top10_cpu']:.0%}; "
                     f"{g['one_dive_strata']} of {g['dived_strata']} dived open nodes got one dive; estimator CPU: "
                     f"setup {_cpu(g['setup_cpu'])}, dives {_cpu(g['dive_cpu'])}")
        if g.get("deepest_well_sampled") is not None:
            lines.append(f"dives going deeper, by depth (reached, went on): " + ", ".join(
                f"{d}: {v[0]}/{v[1]}" for d, v in g["continuation"].items()))
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


def _parts(strata, dives):
    """each stratum's estimated subtree (its exact part and its dives' mean, weighted 1/pi): nodes, CPU"""
    out = {}
    for k, st in strata.items():
        pi = float(st.get("pi", 1.0)) or 1.0
        ds = dives.get(k, [])
        n = st["exact_nodes"] + (statistics.fmean(x["nodes"] for x in ds) if ds else 0.0)
        c = st["exact_cpu"] + (statistics.fmean(x["cpu"] for x in ds) if ds else 0.0)
        out[k] = (n / pi, c / pi)
    return out


def _open_by_depth(u):
    return {int(d): int(v) for d, v in (u.get("by_depth") or {}).items()}


def unit_rest(units, keys, strata, parts):
    """per merged unit: {open depth: [estimated nodes, CPU, open nodes]} of its strata"""
    out = []
    for u in units:
        acc = {}
        for d, n in _open_by_depth(u).items():
            acc.setdefault(d, [0.0, 0.0, 0])[2] += n
        for k in keys.get(u["uid"], []):
            if k in parts:
                d = int(strata[k].get("depth", 0))
                a = acc.setdefault(d, [0.0, 0.0, 0])
                a[0] += parts[k][0]
                a[1] += parts[k][1]
        out.append(acc)
    return out


def rest_from(measured, others):
    """the rest of the tree from the measured units' estimates, the other units' open nodes standing for
    as much, depth by depth, as the measured units' open nodes at that depth (a ratio estimator): its
    (nodes, CPU, imputed nodes, imputed CPU, open nodes at depths no measured unit has)"""
    tot = {}
    for acc in measured:
        for d, (n, c, o) in acc.items():
            t = tot.setdefault(d, [0.0, 0.0, 0])
            t[0] += n; t[1] += c; t[2] += o
    rn = sum(t[0] for t in tot.values())
    rc = sum(t[1] for t in tot.values())
    imp_n = imp_c = 0.0
    unmatched = {}
    for d, o in others.items():
        t = tot.get(d)
        if not t or not t[2]:
            if o:
                unmatched[d] = unmatched.get(d, 0) + o
            continue
        imp_n += o * t[0] / t[2]
        imp_c += o * t[1] / t[2]
    return rn + imp_n, rc + imp_c, imp_n, imp_c, unmatched


def jackknife(accs, pending_open, groups=20, seed=1):
    """delete-a-group jackknife over the merged units (each group's open nodes imputed from the rest
    when it is left out): the standard errors of the rest's nodes and CPU"""
    if len(accs) < 4:
        return None, None
    idx = list(range(len(accs)))
    random.Random(seed).shuffle(idx)
    g = min(groups, len(accs))
    thetas = []
    for j in range(g):
        out = set(idx[j::g])
        others = dict(pending_open)
        for i in out:
            for d, (_, _, o) in accs[i].items():
                others[d] = others.get(d, 0) + o
        n, c, *_ = rest_from([accs[i] for i in idx if i not in out], others)
        thetas.append((n, c))
    mn = statistics.fmean(t[0] for t in thetas)
    mc = statistics.fmean(t[1] for t in thetas)
    f = (g - 1) / g
    return (f * sum((t[0] - mn) ** 2 for t in thetas)) ** 0.5, (f * sum((t[1] - mc) ** 2 for t in thetas)) ** 0.5


def diagnose(strata, dives, e):
    """what says how far to trust an estimate (experiments/tho40_est_diag): the largest dives' share,
    the dived strata with one dive, how dives continue by depth, the estimator's CPU in setup and dives"""
    contrib = []
    for k, ds in dives.items():
        pi = float(strata.get(k, {}).get("pi", 1.0)) or 1.0
        for x in ds:
            contrib.append((x["nodes"] / len(ds) / pi, x["cpu"] / len(ds) / pi))
    rn, rc = e["rest"]["nodes"] or 1.0, e["rest"]["cpu"] or 1.0
    top = lambda i, k: sum(sorted((t[i] for t in contrib), reverse=True)[:k])
    reach, cont = {}, {}
    for ds in dives.values():
        for x in ds:
            deepest = int(x.get("deepest", 0))
            for d in {int(t) for t in (x.get("by_depth") or {})}:
                reach[d] = reach.get(d, 0) + 1
                if deepest > d:
                    cont[d] = cont.get(d, 0) + 1
    setup = sum(st.get("setup_cpu", 0.0) for st in strata.values())
    dive_cpu = sum(x.get("cost", 0.0) for ds in dives.values() for x in ds)
    by = e.get("by_depth", {})
    mass = max(by, key=lambda d: by[d]["rest_nodes"] + by[d]["done_nodes"]) if by else None
    deep = [d for d in sorted(reach) if reach[d] >= 20]
    return {
        "top1_nodes": top(0, 1) / rn, "top10_nodes": top(0, 10) / rn,
        "top1_cpu": top(1, 1) / rc, "top10_cpu": top(1, 10) / rc,
        "dived_strata": len(dives), "one_dive_strata": sum(1 for ds in dives.values() if len(ds) == 1),
        "continuation": {str(d): [reach[d], cont.get(d, 0)] for d in sorted(reach)},
        "deepest_well_sampled": deep[-1] if deep else None,
        "continuation_there": (cont.get(deep[-1], 0) / reach[deep[-1]]) if deep else None,
        "mass_depth": int(mass) if mass is not None else None,
        "setup_cpu": setup, "dive_cpu": dive_cpu,
    }


def flags(e):
    """the warnings a readout leads with: each a short line, empty when none applies"""
    g = e.get("diag", {})
    out = []
    if e.get("units_pending"):
        sh = e["imputed"]["cpu"] / e["rest"]["cpu"] if e["rest"]["cpu"] else 0.0
        out.append(f"PART-WAY: {e['units_pending']} of {e['units']} dives units not merged; their open nodes are "
                   f"imputed, depth by depth, from the merged units' ({sh:.0%} of the rest's CPU)")
    if e.get("slow_pending"):
        out.append(f"SLOW PENDING: {e['slow_pending']} pending units have run longer than 90% of the merged ones took "
                   f"({e['merged_wall_p90']:.0f} s): a unit runs long when its open nodes branch, so these likely hold "
                   f"the large subtrees, and the imputation stands in for them with small ones: read the estimate as "
                   f"a floor until they merge")
    if e.get("imputed", {}).get("unmatched_open"):
        by = e["imputed"].get("unmatched_by_depth", {})
        out.append(f"GAP: {e['imputed']['unmatched_open']} open nodes of pending units are at depths no merged "
                   f"unit has open nodes at yet ({', '.join(f'depth {d}: {o}' for d, o in by.items())}): not "
                   f"counted, so the estimate is low by their share until one is merged")
    if max(g.get("top10_cpu", 0), g.get("top10_nodes", 0)) > 0.5:
        out.append(f"HEAVY TAIL: the 10 largest dives make {g['top10_nodes']:.0%} of the rest's nodes and "
                   f"{g['top10_cpu']:.0%} of its CPU (the largest {g['top1_nodes']:.0%}, {g['top1_cpu']:.0%}): "
                   f"more dives where they come from")
    if g.get("dived_strata") and g["one_dive_strata"] > 0.5 * g["dived_strata"]:
        out.append(f"FEW DIVES: {g['one_dive_strata']} of {g['dived_strata']} dived open nodes got one dive "
                   f"(QAP_EST_MIN_DIVES, a longer dive budget)")
    if g.get("setup_cpu", 0) > g.get("dive_cpu", 0) and g.get("setup_cpu", 0) > 0:
        out.append(f"SETUP-BOUND: bounding the sampled open nodes took {_cpu(g['setup_cpu'])}, the dives "
                   f"{_cpu(g['dive_cpu'])} (sample fewer, or by depth: QAP_EST_SAMPLE_DEPTHS)")
    if g.get("continuation_there") is not None and g["continuation_there"] > 0.05:
        out.append(f"DEEP TAIL: {g['continuation_there']:.0%} of the dives at depth {g['deepest_well_sampled']} "
                   f"went deeper, where few dives reached: the deep levels are thinly sampled")
    if e.get("jackknife_cpu") and e["rest"]["cpu"] and e["jackknife_cpu"] > 0.25 * e["rest"]["cpu"]:
        out.append(f"WIDE: the jackknife over units gives a standard error of "
                   f"{e['jackknife_cpu'] / e['rest']['cpu']:.0%} of the rest's CPU")
    return out


def estimate_run(b, run_id, out=None, log=print, save=True, cache=None, decide=None):
    """`qapctl estimate RUN`: a dives run's estimate (its merged units' dives.jsonl) with its source's
    totals; part-way, the pending units' open nodes imputed depth by depth from the merged units'.
    cache: a dict kept between calls (the estimate-flow's decision checks), so that each merged unit's
    dives.jsonl is read once. decide: the decision point's share (default: the run's, from estimate-flow)"""
    r = b.get_run(run_id)
    if r is None or r.get("config", {}).get("kind") != "dives":
        raise SystemExit(f"{run_id} is not a dives run")
    src = json.loads(b.get_bytes(run_key(run_id, "source_totals.json")))
    tmp = tempfile.mkdtemp(prefix=f"qapest-{run_id}-")
    p = os.path.join(tmp, "src.state")
    sk = ("src", src["totals_key"])
    if cache is not None and sk in cache:
        data = cache[sk]
    else:
        data = b.get_bytes(src["totals_key"]) or b.get_bytes(run_key(run_id, "totals.state"))
        if cache is not None:
            cache[sk] = data
    with open(p, "wb") as f:
        f.write(data)
    st = State(p)
    strata, dives = {}, {}
    samples = []
    merged, keys, pending_open = [], {}, {}
    units = b.all_units(run_id)
    # what the estimator has spent: the merged units' CPU (the decision point; its pilot below)
    spend = {"merged": sum(1 for u in units if u.get("status") == "merged"),
             "merged_cpu": sum(float(u.get("cpu") or 0.0) for u in units if u.get("status") == "merged")}
    for u in units:
        data = None
        if u.get("status") == "merged":
            ck = (u["uid"], u["attempt"])
            if cache is not None and ck in cache:
                data = cache[ck]
            else:
                data = b.get_bytes(unit_key(run_id, u["uid"], f"a{u['attempt']}", "dives.jsonl"))
                if cache is not None and data is not None:
                    cache[ck] = data
        if data is None:
            for d, n in _open_by_depth(u).items():
                pending_open[d] = pending_open.get(d, 0) + n
            continue
        s1, d1 = read_lines(data.decode(), samples)
        strata.update(s1)
        keys[u["uid"]] = list(s1)
        merged.append(u)
        for k, v in d1.items():
            dives.setdefault(k, []).extend(v)
    e = combine(st, strata, dives)
    e["source_run"] = src["src"]
    if samples:
        e["open_nodes"] = sum(x.get("open", 0) for x in samples)
        e["open_sampled"] = sum(x.get("kept", 0) for x in samples)
    pending = len(units) - len(merged)
    e["units"], e["units_pending"] = len(units), pending
    # the pending units still running are not a random part: a unit runs longer when its open nodes
    # branch (setup bounds every child, then dives), so the last to merge hold the large subtrees
    # (tai50b-est-261119 at 06:05: its 7 pending units were the 7 depth-4 nodes that branch, under which
    # the tree goes on to depth 14; the imputation from the 252 merged, nearly all fathomed, read 0.28)
    walls = sorted(float(u.get("wall") or 0) for u in merged if u.get("wall"))
    p90 = walls[min(int(0.9 * len(walls)), len(walls) - 1)] if walls else None
    now = time.time()
    running = [now - float(u["leased_at"]) for u in units
               if u.get("status") in ("leased", "done", "returned") and u.get("leased_at")]
    e["slow_pending"] = sum(1 for t in running if p90 is not None and t > p90)
    e["merged_wall_p90"] = p90
    accs = unit_rest(merged, keys, strata, _parts(strata, dives))
    rn, rc, imp_n, imp_c, unmatched = rest_from(accs, pending_open)
    e["imputed"] = {"nodes": imp_n, "cpu": imp_c, "open": sum(pending_open.values()),
                    "unmatched_open": sum(unmatched.values()),
                    "unmatched_by_depth": {str(d): o for d, o in sorted(unmatched.items())}}
    if pending:
        # the merged part as combine() has it, plus the pending units' share; the bootstrap's band moves with it
        for part in ("rest", "total"):
            e[part]["nodes"] += imp_n
            e[part]["cpu"] += imp_c
            e[part]["nodes_5_95"] = [v + imp_n for v in e[part]["nodes_5_95"]]
            e[part]["cpu_5_95"] = [v + imp_c for v in e[part]["cpu_5_95"]]
    e["jackknife_nodes"], e["jackknife_cpu"] = jackknife(accs, pending_open)
    e["diag"] = diagnose(strata, dives, e)
    e["flags"] = flags(e)
    e["estimator"] = spend
    # the pilot (decision()): its units' CPU, and once all have merged the rest as they alone estimate
    # it, every other unit's open nodes imputed from theirs depth by depth
    pilot = [u for u in units if u.get("pilot")]
    if pilot:
        pm = [u for u in merged if u.get("pilot")]
        e["pilot"] = {"units": len(pilot), "merged": len(pm), "cpu": sum(float(u.get("cpu") or 0.0) for u in pm),
                      "rest_cpu": None, "rest_nodes": None}
        if len(pm) == len(pilot):
            others = {}
            for u in units:
                if not u.get("pilot"):
                    for d, n in _open_by_depth(u).items():
                        others[d] = others.get(d, 0) + n
            pn, pc, *_ = rest_from(unit_rest(pm, keys, strata, _parts(strata, dives)), others)
            e["pilot"].update(rest_cpu=pc, rest_nodes=pn)
    # the decision point: the one taken (a dives run the estimate-flow stopped for the solve), else the
    # check as it stands now, with the run's share (or the one given)
    frac = decide if decide is not None else r.get("config", {}).get("dive", {}).get("decide")
    if r.get("decision"):
        e["decision"] = r["decision"]
    elif not pending:       # the dives ran to the end: what the estimate cost, against what it found
        e["decision"] = final_cost(e, frac)
    else:
        e["decision"] = decision(e, frac) if frac else None
    text = report(e)
    if e["flags"]:
        text = "\n".join(e["flags"]) + "\n" + text
    dt = decision_text(e["decision"])
    if dt:
        text = dt + "\n" + text
    e["text"] = text
    if save:      # (a readout made with a read-only role writes nothing: coord/readout.py)
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


def export_dives(b, run_id, out_dir, log=print):
    """`qapctl dives-export RUN --out DIR`: a dives run's record for analysis, as the user gathered
    tho40-est-dives's by hand (infra/phase1/runs/tho40-est-dives/README.md): dives.jsonl.gz, every merged
    unit's dives.jsonl joined in merge order (the lines coord/dives.py reads), and units.tsv, a row per
    unit, merged or not, with its line range in dives.jsonl.gz. Names no bucket."""
    import gzip
    r = b.get_run(run_id)
    if r is None or r.get("config", {}).get("kind") != "dives":
        raise SystemExit(f"{run_id} is not a dives run")
    os.makedirs(out_dir, exist_ok=True)
    units = sorted(b.all_units(run_id), key=lambda u: (u.get("merged_seq") is None, u.get("merged_seq") or 0, u["uid"]))
    rows, n = [], 0
    tmp = os.path.join(out_dir, "dives.jsonl.gz.tmp")
    with gzip.open(tmp, "wt") as f:
        for u in units:
            data = None
            if u.get("status") == "merged":
                data = b.get_bytes(unit_key(run_id, u["uid"], f"a{u['attempt']}", "dives.jsonl"))
            lines = [l for l in (data.decode().splitlines() if data else []) if l.strip()]
            for l in lines:
                f.write(l + "\n")
            rows.append([u["uid"], u.get("status", ""), u.get("attempt", ""), u.get("lane", ""), u.get("depth_min", ""),
                         u.get("depth_max", ""), u.get("nodes", ""), u.get("merged_seq", ""), u.get("cpu", ""),
                         u.get("wall", ""), n if lines else "", n + len(lines) if lines else "",
                         "yes" if data is not None else "no"])
            n += len(lines)
    os.replace(tmp, os.path.join(out_dir, "dives.jsonl.gz"))
    with open(os.path.join(out_dir, "units.tsv"), "w") as f:
        f.write("uid\tstatus\tattempt\tlane\tdepth_min\tdepth_max\tnodes\tmerged_seq\tcpu_s\twall_s\tline_from\tline_to\thas_dives\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    log(f"dives-export: {run_id}: {n} lines from {sum(1 for r in rows if r[-1] == 'yes')} of {len(rows)} units in {out_dir}")
    return 0


def _decide_solve(b, run_id, new, d, e=None, log=print):
    """the decision point said solve: RUN resumed (first, so that a worker leaving NEW finds it running),
    NEW stopped with the decision, its running units cancelled (their workers kill them at the next
    heartbeat: the lease is gone) and its queues deleted, the readout written. Idempotent: a restarted
    estimate-flow finishes what an interrupted one began."""
    from .coordinator import drop_queues, ensure_queues
    if e is not None:       # the estimate as it stood when the dives stopped, for the readouts
        b.put_bytes(run_key(new, "estimate.json"), json.dumps(e, indent=1).encode())
        b.put_bytes(run_key(new, "estimate.txt"), (e["text"] + "\n").encode())
    src = b.get_run(run_id)
    if src.get("status") == "paused":
        b.update_run(run_id, {"status": "running", "estimate_decision": d})
        n = ensure_queues(b, run_id, log=log, resend=True)
        log(f"estimate-flow: {run_id} resumed for the solve ({n} queued units sent again)")
    elif not src.get("estimate_decision"):
        b.update_run(run_id, {"estimate_decision": d})
    nr = b.get_run(new)
    if nr is not None and nr.get("status") != "stopped":
        b.update_run(new, {"status": "stopped", "stopped_at": time.time(), "decision": d})
    if nr is not None:
        for u in b.all_units(new):
            if u.get("status") == "leased":
                b.update_unit(new, u["uid"], {"status": "cancelled"},
                              expect={"status": "leased", "attempt": u["attempt"]})
        if nr.get("queues"):
            drop_queues(b, new)
    from .readout import write_readout
    write_readout(b, run_id, log=log)


def estimate_flow(b, run_id, solve_wall, dive_budget, dive_b="0.4", poll=2.0, log=print, sample=1.0, weights="unif",
                  env=None, decide=DECIDE, check_every=120.0):
    """`qapctl estimate-flow RUN`: the estimate run, in place of `qapctl coordinator RUN` on the
    coordinating host. RUN is solved until solve_wall seconds after its start, paused (dives_run set
    on it, for the workers that follow it: `qapctl worker --follow`), its units left to return and
    merge, then RUN-dives made from its open units and coordinated to the end (its results.txt is the
    estimate). Restartable at every step: each looks at what the store says is done. RUN stays
    paused and can be resumed (`qapctl resume RUN`) to go on with the solve.

    The decision point (decide, a share; 0 or None: off): while the dives run, at most every
    check_every seconds that a unit has merged, the estimate so far; once every pilot unit has merged,
    if the whole dives run would cost over decide of the rest the pilot estimates (decision()), the
    dives stop and RUN is resumed and coordinated here to the end, its workers coming back from the
    dives run (`worker --follow`); else the dives go on to the end."""
    from .coordinator import Coordinator, PENDING
    new = f"{run_id}-dives"
    r = b.get_run(run_id)
    if r is None:
        raise SystemExit(f"no run {run_id}")
    if r.get("status") == "finished":
        log(f"estimate-flow: {run_id} is finished: the tree closed, the solve is the answer")
        return b.get_bytes(run_key(run_id, "results.txt")).decode()
    if r.get("estimate_decision"):
        # decided before (a restart): finish the switch, then coordinate the solve
        _decide_solve(b, run_id, new, r["estimate_decision"], log=log)
        log(f"estimate-flow: {run_id}'s estimate decided to solve: coordinating the solve")
        return Coordinator(b, run_id, log=log).loop(poll=poll)
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
    if b.get_run(new) is None and not b.units_by_status(run_id, "queued"):
        # the tree closed while its units returned: nothing is open, the solve is the answer
        log(f"estimate-flow: {run_id} has nothing open after the pause: the tree is closed")
        return c.finish()
    if b.get_bytes(run_key(run_id, "results_partial.txt")) is None:
        # the paused solve's RESULTS block, what the estimate's "done so far" is: in git with the estimate
        # (collect.sh), as a finished solve's results.txt is. First its status as the pause left it: the
        # last one was written before the returning units merged (nug17 in dives_test: 20,068 nodes against
        # the 35,890 merged), and the readout reads it
        try:
            c.status(b.count_by_status(run_id), force=True)
            c.write_partial(f"the solve paused after {solve_wall:g} s of wall for the estimate {new}")
        except Exception as ex:       # a report must never stop the estimate
            log(f"estimate-flow: results_partial.txt failed: {ex!r}")
    if b.get_run(new) is None:
        init_dives(b, run_id, new, dive_b=dive_b, budget=dive_budget, log=log, sample=sample, weights=weights, env=env,
                   decide=decide)
    frac = (b.get_run(new) or {}).get("config", {}).get("dive", {}).get("decide")   # the run's, on a restart
    cache, seen = {}, {"merged": -1, "t": 0.0, "decided": False}

    def check(counts):
        # once: when every pilot unit has merged (decision()); the dives go on after a 'dive'
        if not frac or seen["decided"]:
            return False
        n = counts.get("merged", 0)
        if n == seen["merged"] or time.time() - seen["t"] < check_every:
            return False
        seen.update(merged=n, t=time.time())
        try:
            e = estimate_run(b, new, log=lambda *_: None, save=False, cache=cache)
        except Exception as ex:      # a check must never stop the dives
            log(f"estimate-flow: the decision check failed: {ex!r}")
            return False
        d = e.get("decision") or {}
        if d.get("verdict") in ("dive", "wait"):
            log(f"estimate-flow: {decision_text(d) or 'decision check: no pilot units (a dives run from before Nov 19)'}")
        if d.get("verdict") == "dive" or not e.get("pilot"):
            seen["decided"] = True
        if d.get("verdict") != "solve":
            return False
        log(f"estimate-flow: {decision_text(d)}")
        _decide_solve(b, run_id, new, d, e, log=log)     # (e's text leads with the decision)
        return True

    out = Coordinator(b, new, log=log).loop(poll=poll, check=check)
    if out is None and (b.get_run(run_id) or {}).get("estimate_decision"):
        log(f"estimate-flow: the dives stopped; coordinating the solve {run_id} to the end")
        return Coordinator(b, run_id, log=log).loop(poll=poll)
    return out
