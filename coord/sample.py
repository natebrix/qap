"""Sampled runs (Nov 19; experiments/estimates_nov19, "Tuning on paired subtrees"): a run whose first
units are a sample of another run's open nodes, so that two or more runs, the arms of an A/B with
their own params, solve the same nodes, each to the end, side by side, and are compared node by node.

The source SRC is a paused run with queued units (an estimate run: tho40-est2-261119), and the size
measure is a dives run of it (default SRC's own, SRC-dives): every open node the dives run sampled (a
stratum, with its probability pi of being sampled) has an estimate of its subtree's CPU, its exact
part plus its dives' mean. tho40's subtrees are uneven (at depth 3 a median of 0.03 CPU h, a mean
of 0.40, the largest 1% holding 35%), so a uniform sample's ratio between arms would ride on its few
largest subtrees; the draw favours the large ones instead, in proportion to a power of that estimate
(power 0.5 by default, qapctl --sample-power: for a ratio between arms whose effect is about proportional to a subtree's
size, the variance at a fixed cost is least with pi ~ sqrt(size), the cost of a node being its size;
in proportion to the size itself the sample holds every large subtree and costs several times as
much: tho40-est2's depth-3 strata, 50 nodes, about 300 CPU h an arm against 113, for a ratio known to
2.3% against 3.8% at 20% of node-to-node variation; experiments/estimates_nov19):

  at depth d, with n_d wanted out of the N_d strata there, s_i the estimate to the power:
  q_i = (1 - mix) s_i / sum s + mix / N_d,
  pi_i = min(1, n_d q_i) (the excess of the capped ones spread over the rest), each stratum taken
  with probability pi_i (Poisson sampling, seeded, in the order of the strata's keys)

and every node taken stands for w_i = 1 / (pi_dives_i x pi_i) of SRC's open nodes at its depth. The
same SRC, dives run, plan, mix and seed give the same sample on every arm; each run keeps it as
runs/RUN/sample.json (experiments/estimates_nov19/paired.py compares the arms with it).

The frontier: the sampled nodes' lines and their parents' warm starts, out of SRC's unit state files
(read, never changed: SRC can still be resumed), the parents' DNN states fetched from SRC's store,
the statistics zeroed, and the settings record the run's own (an arm's params differ from SRC's;
the record's hash is what the solver checks on -r)."""
import json
import os
import random
import subprocess

from . import warm
from .common import State, run_key, unit_key


def node_key(text, n):
    """a node line's fixes in the order they were made, as the estimator's strata and the node trace
    name a node: "facility:location,..." (latest/estimate.C, fixes_of)"""
    f = text.split()
    nfix = int(f[4])
    p, hist = f[5:5 + n], f[5 + n:5 + 2 * n]
    return ",".join(f"{hist[k]}:{p[int(hist[k])]}" for k in range(nfix))


def strata_sizes(b, dv):
    """{key: {"depth", "pi", "cpu", "nodes", "uid"}} for every stratum of dives run DV's merged units"""
    out = {}
    for u in b.all_units(dv):
        if u.get("status") != "merged":
            continue
        data = b.get_bytes(unit_key(dv, u["uid"], f"a{u['attempt']}", "dives.jsonl"))
        if not data:
            continue
        strata, dives = {}, {}
        for line in data.decode().splitlines():
            if not line.strip():
                continue
            x = json.loads(line)
            if x["type"] == "stratum":
                strata[x["key"]] = x
            elif x["type"] == "dive":
                dives.setdefault(x["key"], []).append(x)
        for k, s in strata.items():
            ds = dives.get(k, [])
            cpu = s["exact_cpu"] + (sum(x["cpu"] for x in ds) / len(ds) if ds else 0.0)
            nodes = s["exact_nodes"] + (sum(x["nodes"] for x in ds) / len(ds) if ds else 0.0)
            out[k] = {"depth": int(s["depth"]), "pi": float(s.get("pi", 1.0)) or 1.0, "cpu": cpu, "nodes": nodes,
                      "uid": u["uid"]}
    return out


def parse_plan(plan):
    """"3:150,4:150" -> {3: 150, 4: 150}"""
    out = {}
    for item in str(plan).split(","):
        d, sep, n = item.partition(":")
        if not sep or not d.strip().isdigit() or not n.strip().isdigit() or int(n) < 1:
            raise SystemExit(f"sample: {item!r} is not DEPTH:COUNT (e.g. 3:150,4:150)")
        out[int(d)] = int(n)
    return out


def inclusion(sizes, n, mix):
    """the inclusion probabilities of a size-proportional sample of n out of len(sizes)"""
    N = len(sizes)
    if N == 0:
        return []
    tot = sum(sizes)
    q = [(1.0 - mix) * (s / tot if tot > 0 else 1.0 / N) + mix / N for s in sizes]
    pi = [0.0] * N
    capped = [False] * N
    for _ in range(50):          # the capped ones (pi = 1) taken out, the rest scaled to the remainder
        m = n - sum(1 for c in capped if c)
        free = sum(q[i] for i in range(N) if not capped[i])
        more = False
        for i in range(N):
            if capped[i]:
                pi[i] = 1.0
                continue
            pi[i] = m * q[i] / free if free > 0 and m > 0 else 0.0
            if pi[i] >= 1.0:
                capped[i] = True
                more = True
        if not more:
            break
    return [min(1.0, x) for x in pi]


def draw(strata, plan, seed=1, mix=0.2, power=0.5):
    """the sample: a list of {key, depth, uid, pi_dives, pi, w, est_cpu, est_nodes}, in key order"""
    picked = []
    for d in sorted(plan):
        keys = sorted(k for k, s in strata.items() if s["depth"] == d)
        if not keys:
            raise SystemExit(f"sample: the dives run has no open node at depth {d}")
        pis = inclusion([max(strata[k]["cpu"], 0.0) ** power for k in keys], min(plan[d], len(keys)), mix)
        rng = random.Random(f"{seed}:{d}")
        for k, pi in zip(keys, pis):
            if rng.random() < pi:
                s = strata[k]
                picked.append({"key": k, "depth": d, "uid": s["uid"], "pi_dives": s["pi"], "pi": pi,
                               "w": 1.0 / (s["pi"] * pi), "est_cpu": s["cpu"], "est_nodes": s["nodes"]})
    return picked


def settings_record(solver, sdir, params, flags, env):
    """the run's settings record as the solver writes it in a state file: its hash (what the solver
    checks on -r), then the policy and experiment settings away from their defaults"""
    r = subprocess.run([solver, "--print-settings", "-p", params] + flags, cwd=sdir, capture_output=True, text=True,
                       env=env)
    rec = None
    for line in r.stdout.splitlines():
        f = line.split()
        if f[:1] == ["settings"] and len(f) >= 4 and f[2] == "hash":
            rec = [f[3]]
        elif f[:1] == ["setting"] and rec is not None and len(f) >= 5 and f[-2] != "default" and \
                f[-1] in ("policy", "experiment"):
            rec.append(f"{f[1]}={' '.join(f[2:-2]).strip(chr(34))}")
    if rec is None:
        raise SystemExit(f"sample: no settings hash from {solver} --print-settings:\n{r.stdout[-500:]}{r.stderr[-500:]}")
    return rec


def sample_frontier(b, src, dv, plan, seed, mix, tmp, record, log=print, power=0.5):
    """write tmp/frontier.state (its parents in tmp/warm) with the sampled open nodes of SRC; returns the
    sample (draw()) with each node's line found"""
    s = b.get_run(src)
    if s is None:
        raise SystemExit(f"sample: no run {src}")
    busy = {k: v for k, v in b.count_by_status(src).items() if k in ("leased", "done", "returned", "failed") and v}
    if s.get("status") not in ("paused", "stopped") or busy:
        raise SystemExit(f"sample: {src} is {s.get('status')} with {busy or 'nothing'} in flight: a paused run, its units back")
    dv = dv or s.get("dives_run")
    if not dv or b.get_run(dv) is None:
        raise SystemExit(f"sample: no dives run of {src} to size its open nodes by (--sample-by)")
    strata = strata_sizes(b, dv)
    picked = draw(strata, plan, seed, mix, power)
    log(f"sample: {len(picked)} of {src}'s open nodes, by the estimates of {dv} ({len(strata)} strata), plan "
        + ", ".join(f"depth {d}: {n} wanted, {sum(1 for x in picked if x['depth'] == d)} taken" for d, n in sorted(plan.items()))
        + f", in proportion to the estimate^{power}, mix {mix}, seed {seed}")
    units = {u["uid"]: u for u in b.units_by_status(src, "queued")}
    want = {}
    for x in picked:
        want.setdefault(x["uid"], set()).add(x["key"])
    nodes, template, found = [], None, set()
    os.makedirs(tmp, exist_ok=True)
    for uid in sorted(want):
        u = units.get(uid)
        if u is None:
            raise SystemExit(f"sample: unit {uid} of {src} is not queued (was {src} resumed?)")
        path = os.path.join(tmp, f"src-{uid}.state")
        if not b.get_file(u["in_key"], path):
            raise SystemExit(f"sample: no {u['in_key']}")
        warm.fetch(b, src, path, log=log)          # its parents into tmp/warm, beside the frontier
        st = State(path)
        template = template or st
        for text, w in st.nodes:
            k = node_key(text, st.n)
            if k in want[uid]:
                nodes.append((text, w))
                found.add(k)
        os.remove(path)
    lost = [x["key"] for x in picked if x["key"] not in found]
    if lost:
        raise SystemExit(f"sample: {len(lost)} sampled nodes not in their units' state files (e.g. {lost[0]})")
    template.zero_stats()
    template.settings = record
    inc = s.get("incumbent") or {}
    if inc.get("obj") is not None and float(inc["obj"]) < template.obj:
        template.obj = float(inc["obj"])
        if inc.get("perm"):
            template.perm = [str(x) for x in inc["perm"]]
    frontier = os.path.join(tmp, "frontier.state")
    template.write(frontier, nodes)
    return frontier, picked, dv


def save_sample(b, run_id, picked, meta):
    b.put_bytes(run_key(run_id, "sample.json"), json.dumps({"meta": meta, "nodes": picked}, indent=1).encode())
