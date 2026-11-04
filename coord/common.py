"""Shared pieces of the coordinator and the worker: state files (the State
class of bench/parallel.py), the instance, incumbent verification, unit
documents and the split of open nodes into units."""
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
from parallel import State  # noqa: E402  (the solver's state-file format, latest/state.h)

SOLVER = os.path.join(ROOT, "latest", "bnbdrv")
LANES = ("large", "small")


def run_key(run_id, *parts):
    return "/".join(("runs", run_id) + parts)


def unit_key(run_id, uid, *parts):
    return run_key(run_id, "units", uid, *parts)


# ---------------------------------------------------------------- instance

def read_instance(path_or_text):
    text = path_or_text
    if os.path.exists(path_or_text):
        with open(path_or_text) as f:
            text = f.read()
    t = text.split()
    n = int(t[0])
    vals = [int(float(x)) for x in t[1:]]
    A = [vals[i * n:(i + 1) * n] for i in range(n)]
    B = [vals[n * n + i * n:n * n + (i + 1) * n] for i in range(n)]
    C = None
    if len(vals) >= 3 * n * n:
        C = [vals[2 * n * n + i * n:2 * n * n + (i + 1) * n] for i in range(n)]
    return n, A, B, C


def objective(inst, perm):
    """exact integer objective sum_ij A[i][j] B[p_i][p_j] (+ C[i][p_i]) of the
    original data, the solver's convention for the state file's incumbent"""
    n, A, B, C = inst
    v = sum(A[i][j] * B[perm[i]][perm[j]] for i in range(n) for j in range(n))
    if C:
        v += sum(C[i][perm[i]] for i in range(n))
    return v


def verify_incumbent(inst, obj, perm):
    """(ok, message): perm a permutation of 0..n-1 whose objective equals obj"""
    n = inst[0]
    try:
        p = [int(x) for x in perm]
    except (TypeError, ValueError):
        return False, "permutation is not integers"
    if sorted(p) != list(range(n)):
        return False, "not a permutation"
    v = objective(inst, p)
    if abs(v - obj) > 1e-9 * max(1, abs(obj)):
        return False, f"objective {v} recomputed, {obj} reported"
    return True, "ok"


# ---------------------------------------------------------------- nodes and units

def node_depth(node):
    return int(node[0].split()[1])


def node_bound(node):
    return float(node[0].split()[2])


def est_node_cost(totals, depth, default=60.0):
    """mean CPU seconds per node at this depth so far (t_level / n_level)"""
    try:
        n = totals.arrays["n_level"][depth]
        t = totals.arrays["t_level"][depth]
        if n > 0 and t > 0:
            return t / n
    except (KeyError, IndexError):
        pass
    return default


def plan_split(nodes, totals, cfg, queued):
    """chunks of consecutive nodes (stack order keeps siblings, which share a
    parent's warm start, together): as many units as the CPU target asks for,
    more when the queue holds fewer than reserve x slots units"""
    N = len(nodes)
    if N == 0:
        return []
    cost = sum(est_node_cost(totals, node_depth(x)) for x in nodes)
    by_cost = math.ceil(cost / max(1.0, cfg["target_unit_cpu"]))
    need = max(0, int(cfg["reserve"] * cfg["slots"]) - queued)
    units = max(1, min(N, max(by_cost, min(N, need))))
    k = math.ceil(N / units)
    return [nodes[i:i + k] for i in range(0, N, k)]


def new_unit(run_id, uid, parent, nodes, cfg, now=None):
    depths = [node_depth(x) for x in nodes]
    return {"uid": uid, "run": run_id, "parent": parent, "status": "queued", "attempt": 0,
            "nodes": len(nodes), "depth_min": min(depths), "depth_max": max(depths),
            "min_bound": min(node_bound(x) for x in nodes),
            "lane": "large" if min(depths) <= cfg.get("large_depth", 4) else "small",
            "budget": cfg["budget"], "in_key": unit_key(run_id, uid, "in.state"),
            "created": now or time.time()}


def write_unit_state(template, nodes, path):
    template.write(path, nodes)


def fmt_secs(s):
    s = int(s)
    return f"{s // 3600}h{s % 3600 // 60:02d}m{s % 60:02d}s"
