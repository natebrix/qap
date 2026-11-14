"""The parent-state store of a coordinator run (docs/parent_state_store.md, the
coordinator's half). A child's DNN bound starts from its parent's final state
(two (m^2+1)^2 matrices; ship.bits sets their precision). Written inline, every
open node of a state file carries its parent's copy, so siblings repeat it: at
m = 40 and ship.bits = 64 that was 55 MB of base64 a node, and tho40-a's
checkpoints filled a 32 GB disk in 30 minutes (Nov 14). With the store each
parent is written once, to runs/RUN/warm/<fnv64>.dnn, and a state file's warm
line refers to it as "@warm/<name>":

- the solver writes the store when QAP_DNN_WARM_DIR is set, the reference
  relative to the state file (latest/state.C dnn_field); the frontier pass and
  every unit attempt run with QAP_DNN_WARM_DIR=<the state file's dir>/warm;
- a worker fetches the parents a state file names before the solver reads it
  (fetch) and uploads the ones a file it uploads names (push), each once;
- the coordinator never reads a parent: its state files keep the references
  (rebase makes them "@warm/<name>" wherever the file was read), and each unit
  records the parents its nodes start from (node_refs; a checkpoint's are added
  when it is uploaded);
- gc deletes the parents no unit not yet merged refers to and no attempt still
  running can have written. A parent that is missing all the same costs a cold
  start, never the unit: fetch rewrites its references to "-".
"""
import os
import re
import time

from .common import run_key

REF = re.compile(rb"@warm/([0-9a-f]{16}\.dnn)")


def warm_key(run_id, name):
    return run_key(run_id, "warm", name)


def refs_in(path):
    """the store files a state file refers to"""
    try:
        with open(path, "rb") as f:
            return sorted({m.decode() for m in REF.findall(f.read())})
    except FileNotFoundError:
        return []


def fetch(b, run_id, state_path, log=print):
    """bring the parents a state file refers to into warm/ beside it; a parent
    the store has lost is dropped from the file (its nodes start cold)"""
    warm_dir = os.path.join(os.path.dirname(os.path.abspath(state_path)), "warm")
    os.makedirs(warm_dir, exist_ok=True)
    missing = []
    for name in refs_in(state_path):
        local = os.path.join(warm_dir, name)
        if not os.path.exists(local) and not b.get_file(warm_key(run_id, name), local):
            missing.append(name)
    if missing:
        with open(state_path, "rb") as f:
            data = f.read()
        for name in missing:
            data = data.replace(b"@warm/" + name.encode(), b"-")
        tmp = state_path + ".tmp"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, state_path)
        log(f"warm: {len(missing)} parent state(s) missing from the store for {state_path}: those nodes start cold")
    return missing


def push(b, run_id, state_path, pushed, log=print):
    """upload the parents a state file refers to that the store lacks (pushed:
    the names this process already knows are there)"""
    warm_dir = os.path.join(os.path.dirname(os.path.abspath(state_path)), "warm")
    names = refs_in(state_path)
    for name in names:
        if name in pushed:
            continue
        local = os.path.join(warm_dir, name)
        key = warm_key(run_id, name)
        if os.path.exists(local):
            if not b.exists(key):
                b.put_file(key, local)
        elif not b.exists(key):
            log(f"warm: {state_path} refers to {name}, which is neither here nor in the store")
            continue
        pushed.add(name)
    return names


def rebase(state, base_dir):
    """make a read state file's references point into base_dir/warm, so that
    written in base_dir they read "@warm/<name>" again"""
    for _, w in state.nodes:
        if w is not None and w.ref:
            w.ref = os.path.join(base_dir, "warm", os.path.basename(w.ref))
    return state


def node_refs(nodes):
    """the store files a list of nodes starts from"""
    return sorted({os.path.basename(w.ref) for _, w in nodes if w is not None and w.ref})


def gc(b, run_id, units, now=None, margin=300):
    """delete the parents nothing can need: kept are those a unit not yet merged
    starts from or its last uploaded checkpoint names, and those written since
    the oldest attempt still running or waiting to be merged began (its output
    may name them), less a margin. Returns (kept, deleted)."""
    now = now or time.time()
    live, since = set(), now
    for u in units:
        if u.get("status") == "merged":
            continue
        live.update(u.get("warm", []))
        live.update(u.get("ckpt_warm", []))
        if u.get("status") in ("leased", "done", "returned") and u.get("leased_at"):
            since = min(since, float(u["leased_at"]))
    kept = deleted = 0
    for key, mtime in b.list_meta(run_key(run_id, "warm") + "/"):
        if os.path.basename(key) in live or mtime >= since - margin:
            kept += 1
            continue
        b.delete(key)
        deleted += 1
    return kept, deleted
