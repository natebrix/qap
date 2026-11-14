#!/usr/bin/env python3
"""coord/warm.py, the coordinator's parent-state store: what gc keeps and
deletes, and fetch's cold start for a parent the store has lost.

  python3 tests/warm_test.py
"""
import os
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from coord import warm  # noqa: E402
from coord.backend import open_backend  # noqa: E402


def main():
    tmp = tempfile.mkdtemp(prefix="warmtest-")
    b = open_backend(f"local:{os.path.join(tmp, 'store')}")
    run = "r"
    names = {k: f"{i:016x}.dnn" for i, k in enumerate("abcdefg")}
    now = time.time()
    for k, name in names.items():
        b.put_bytes(warm.warm_key(run, name), b"dnn32:1:1:0:AAAA\n")
        old = now - 3600 if k != "g" else now - 30          # g: written a moment ago
        os.utime(b._path(warm.warm_key(run, name)), (old, old))
    units = [
        {"uid": "u0", "status": "merged", "warm": [names["a"]]},                  # merged: a is free
        {"uid": "u1", "status": "queued", "warm": [names["b"]]},                  # queued: b is kept
        {"uid": "u2", "status": "leased", "warm": [names["c"]], "ckpt_warm": [names["d"]],
         "leased_at": now - 7200},                                               # c, d kept; and since 2 h ago
        {"uid": "u3", "status": "done", "warm": [], "leased_at": now - 600},      # e, f written within its attempt
    ]
    bad = []
    kept, deleted = warm.gc(b, run, units, now=now)
    left = {os.path.basename(k) for k in b.list(f"runs/{run}/warm/")}
    # the leased unit began 2 h ago, so everything written since then is kept
    if left != set(names.values()):
        bad.append(f"with an attempt running since 2 h ago, kept {sorted(left)}")
    units[2]["leased_at"] = now - 60                                              # now only a minute ago
    kept, deleted = warm.gc(b, run, units, now=now)
    left = {os.path.basename(k) for k in b.list(f"runs/{run}/warm/")}
    want = {names[k] for k in "bcdg"}                                             # e, f older than u3's lease less 5 min
    if left != want:
        bad.append(f"kept {sorted(left)}, want {sorted(want)}")
    # fetch: a lost parent is dropped from the state file, the others fetched
    sdir = os.path.join(tmp, "unit")
    os.makedirs(sdir)
    st = os.path.join(sdir, "in.state")
    with open(st, "w") as f:
        f.write(f"warm 3 0 1 2 0 0 - @warm/{names['b']}\nwarm 3 0 1 2 0 0 - @warm/{names['a']}\n")
    missing = warm.fetch(b, run, st, log=lambda *_: None)
    text = open(st).read()
    if missing != [names["a"]] or f"@warm/{names['a']}" in text or f"@warm/{names['b']}" not in text \
            or not os.path.exists(os.path.join(sdir, "warm", names["b"])):
        bad.append(f"fetch: missing {missing}, file {text!r}")
    print("warm: " + ("all checks passed" if not bad else "FAIL\n  " + "\n  ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
