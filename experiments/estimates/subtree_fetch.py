#!/usr/bin/env python3
"""subtree_fetch.py BACKEND RUN DIR: copy a finished subtree run's results.txt and
its units' trace lines out of a coordinator backend into DIR (results.txt,
trace.txt), as subtree.sh does from a local backend's directory. For a run on
another backend (subtree.sh with EST_BACKEND, e.g. aws:BUCKET?region=R).
Exit status 1 when the run has no results yet (it has not ended)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from coord.backend import open_backend
from coord.common import run_key

backend, run, d = sys.argv[1:4]
b = open_backend(backend)
res = b.get_bytes(run_key(run, "results.txt"))
if res is None:
    sys.exit(1)
with open(os.path.join(d, "results.txt"), "wb") as f:
    f.write(res)
with open(os.path.join(d, "trace.txt"), "wb") as f:
    for k in b.list(run_key(run, "units") + "/"):
        if k.endswith("/trace.txt"):
            f.write(b.get_bytes(k) or b"")
