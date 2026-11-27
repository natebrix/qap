#!/usr/bin/env python3
"""qapctl unpublish (coord/publish.py) on a local store, without the solver: the run leaves the
catalog (both copies) and loses its public files, is marked private and keeps its record; another
run is untouched; a running run is refused; an unknown run is an error.

  python3 tests/unpublish_test.py [BACKEND]     (BACKEND: a set-up store, e.g. aws:... on a moto server;
                                                 default: a new local one)
"""
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from coord.backend import open_backend          # noqa: E402
from coord.publish import CATALOG, PUBLIC, unpublish   # noqa: E402
from coord import qapctl                        # noqa: E402


def main():
    spec = sys.argv[1] if len(sys.argv) > 1 else "local:" + tempfile.mkdtemp(prefix="qap-unpublish-")
    b = open_backend(spec)
    for run, status in (("stray", "finished"), ("stray-2", "finished"), ("live", "running")):
        b.put_run(run, {"run_id": run, "status": status, "visibility": "public"})
        b.put_bytes(f"{PUBLIC}/{run}/status.json", b"{}")
        b.put_bytes(f"{PUBLIC}/{run}/tree/meta.json", b"{}")
        b.put_bytes(f"runs/{run}/results.txt", b"kept")
    body = json.dumps({"schema": 1, "generated_at": "x", "runs": [{"run": r} for r in ("live", "stray-2", "stray")]}).encode()
    b.put_bytes(CATALOG, body)
    b.put_bytes(f"{PUBLIC}/index.json", body)
    problems = []
    said = []
    if qapctl.main(["unpublish", "stray", "--backend", spec]) != 0:
        problems.append("unpublish stray: not 0")
    for key in (CATALOG, f"{PUBLIC}/index.json"):
        names = [r["run"] for r in json.loads(b.get_bytes(key))["runs"]]
        if names != ["live", "stray-2"]:
            problems.append(f"{key}: {names}")
    if b.list(f"{PUBLIC}/stray/"):
        problems.append(f"public files left: {b.list(f'{PUBLIC}/stray/')}")
    if len(b.list(f"{PUBLIC}/stray-2/")) != 2 or len(b.list(f"{PUBLIC}/live/")) != 2:
        problems.append("another run's public files were touched")
    if b.get_bytes("runs/stray/results.txt") != b"kept":
        problems.append("the run's record is gone")
    if b.get_run("stray").get("visibility") != "private" or b.get_run("stray-2").get("visibility") != "public":
        problems.append("visibility")
    if unpublish(b, "live", log=said.append) != 1 or [r["run"] for r in json.loads(b.get_bytes(CATALOG))["runs"]] != ["live", "stray-2"]:
        problems.append(f"a running run was not refused: {said}")
    if unpublish(b, "nosuch", log=said.append) != 1:
        problems.append("an unknown run was not an error")
    if unpublish(b, "stray", log=said.append) != 0:      # again: nothing left to do, still fine
        problems.append("a second unpublish failed")
    print("\n".join(problems) if problems else "unpublish: PASS")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
