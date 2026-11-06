#!/usr/bin/env python3
"""Audit complete raw results and, when locally present, cached state bytes."""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def check(folder):
    m = read(folder / "manifest.json")
    done = read(folder / "completed.json")
    s = read(folder / "summary.json")
    assert sha(folder / "summary.json") == done["summary_sha256"]
    instance = ROOT / m["instance"]
    assert sha(instance) == m["instance_sha256"]
    binary = ROOT / "latest/dnnconditioned"
    if binary.exists():
        assert sha(binary) == m["binary_sha256"]
    locations = m["locations"]
    assert len(set(locations)) == len(locations) == s["children"]
    parent = read(folder / "parent.json")
    probes = {j: read(folder / f"probe_{j}.json") for j in locations}
    survivors = {j for j in locations if probes[j]["lb"] <= m["cutoff"]}
    assert {int(p.stem.split("_")[1]) for p in folder.glob("resume_*.json")} == survivors
    count = 0
    for mode, js in [("probe", locations), ("fresh", locations), ("resume", sorted(survivors))]:
        for j in js:
            path = folder / f"{mode}_{j}.json"
            r = read(path)
            raw = [json.loads(x) for x in path.with_suffix(".log").read_text().splitlines() if x.startswith("{")]
            assert len(raw) == 1 and all(r[k] == v for k, v in raw[0].items())
            assert r["mode"] == mode and r["facility"] == m["facility"] and r["location"] == j
            budget = m["probe"] if mode == "probe" else (m["rest"] if mode == "resume" else m["probe"] + m["rest"])
            cmd = r["command"]
            assert len(cmd) == 11 and cmd[1] == mode and cmd[3] == m["fixes"]
            assert Path(cmd[2]).name == instance.name
            assert list(map(int, cmd[4:6])) == [m["facility"], j]
            assert int(cmd[7]) == budget and int(cmd[8]) == m["eig"] and float(cmd[9]) == m["cutoff"]
            assert 0 < r["iterations"] <= budget and r["iterations"] % 10 == 0
            assert all(math.isfinite(r[k]) for k in ["lb", "raw_lb", "inherited_lb", "cpu", "wall", "margin"])
            assert r["lb"] == max(r["raw_lb"], r["inherited_lb"])
            assert r["lb"] >= parent["lb"]
            assert r["status"] in ["fathomed", "iter_limit"]
            assert (r["status"] == "fathomed") == (r["raw_lb"] > m["cutoff"])
            if mode == "resume":
                assert r["inherited_lb"] == probes[j]["lb"]
            if mode == "fresh" and j not in survivors:
                assert r["iterations"] == probes[j]["iterations"] and r["lb"] == probes[j]["lb"]
            count += 1
    for mode in ["parent", "control"]:
        r = read(folder / f"{mode}.json")
        raw = [json.loads(x) for x in (folder / f"{mode}.log").read_text().splitlines() if x.startswith("{")]
        assert len(raw) == 1 and all(r[k] == v for k, v in raw[0].items())
    control = read(folder / "control.json")
    assert control["status"] == "time_limit"
    assert float(control["command"][6]) == s["probe"]["cpu_including_parent"]
    assert control["cpu"] > float(control["command"][6])
    assert s["probe"]["lb"] == min(r["lb"] for r in probes.values())
    assert {r["location"] for r in s["pairs"]} == set(locations)
    for row in s["pairs"]:
        j = row["location"]
        assert row["probe"] == probes[j]["lb"]
        assert row["fresh"] == read(folder / f"fresh_{j}.json")["lb"]
        expected = read(folder / f"resume_{j}.json")["lb"] if j in survivors else probes[j]["lb"]
        assert row["reuse"] == expected
    cache_count = 0
    for j, record in read(folder / "cache_manifest.json").items():
        path = folder / "cache" / f"child_{j}.bin"
        if path.exists():
            assert path.stat().st_size == record["bytes"] and sha(path) == record["sha256"]
            cache_count += 1
    assert cache_count in [0, len(locations)]
    return {"case": folder.name, "children": len(locations), "raw_child_results_checked": count,
            "parent_and_control_logs_checked": 2, "verified_child_cache_hashes": cache_count,
            "binary_hash_checked": binary.exists(), "summary_hash_checked": True,
            "pruning_disagreements": [r["location"] for r in s["pairs"]
                if (r["reuse"] > m["cutoff"]) != (r["fresh"] > m["cutoff"])]}


if __name__ == "__main__":
    cases = [check(HERE / "results" / name) for name in ["nug20", "wil50_node9"]]
    result = {"passed": True, "cases": cases}
    (HERE / "results" / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
