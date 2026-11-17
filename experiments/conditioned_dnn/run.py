#!/usr/bin/env python3
"""Run one complete split, actual on-disk reuse, and charged control policies.

Every raw result is retained. A failed or partial family never produces a
conditioned bound. No production solver defaults are changed.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "latest/dnnconditioned"


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def run_one(out, name, args):
    command = [str(BIN), *map(str, args)]
    start = time.monotonic()
    with (out / (name + ".log")).open("w") as log:
        proc = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                              env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
                              timeout=7200)
    if proc.returncode:
        raise RuntimeError(f"{name} exited {proc.returncode}; see retained log")
    rows = [json.loads(s) for s in (out / (name + ".log")).read_text().splitlines() if s.startswith("{")]
    if len(rows) != 1:
        raise RuntimeError(f"{name}: expected one result")
    row = rows[0]
    row.update(command=command, process_wall=time.monotonic() - start)
    write(out / (name + ".json"), row)
    print(f"{name}: lb={row['lb']:.6f} it={row['iterations']} cpu={row['cpu']:.2f}", flush=True)
    return row


def summarize(out):
    m = json.loads((out / "manifest.json").read_text())
    locations = m["locations"]
    def read(name):
        return json.loads((out / (name + ".json")).read_text())
    parent = read("parent")
    probes = {j: read(f"probe_{j}") for j in locations}
    # Coverage is explicit, including original (not reduced) location labels.
    assert len(probes) == parent["m"] == len(set(locations))
    assert all(r["location"] == j and r["facility"] == m["facility"] for j, r in probes.items())
    cutoff = m["cutoff"]
    survivors = [j for j in locations if probes[j]["lb"] <= cutoff]
    resumed = {j: read(f"resume_{j}") for j in survivors}
    fresh = {j: read(f"fresh_{j}") for j in locations}
    reuse_lb = {j: resumed[j]["lb"] if j in resumed else probes[j]["lb"] for j in locations}
    discard_lb = {j: max(probes[j]["lb"], fresh[j]["lb"]) if j in resumed else probes[j]["lb"] for j in locations}
    for j in locations:
        assert probes[j]["lb"] >= parent["lb"]
        assert reuse_lb[j] >= probes[j]["lb"]
        assert fresh[j]["lb"] >= parent["lb"]
    p_cpu = sum(r["cpu"] for r in probes.values())
    r_cpu = sum(r["cpu"] for r in resumed.values())
    d_cpu = sum(fresh[j]["cpu"] for j in survivors)
    f_cpu = sum(r["cpu"] for r in fresh.values())
    def policy(bounds, cpu, iterations):
        return {"lb": min(bounds.values()), "bottleneck_locations": [j for j, b in bounds.items() if b == min(bounds.values())],
                "children_fathomed": sum(b > cutoff for b in bounds.values()),
                "cpu_including_parent": parent["cpu"] + cpu, "child_iterations": iterations}
    pit = sum(r["iterations"] for r in probes.values())
    summary = {
        "complete": True, "children": len(locations), "cutoff": cutoff,
        "parent": parent,
        "probe": policy({j: r["lb"] for j, r in probes.items()}, p_cpu, pit),
        "reuse": policy(reuse_lb, p_cpu+r_cpu, pit+sum(r["iterations"] for r in resumed.values())),
        "discard": policy(discard_lb, p_cpu+d_cpu, pit+sum(fresh[j]["iterations"] for j in survivors)),
        "direct": policy({j: r["lb"] for j, r in fresh.items()}, f_cpu, sum(r["iterations"] for r in fresh.values())),
        "matched_parent_control": read("control"),
        "saved_child_state_bytes": sum(r["state_bytes"] for r in probes.values()),
        "probe_snapshot_write_wall": sum(r["save_wall"] for r in probes.values()),
        "cert_failures": sum(r["cert_failures"] for r in [parent,*probes.values(),*resumed.values(),*fresh.values()]),
        "pairs": [{"location": j, "probe": probes[j]["lb"], "reuse": reuse_lb[j],
                   "fresh": fresh[j]["lb"], "reuse_minus_fresh": reuse_lb[j]-fresh[j]["lb"],
                   "probe_iterations": probes[j]["iterations"],
                   "resume_iterations": resumed[j]["iterations"] if j in resumed else 0,
                   "fresh_iterations": fresh[j]["iterations"]} for j in locations],
    }
    write(out / "summary.json", summary)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--instance", required=True)
    p.add_argument("--fixes", default="-")
    p.add_argument("--facility", type=int, required=True)
    p.add_argument("--cutoff", type=float, required=True)
    p.add_argument("--parent", type=int, default=200)
    p.add_argument("--probe", type=int, default=50)
    p.add_argument("--rest", type=int, default=100)
    p.add_argument("--eig", type=int, default=2)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "manifest.json").exists():
        raise RuntimeError("use a new output directory; never overwrite an attempt")
    instance = (ROOT / a.instance).resolve()
    n = int(instance.read_text().split()[0])
    fixed = [] if a.fixes == "-" else [tuple(map(int, x.split(":"))) for x in a.fixes.split(",")]
    assert a.facility not in [i for i, j in fixed]
    locations = [j for j in range(n) if j not in [k for i, k in fixed]]
    cache = out / "cache"
    cache.mkdir()
    manifest = {**vars(a), "out": str(out), "locations": locations, "n": n,
                "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
                "binary_sha256": digest(BIN), "instance_sha256": digest(instance),
                "platform": platform.platform(), "start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "control": "cold parent; CPU budget equals parent+complete probe solver CPU; checkpoint overshoot reported",
                "settings": "PRSM gamma1=gamma=1 beta=.3 refresh=40 degree=3 tight_lam=2 exact LAP check=10; tol/hopeless/stall disabled; one BLAS thread"}
    write(out / "manifest.json", manifest)
    common = [instance, a.fixes]
    parent = run_one(out, "parent", ["parent", *common, a.parent, a.eig, a.cutoff, cache / "parent.bin"])
    parent_hash = digest(cache / "parent.bin")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        jobs = {pool.submit(run_one, out, f"probe_{j}", ["probe", *common, a.facility, j, cache / "parent.bin",
                a.probe, a.eig, a.cutoff, cache / f"child_{j}.bin"]): j for j in locations}
        probes = {jobs[f]: f.result() for f in cf.as_completed(jobs)}
    assert digest(cache / "parent.bin") == parent_hash
    write(out / "cache_manifest.json", {str(j): {"sha256": digest(cache / f"child_{j}.bin"),
         "bytes": (cache / f"child_{j}.bin").stat().st_size} for j in locations})
    probe_lb = min(r["lb"] for r in probes.values())
    budget = parent["cpu"] + sum(r["cpu"] for r in probes.values())
    write(out / "probe_complete.json", {"lb": probe_lb, "parent_lb": parent["lb"], "cpu": budget,
          "children": len(locations), "fathomed": sum(r["lb"] > a.cutoff for r in probes.values())})
    print(f"COMPLETE PROBE: lb={probe_lb:.6f} gain={probe_lb-parent['lb']:.6f} CPU={budget:.2f}", flush=True)
    # Control and child comparisons use at most jobs+1 one-thread solver processes.
    # CPU totals, not elapsed time of this concurrent experiment, price policies.
    with cf.ThreadPoolExecutor(max_workers=1) as control_pool, cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        control = control_pool.submit(run_one, out, "control", ["control", *common, 100000, a.eig, budget, "-"])
        def pair(j):
            modes = ["resume", "fresh"] if j%2==0 else ["fresh", "resume"]
            for mode in modes:
                if mode=="resume" and probes[j]["lb"]>a.cutoff:
                    continue
                source = cache / (f"child_{j}.bin" if mode=="resume" else "parent.bin")
                iters = a.rest if mode=="resume" else a.probe+a.rest
                run_one(out, f"{mode}_{j}", [mode, *common, a.facility, j, source, iters, a.eig, a.cutoff, "-"])
        for f in cf.as_completed([pool.submit(pair,j) for j in locations]):
            f.result()
        control.result()
    summary = summarize(out)
    write(out / "completed.json", {"end_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                   "summary_sha256": digest(out / "summary.json")})
    for k in ["probe", "reuse", "discard", "direct"]:
        print(k, json.dumps(summary[k]), flush=True)


if __name__ == "__main__":
    main()
