"""Run two full cut-phase A/B pairs; retain logs, process CPU, and provenance.

Run from any directory. The output directory must not already exist.
Four single-thread processes run concurrently. Linux runs pin distinct physical
cores; macOS runs leave placement to the operating system.
The only solver option differing within each pair is -warmY.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
NODES = [
    ("wil50_n13", "wil50", "21:28,8:41,28:19", "5:10", 48816),
    ("sko49_d4_13", "sko49", "7:48,26:41,2:12,44:30", "7:7", 23386),
]
COMMON = ["-side", "0", "-T", "1", "-k", "1000", "-v", "-I", "400",
          "-i", "300", "-r", "12", "-delta", "50", "-nocap", "-freeze",
          "0:1500", "-stall", "0", "-e", "2", "-d", "4", "-g", "1", "-G", "1"]


def parse_log(path):
    text = path.read_text()
    rounds = []
    for line in text.splitlines():
        m = re.match(r"(round\s+\d+|quiet round)\s+lb (\S+).*?status (\S+)\s+iters (\d+)\s+cpu (\S+)", line)
        if m:
            rounds.append(dict(round=m[1], lb=float(m[2]), status=m[3],
                               iterations=int(m[4]), dnn_cpu=float(m[5])))
    final = re.search(r"^FINAL best_lb (\S+) first_lb (\S+) gain (\S+)", text, re.M)
    return dict(rounds=rounds, iterations=sum(r["iterations"] for r in rounds),
                dnn_cpu=sum(r["dnn_cpu"] for r in rounds),
                fathomed=any(r["status"] == "fathomed" for r in rounds),
                solver_error=any(r["status"] == "error" for r in rounds),
                final=None if not final else dict(zip(("best_lb", "first_lb", "gain"), map(float, final.groups()))))


def physical_cores():
    if not hasattr(os, "sched_getaffinity"):
        return [None] * 4
    selected, seen = [], set()
    for cpu in sorted(os.sched_getaffinity(0)):
        topo = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key = ((topo / "physical_package_id").read_text().strip(),
               (topo / "core_id").read_text().strip())
        if key not in seen:
            selected.append(cpu)
            seen.add(key)
    if len(selected) < 4:
        raise RuntimeError("Four distinct physical cores are required")
    return selected[:4]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=HERE / "results" / "full_cut_phase_ab")
    ap.add_argument("--hours", type=float, default=3.75, help="wall deadline; max 3.75 hours")
    args = ap.parse_args()
    if not 0 < args.hours <= 3.75:
        ap.error("--hours must be in (0, 3.75]")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = ROOT / "experiments/dnn_cuts/dnncuts"
    cores = physical_cores()
    start = time.monotonic()
    manifest = dict(commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    wall_limit_seconds=args.hours * 3600, jobs=[],
                    timing_note="Four concurrent single-thread processes; pinned to distinct physical cores when CPU IDs are recorded, otherwise OS-scheduled. Process CPU includes cut work. Shared-memory contention and core differences remain timing confounders.")
    running = []
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

    def persist():
        manifest["elapsed_wall_seconds"] = time.monotonic() - start
        temp = out / "manifest.json.tmp"
        temp.write_text(json.dumps(manifest, indent=2) + "\n")
        temp.replace(out / "manifest.json")

    try:
        for tag, instance, fixes, grid, ub in NODES:
            for arm in ("Z_only", "YZ"):
                cpu = cores[len(running)]
                command = [str(binary), "-x", fixes, "-u", str(ub), "-grid", grid, *COMMON]
                if arm == "YZ":
                    command.append("-warmY")
                command.append(f"qaplib/{instance}.dat")
                log = out / f"{tag}_{arm}.log"
                handle = log.open("w")
                p = subprocess.Popen(command, cwd=ROOT, env=env, stdout=handle,
                                     stderr=subprocess.STDOUT, start_new_session=True,
                                     preexec_fn=None if cpu is None else lambda cpu=cpu: os.sched_setaffinity(0, {cpu}))
                job = dict(node=tag, arm=arm, cpu=cpu, command=command, log=log.name,
                           state="running", pid=p.pid)
                manifest["jobs"].append(job)
                running.append((p, handle, log, job, time.monotonic()))
                print(f"START {tag} {arm} cpu={cpu}", flush=True)
        persist()
        last_report = 0
        while running:
            for item in list(running):
                p, handle, log, job, t0 = item
                pid, status, usage = os.wait4(p.pid, os.WNOHANG)
                if pid:
                    p.returncode = os.waitstatus_to_exitcode(status)
                    handle.close()
                    result = parse_log(log)
                    job.update(result, exit_code=p.returncode, wall_seconds=time.monotonic()-t0,
                               process_cpu_seconds=usage.ru_utime+usage.ru_stime,
                               max_rss_kib=usage.ru_maxrss / (1024 if sys.platform == "darwin" else 1),
                               state="completed" if p.returncode == 0 and result["final"] and not result["solver_error"] else "failed")
                    running.remove(item)
                    print(f"DONE {job['node']} {job['arm']} {job['state']} {job['final']}", flush=True)
                    persist()
            now = time.monotonic()
            if running and now-start >= args.hours*3600:
                raise TimeoutError("Experiment reached its declared wall deadline")
            if now-last_report >= 60:
                for _, _, log, job, _ in running:
                    result = parse_log(log)
                    job["progress"] = result
                    print(f"PROGRESS {job['node']} {job['arm']} iterations={result['iterations']} rounds={len(result['rounds'])}", flush=True)
                persist()
                last_report = now
            if running:
                time.sleep(2)
    finally:
        for p, handle, log, job, t0 in running:
            os.killpg(p.pid, signal.SIGKILL)
            _, status, usage = os.wait4(p.pid, 0)
            p.returncode = os.waitstatus_to_exitcode(status)
            handle.close()
            job.update(parse_log(log), state="terminated", exit_code=p.returncode,
                       wall_seconds=time.monotonic()-t0,
                       process_cpu_seconds=usage.ru_utime+usage.ru_stime,
                       max_rss_kib=usage.ru_maxrss / (1024 if sys.platform == "darwin" else 1))
        persist()
    if any(j["state"] != "completed" for j in manifest["jobs"]):
        raise RuntimeError("An arm failed; inspect the retained logs")
    print(f"COMPLETE {out}", flush=True)


if __name__ == "__main__":
    main()
