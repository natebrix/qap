"""Predeclared six-node, five-arm continuation-budget experiment (Linux)."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from run_cut_phase_ab import parse_log

CUTOFF = 48814.05
BASE = "4095416756ef487251912396a727fac7ceacacaf"
NODES = [
    dict(id="n14", fixes="47:23,10:44,40:34", near=True, line=14),
    dict(id="n6", fixes="20:12,40:28,17:25", near=True, line=6),
    dict(id="n9", fixes="26:27,8:38,34:28", near=True, line=9),
    dict(id="n13", fixes="21:28,8:41,28:19", near=False, line=13),
    dict(id="d2", fixes="39:27,30:48,25:49,26:5", near=False, line=None),
    dict(id="n10", fixes="39:27,30:48,25:49", near=False, line=10),
]
ARMS = dict(B=(12, 300, 1500, False), Y12=(12, 300, 1500, True),
            Y6=(6, 300, 1500, True), Y8s=(8, 300, 600, True),
            Y12s=(12, 200, 1000, True))
COMMON = ["-grid", "5:10", "-u", "48816", "-T", "1", "-k", "1000",
          "-v", "-I", "400", "-delta", "50", "-nocap"]


def verify_nodes():
    lines = (ROOT / "experiments/wil50_feasibility/tuning/giveup/all32.txt").read_text().splitlines()
    for node in NODES:
        if node["line"] is not None:
            assert lines[node["line"]-1].strip() == node["fixes"], node
    row = (ROOT / "experiments/wil50_feasibility/stronger_bounds/depth4/nodes.txt").read_text().splitlines()[1].split()
    assert row[0] == "2" and row[2] == NODES[4]["fixes"]


def plan():
    verify_nodes()
    # Each wave covers all nodes and rotates arm assignment. Shuffle only
    # within each wave, before any outcomes are seen; no arm runs exclusively
    # at the start or end of the study.
    rng = random.Random(20261102)
    jobs = []
    names = list(ARMS)
    for wave in range(5):
        block = []
        for i, node in enumerate(NODES):
            arm = names[(i+wave) % 5]
            rounds, iters, quiet, warm = ARMS[arm]
            command = ["experiments/dnn_cuts/dnncuts", *COMMON, "-x", node["fixes"],
                       "-r", str(rounds), "-i", str(iters), "-freeze", f"0:{quiet}"]
            if warm:
                command.append("-warmY")
            command.append("qaplib/wil50.dat")
            block.append(dict(node=node["id"], near=node["near"], fixes=node["fixes"],
                              arm=arm, wave=wave, command=command,
                              max_iterations=400+rounds*iters+quiet))
        rng.shuffle(block)
        jobs.extend(block)
    assert len({(j['node'], j['arm']) for j in jobs}) == 30
    return dict(base_cognac=BASE, cutoff=CUTOFF, schedule_seed=20261102,
                max_workers=8, wall_limit_hours=6, nodes=NODES, arms=ARMS, jobs=jobs)


def parse_run(path, near):
    result = parse_log(path)
    cumulative = 0
    first_cross = None
    polyak = []
    quiet = None
    for r in result["rounds"]:
        cumulative += r["iterations"]
        r["cumulative_iterations"] = cumulative
        if first_cross is None and (r["status"] == "fathomed" or r["lb"] > CUTOFF):
            first_cross = cumulative
        if r["round"] == "quiet round":
            quiet = r["lb"]
        else:
            index = int(r["round"].split()[-1])
            if index > 0:
                polyak.append((r["lb"], -index))
    best = max(polyak) if polyak else None
    reported = result.pop("fathomed")
    crossed = reported or first_cross is not None or (result["final"] is not None and result["final"]["best_lb"] > CUTOFF)
    result.update(best_polyak_round=None if best is None else -best[1],
                  best_polyak_bound=None if best is None else best[0], quiet_bound=quiet,
                  first_crossing_iteration=first_cross, raw_reported_fathom=reported,
                  validity_question=bool(crossed and not near),
                  outcome="validity_question" if crossed and not near else "fathomed" if crossed else "short")
    return result


def cores(count):
    selected, seen = [], set()
    for cpu in sorted(os.sched_getaffinity(0)):
        p = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key = ((p / "physical_package_id").read_text(), (p / "core_id").read_text())
        if key not in seen:
            selected.append(cpu)
            seen.add(key)
    if len(selected) < count:
        raise RuntimeError("Insufficient distinct physical cores")
    quota = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    if quota[0] != 'max' and count > int(quota[0])/int(quota[1]):
        raise RuntimeError("Worker count exceeds CPU quota")
    return selected[:count]


def atomic_json(path, data):
    def clean(value):
        if isinstance(value, float) and not math.isfinite(value):
            return repr(value)
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(v) for v in value]
        return value
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(clean(data), indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plan", action="store_true", help="write and print the fixed plan; launch nothing")
    ap.add_argument("--out", type=Path, default=HERE / "results" / "runs")
    ap.add_argument("--workers", type=int, default=8, choices=range(1, 9))
    args = ap.parse_args()
    definition = plan()
    if args.plan:
        atomic_json(HERE / "plan.json", definition)
        print(json.dumps([(j['node'], j['arm']) for j in definition['jobs']]))
        return
    assert json.loads((HERE / "plan.json").read_text()) == json.loads(json.dumps(definition)), "Predeclared plan changed"
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    git = lambda *a: subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
    subprocess.run(["git", "merge-base", "--is-ancestor", BASE, "HEAD"], cwd=ROOT, check=True)
    assert not git("diff", "HEAD", "--", "latest", "experiments/dnn_cuts"), "Dirty solver source"
    cpus = cores(args.workers)
    binary = ROOT / "experiments/dnn_cuts/dnncuts"
    started = time.monotonic()
    manifest = dict(source_commit=git("rev-parse", "HEAD"), base_cognac=BASE,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    instance_sha256=hashlib.sha256((ROOT/'qaplib/wil50.dat').read_bytes()).hexdigest(),
                    plan_sha256=hashlib.sha256((HERE / 'plan.json').read_bytes()).hexdigest(),
                    started_utc=datetime.now(timezone.utc).isoformat(), workers=args.workers,
                    cores=cpus, wall_limit_seconds=21600, jobs=[],
                    timing_note="Single observations under a balanced predeclared queue, at most eight concurrent single-thread processes pinned to distinct physical cores. Shared bandwidth, core differences and tail concurrency affect CPU comparisons.")
    for spec in definition['jobs']:
        manifest['jobs'].append(dict(spec, state="queued", log=f"{spec['node']}_{spec['arm']}.log"))
    running = {}
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    events = (out / "events.jsonl").open("w", buffering=1)

    def persist():
        manifest['elapsed_wall_seconds'] = time.monotonic()-started
        atomic_json(out / "manifest.json", manifest)

    def event(kind, job, **extra):
        events.write(json.dumps(dict(event=kind, elapsed=time.monotonic()-started,
                                     node=job['node'], arm=job['arm'], **extra))+'\n')

    def finish(pid, status, usage):
        p, handle, job, t0 = running.pop(pid)
        p.returncode = os.waitstatus_to_exitcode(status)
        handle.close()
        result = parse_run(out / job['log'], job['near'])
        good = p.returncode == 0 and result['final'] is not None and not result['solver_error']
        if good:
            good = all(math.isfinite(v) for v in result['final'].values())
        job.pop('progress', None)
        job.update(result, state="completed" if good else "failed", exit_code=p.returncode,
                   process_cpu_seconds=usage.ru_utime+usage.ru_stime,
                   wall_seconds=time.monotonic()-t0, max_rss_kib=usage.ru_maxrss)
        event("finished", job, outcome=job['outcome'], exit_code=p.returncode)
        print(f"DONE {job['node']} {job['arm']} {job['state']} {job['outcome']} {job['final']}", flush=True)

    try:
        last_report = 0
        while running or any(j['state']=='queued' for j in manifest['jobs']):
            elapsed = time.monotonic()-started
            if elapsed >= 21600:
                raise TimeoutError("Six-hour predeclared runtime limit reached")
            used = {v[2]['cpu'] for v in running.values()}
            for cpu in cpus:
                if cpu in used:
                    continue
                job = next((j for j in manifest['jobs'] if j['state']=='queued'), None)
                if job is None:
                    break
                handle = (out / job['log']).open('w')
                p = subprocess.Popen(job['command'], cwd=ROOT, env=env, stdout=handle,
                                     stderr=subprocess.STDOUT, start_new_session=True,
                                     preexec_fn=lambda cpu=cpu: os.sched_setaffinity(0, {cpu}))
                job.update(state='running', cpu=cpu, started_elapsed=elapsed)
                running[p.pid] = (p, handle, job, time.monotonic())
                event('started', job, cpu=cpu)
                print(f"START {job['node']} {job['arm']} core={cpu}", flush=True)
            for pid in list(running):
                got, status, usage = os.wait4(pid, os.WNOHANG)
                if got:
                    finish(pid, status, usage)
            for p, _, job, _ in running.values():
                result = parse_run(out / job['log'], job['near'])
                old = len(job.get('progress', {}).get('rounds', []))
                if len(result['rounds']) > old:
                    for r in result['rounds'][old:]:
                        event('round_observed', job, round=r)
                    job['progress'] = result
                    if result['validity_question']:
                        print(f"VALIDITY QUESTION: {job['node']} {job['arm']} crossed cutoff", flush=True)
            persist()
            if elapsed-last_report >= 60:
                done = sum(j['state'] in ('completed', 'failed') for j in manifest['jobs'])
                print(f"PROGRESS minutes={elapsed/60:.1f} done={done}/30 running={len(running)}", flush=True)
                last_report = elapsed
            if running:
                time.sleep(5)
    finally:
        for pid in list(running):
            os.killpg(pid, signal.SIGKILL)
            _, status, usage = os.wait4(pid, 0)
            finish(pid, status, usage)
            manifest['interrupted'] = True
        persist()
        events.close()
    if any(j['state'] != 'completed' for j in manifest['jobs']):
        raise RuntimeError("Some runs failed or did not finish; retained data cannot support an all-node gate")
    print(f"COMPLETE {out}", flush=True)


if __name__ == '__main__':
    main()
