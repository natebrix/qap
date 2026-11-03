"""Complete only the four executor-interrupted C attempts, without changing inputs."""
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time

from run import HERE, ROOT, atomic_json, digest, parse, physical_cores, process_cpu


def main():
    definition = json.loads((HERE/'recovery_plan.json').read_text())
    out = HERE/'results/runs'
    original = HERE/definition['original_snapshot']
    manifest = json.loads((original/'manifest.json').read_text())
    assert digest(original/'manifest.json') == definition['original_manifest_sha256']
    assert digest(ROOT/'experiments/dnn_cuts/dnncuts') == definition['binary_sha256']
    assert digest(out/'parent.pool') == definition['pool_sha256']
    assert all(j['cpu'] in physical_cores(8) for j in definition['jobs'])
    git = lambda *args: subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()
    assert not git('status', '--porcelain'), 'Publish the recovery plan before launch'
    epoch = datetime.fromisoformat(manifest['started_utc'])
    offset = (datetime.now(timezone.utc)-epoch).total_seconds()
    start = time.monotonic()
    elapsed = lambda: offset+time.monotonic()-start
    remaining = definition['wall_limit_seconds']-offset
    assert remaining > 0
    manifest['recovery'] = dict(plan_sha256=digest(HERE/'recovery_plan.json'),
        source_commit=git('rev-parse', 'HEAD'), started_utc=datetime.now(timezone.utc).isoformat(),
        original_snapshot=definition['original_snapshot'],
        interrupted_whole_process_CPU_unknown=True, remaining_wall_seconds=remaining)
    selected = {(x['node'], x['arm']) for x in definition['jobs']}
    jobs = [j for j in manifest['jobs'] if (j['node'], j['arm']) in selected]
    assert len(jobs) == 4
    for j in jobs:
        assert j['state'] == 'running'
        j['interrupted_first_started_elapsed'] = j['started_elapsed']
        j.update(state='queued', attempt=2, stage_observations=[])
        j.pop('progress', None)
    events = (out/'recovery_events.jsonl').open('x', buffering=1)
    running = {}
    env = dict(os.environ, **manifest['environment'])
    def persist():
        manifest['elapsed_wall_seconds'] = elapsed()
        atomic_json(out/'manifest.json', manifest)
    def event(kind, job, **kw):
        events.write(json.dumps(dict(event=kind, elapsed=elapsed(), node=job['node'],
                                    arm=job['arm'], attempt=2, **kw))+'\n')
    def observe(job, result, cpu, at_exit=False):
        for row in result['rounds'][len(job['stage_observations']):]:
            item = dict(row, observed_process_cpu_seconds=cpu, observed_at_exit=at_exit,
                        elapsed=elapsed())
            job['stage_observations'].append(item)
            event('round_observed', job, round=item)
    def finish(pid, status, usage):
        proc, handle, job, t0 = running.pop(pid)
        proc.returncode = os.waitstatus_to_exitcode(status)
        handle.close()
        result = parse(out/job['log']); cpu = usage.ru_utime+usage.ru_stime
        observe(job, result, cpu, True)
        previous = parse(original/job['log'])['rounds']
        fields = ('round', 'lb', 'primal', 'iterations', 'status')
        matches = [[r[k] for k in fields] for r in result['rounds'][:len(previous)]] == [
                  [r[k] for k in fields] for r in previous]
        good = proc.returncode == 0 and result['final'] is not None and not result['solver_error']
        if good:
            good = all(math.isfinite(x) for x in result['final'].values())
        assert digest(out/'parent.pool') == definition['pool_sha256']
        job.pop('progress', None)
        job.update(result, state='completed' if good else 'failed', exit_code=proc.returncode,
            process_cpu_seconds=cpu, wall_seconds=time.monotonic()-t0, max_rss_kib=usage.ru_maxrss,
            interrupted_prefix_matches=matches, interrupted_prefix_rounds=len(previous))
        event('finished', job, exit_code=proc.returncode, process_cpu_seconds=cpu,
              interrupted_prefix_matches=matches)
        print('DONE', job['node'], job['arm'], job['state'], job['final'], 'prefix_match', matches, flush=True)
        if not good:
            raise RuntimeError('Recovery run failed; preserve all attempts')
    try:
        for job in jobs:
            handle = (out/job['log']).open('w')
            cpu = job['cpu']
            proc = subprocess.Popen(job['command'], cwd=ROOT, env=env, stdout=handle,
                stderr=subprocess.STDOUT, start_new_session=True,
                preexec_fn=lambda cpu=cpu: os.sched_setaffinity(0, {cpu}))
            job.update(state='running', started_elapsed=elapsed())
            running[proc.pid] = (proc, handle, job, time.monotonic())
            event('started', job, cpu=cpu)
            print('START', job['node'], job['arm'], 'attempt 2 core', cpu, flush=True)
        persist(); last_report = 0
        while running:
            if elapsed() >= definition['wall_limit_seconds']:
                raise TimeoutError('Original six-hour study deadline')
            for pid in list(running):
                got, status, usage = os.wait4(pid, os.WNOHANG)
                if got:
                    finish(pid, status, usage)
            for pid, (_, _, job, _) in running.items():
                result = parse(out/job['log'])
                if len(result['rounds']) > len(job['stage_observations']):
                    try:
                        cpu = process_cpu(pid)
                    except (OSError, ValueError):
                        cpu = None
                    observe(job, result, cpu)
            persist()
            if elapsed()-last_report >= 60:
                print('PROGRESS total_minutes', round(elapsed()/60, 1), 'recovery_running', len(running), flush=True)
                last_report = elapsed()
            if running:
                time.sleep(5)
    finally:
        for pid in list(running):
            os.killpg(pid, signal.SIGKILL)
            _, status, usage = os.wait4(pid, 0)
            try:
                finish(pid, status, usage)
            except RuntimeError:
                pass
            manifest['interrupted'] = True
        persist(); events.close()
    assert all(j['state'] == 'completed' for j in manifest['jobs'])
    print('COMPLETE recovery; 19 completed runs, 4 retained interrupted attempts', flush=True)


if __name__ == '__main__':
    main()
