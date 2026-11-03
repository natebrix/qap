"""Predeclared parent-pool transfer study; Linux resource accounting and affinity."""
import argparse
from datetime import datetime, timezone
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / 'same_node_continuation'))
from run_cut_phase_ab import parse_log

BASE = '0be4e7d'  # expanded in plan.json before launch
CUTOFF = 48814.05
PARENT = '39:27,30:48,25:49'
SAMPLE_SEED = 20261103
QUEUE_SEED = 20261104
DEFAULT_OUT = 'experiments/pool_inheritance/results/runs'
POOL = DEFAULT_OUT + '/parent.pool'
COMMON = ['-grid', '5:10', '-side', '0', '-u', '48816', '-T', '1',
          '-v', '-delta', '50', '-nocap']


def plan():
    all32 = (ROOT/'experiments/wil50_feasibility/tuning/giveup/all32.txt').read_text().splitlines()
    depth4 = (ROOT/'experiments/wil50_feasibility/stronger_bounds/depth4/nodes.txt').read_text().splitlines()
    assert all32[9].strip() == PARENT
    for row, facility in ((1, 7), (2, 26), (4, 47)):
        assert depth4[row-1].split()[2] == PARENT+f',{facility}:5'
    candidates = [a for a in range(50) if a not in (39, 30, 25, 7, 26, 47)]
    selected = random.Random(SAMPLE_SEED).sample(candidates, 3)
    nodes = [dict(id=name, facility=a, fixes=PARENT+f',{a}:5', target=target)
             for name, a, target in [('d1', 7, CUTOFF), ('d2', 26, 48789.9),
                                     ('d4', 47, CUTOFF)] +
             [(f's{a}', a, None) for a in selected]]
    def command(fixes, options):
        return ['experiments/dnn_cuts/dnncuts', *COMMON, '-x', fixes,
                *options, 'qaplib/wil50.dat']
    phase = ['-I', '400', '-r', '12', '-i', '300', '-k', '1000', '-freeze', '0:1500']
    jobs = [dict(node='n10', arm='P', fixes=PARENT, dependency=None, max_iterations=5500,
                 command=command(PARENT, phase+['-savepool', POOL]))]
    rng = random.Random(QUEUE_SEED)
    order = nodes.copy(); rng.shuffle(order)
    for node in order:
        jobs.append(dict(node=node['id'], arm='A', fixes=node['fixes'], dependency=None,
                         max_iterations=1500, command=command(node['fixes'],
                         ['-I', '1500', '-r', '0', '-i', '300', '-k', '0'])))
    order = nodes.copy(); rng.shuffle(order)
    for node in order:
        for arm in ('B', 'C'):
            options = ['-I', '0', '-r', '0', '-i', '300', '-k', '0', '-freeze', '0:1500'] if arm == 'B' else phase
            jobs.append(dict(node=node['id'], arm=arm, fixes=node['fixes'], dependency='n10/P',
                             max_iterations=1500 if arm == 'B' else 5500,
                             command=command(node['fixes'], options+['-loadpool', POOL])))
    return dict(base_cognac=subprocess.check_output(['git', 'rev-parse', BASE], cwd=ROOT, text=True).strip(),
                cutoff=CUTOFF, sample_seed=SAMPLE_SEED, queue_seed=QUEUE_SEED,
                sampled_facilities=selected, nodes=nodes, jobs=jobs, max_workers=8,
                wall_limit_seconds=21600, warmY=False, parent_Z_transfer=False)


def parse(path):
    data = parse_log(path)
    total = 0
    for row in data['rounds']:
        total += row['iterations']; row['cumulative_iterations'] = total
    data['first_fathom_iteration'] = next((r['cumulative_iterations'] for r in data['rounds'] if r['status']=='fathomed'), None)
    text = path.read_text(errors='replace') if path.exists() else ''
    primal = re.findall(r'^(round\s+\d+|quiet round).*?primal (\S+)', text, flags=re.MULTILINE)
    for row, (_, value) in zip(data['rounds'], primal): row['primal'] = float(value)
    loaded = re.search(r'POOL loaded kept (\d+) dropped (\d+) source_round (\d+) source_quiet (\d+) source_value (\S+)', text)
    data['pool'] = None if not loaded else dict(kept=int(loaded[1]), dropped=int(loaded[2]),
                  source_round=int(loaded[3]), source_quiet=bool(int(loaded[4])), source_value=float(loaded[5]))
    saved = re.search(r'POOL saved cuts (\d+) source_round (\d+) source_quiet (\d+) source_value (\S+)', text)
    data['saved_pool'] = None if not saved else dict(cuts=int(saved[1]), source_round=int(saved[2]),
                        source_quiet=bool(int(saved[3])), source_value=float(saved[4]))
    return data


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path, data):
    def clean(x):
        if isinstance(x, float) and not math.isfinite(x): return repr(x)
        if isinstance(x, dict): return {k: clean(v) for k, v in x.items()}
        if isinstance(x, list): return [clean(v) for v in x]
        return x
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(clean(data), indent=2, allow_nan=False)+'\n'); temp.replace(path)


def physical_cores(count):
    selected, seen = [], set()
    for cpu in sorted(os.sched_getaffinity(0)):
        root = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        key = ((root/'physical_package_id').read_text(), (root/'core_id').read_text())
        if key not in seen: seen.add(key); selected.append(cpu)
    quota = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    assert len(selected) >= count
    assert quota[0] == 'max' or count <= int(quota[0])/int(quota[1])
    return selected[:count]


def process_cpu(pid):
    # Child /proc entries are unavailable in this executor. Record the round
    # without a CPU observation; wait4 supplies exact terminal accounting.
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', action='store_true')
    ap.add_argument('--out', type=Path, default=ROOT/DEFAULT_OUT)
    ap.add_argument('--workers', type=int, default=8, choices=range(1,9))
    args = ap.parse_args(); definition = plan()
    if args.plan:
        atomic_json(HERE/'plan.json', definition)
        print([(j['node'], j['arm']) for j in definition['jobs']]); return
    assert json.loads((HERE/'plan.json').read_text()) == definition, 'Predeclared plan changed'
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    git = lambda *a: subprocess.check_output(['git', *a], cwd=ROOT, text=True).strip()
    assert not git('status', '--porcelain'), 'Commit protocol and source before launch'
    cpus = physical_cores(args.workers)
    started = time.monotonic()
    manifest = dict(source_commit=git('rev-parse', 'HEAD'), base_cognac=definition['base_cognac'],
        plan_sha256=digest(HERE/'plan.json'), binary_sha256=digest(ROOT/'experiments/dnn_cuts/dnncuts'),
        instance_sha256=digest(ROOT/'qaplib/wil50.dat'), started_utc=datetime.now(timezone.utc).isoformat(),
        workers=args.workers, cores=cpus, wall_limit_seconds=21600, jobs=[], parent_pool_sha256=None,
        environment={k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
    for spec in definition['jobs']:
        job = dict(spec, state='queued', log=f"{spec['node']}_{spec['arm']}.log", stage_observations=[])
        job['command'] = [str(out/'parent.pool') if x==POOL else x for x in job['command']]
        manifest['jobs'].append(job)
    events = (out/'events.jsonl').open('w', buffering=1); running = {}
    env = dict(os.environ, **manifest['environment'])
    def persist():
        manifest['elapsed_wall_seconds'] = time.monotonic()-started
        atomic_json(out/'manifest.json', manifest)
    def event(kind, job, **extra):
        events.write(json.dumps(dict(event=kind, elapsed=time.monotonic()-started,
                                     node=job['node'], arm=job['arm'], **extra))+'\n')
    def observe(job, result, cpu, at_exit=False):
        for row in result['rounds'][len(job['stage_observations']):]:
            observation = dict(row, observed_process_cpu_seconds=cpu, observed_at_exit=at_exit,
                               elapsed=time.monotonic()-started)
            job['stage_observations'].append(observation); event('round_observed', job, round=observation)
    def finish(pid, status, usage):
        proc, handle, job, t0 = running.pop(pid)
        proc.returncode = os.waitstatus_to_exitcode(status); handle.close()
        result = parse(out/job['log']); cpu = usage.ru_utime+usage.ru_stime
        observe(job, result, cpu, True)
        good = proc.returncode==0 and result['final'] is not None and not result['solver_error']
        if good: good = all(math.isfinite(x) for x in result['final'].values())
        if good and job['arm']=='P':
            good = result['saved_pool'] is not None and (out/'parent.pool').is_file()
            if good: manifest['parent_pool_sha256'] = digest(out/'parent.pool')
        if good and job['arm'] in ('B','C'):
            good = result['pool'] is not None and digest(out/'parent.pool')==manifest['parent_pool_sha256']
        job.pop('progress', None)
        job.update(result, state='completed' if good else 'failed', exit_code=proc.returncode,
                   process_cpu_seconds=cpu, wall_seconds=time.monotonic()-t0, max_rss_kib=usage.ru_maxrss)
        event('finished', job, exit_code=proc.returncode, process_cpu_seconds=cpu)
        print(f"DONE {job['node']} {job['arm']} {job['state']} {job['final']}", flush=True)
        if not good: raise RuntimeError('Failed run; preserve logs, no gate decision')
    try:
        last_report = 0
        while running or any(j['state']=='queued' for j in manifest['jobs']):
            elapsed = time.monotonic()-started
            if elapsed >= 21600: raise TimeoutError('Six-hour study deadline')
            completed = {j['node']+'/'+j['arm'] for j in manifest['jobs'] if j['state']=='completed'}
            used = {v[2]['cpu'] for v in running.values()}
            for cpu in cpus:
                if cpu in used: continue
                job = next((j for j in manifest['jobs'] if j['state']=='queued' and
                            (j['dependency'] is None or j['dependency'] in completed)), None)
                if job is None: break
                handle = (out/job['log']).open('w')
                proc = subprocess.Popen(job['command'], cwd=ROOT, env=env, stdout=handle,
                    stderr=subprocess.STDOUT, start_new_session=True,
                    preexec_fn=lambda cpu=cpu: os.sched_setaffinity(0, {cpu}))
                job.update(state='running', cpu=cpu, started_elapsed=elapsed)
                running[proc.pid] = (proc, handle, job, time.monotonic()); event('started', job, cpu=cpu)
                print(f"START {job['node']} {job['arm']} core={cpu}", flush=True)
            for pid in list(running):
                got, status, usage = os.wait4(pid, os.WNOHANG)
                if got: finish(pid, status, usage)
            for pid, (_, _, job, _) in running.items():
                result = parse(out/job['log'])
                if len(result['rounds']) > len(job['stage_observations']):
                    try: cpu = process_cpu(pid)
                    except (OSError, ValueError): cpu = None
                    observe(job, result, cpu); job['progress'] = result
            persist()
            if elapsed-last_report >= 60:
                print(f"PROGRESS minutes={elapsed/60:.1f} done={sum(j['state']=='completed' for j in manifest['jobs'])}/19 running={len(running)}", flush=True)
                last_report = elapsed
            if running: time.sleep(5)
    finally:
        for pid in list(running):
            os.killpg(pid, signal.SIGKILL)
            _, status, usage = os.wait4(pid, 0)
            try: finish(pid, status, usage)
            except RuntimeError: pass
            manifest['interrupted'] = True
        persist(); events.close()
    if any(j['state']!='completed' for j in manifest['jobs']): raise RuntimeError('Incomplete study')
    print(f'COMPLETE {out}', flush=True)


if __name__=='__main__':
    main()
