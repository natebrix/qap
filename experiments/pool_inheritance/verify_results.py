"""Audit the complete study against its prelaunch plan and retained raw logs."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess

from run import ROOT, HERE, POOL, digest, parse
from summarize import evaluate, prefix_cpu_interval


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    out = args.directory.resolve()
    manifest = json.loads((out/'manifest.json').read_text())
    plan = json.loads((HERE/'plan.json').read_text())
    assert not manifest.get('interrupted')
    for path, key in ((HERE/'plan.json', 'plan_sha256'),
                      (ROOT/'qaplib/wil50.dat', 'instance_sha256'),
                      (out/'parent.pool', 'parent_pool_sha256')):
        assert digest(path) == manifest[key], (path, key)
    binary = ROOT/'experiments/dnn_cuts/dnncuts'
    if binary.exists():
        assert digest(binary) == manifest['binary_sha256'], 'Different executable; audit the pinned binary'
    source = manifest['source_commit']
    for path in ('latest/dnncut.h', 'latest/dnncut.C', 'experiments/dnn_cuts/dnncuts.C',
                 'experiments/pool_inheritance/plan.json'):
        pinned = subprocess.check_output(['git', 'show', source+':'+path], cwd=ROOT)
        assert pinned == (ROOT/path).read_bytes(), 'Measured source changed: '+path
    jobs = manifest['jobs']
    assert len(jobs) == len(plan['jobs']) == 19
    by_id = {(j['node'], j['arm']): j for j in jobs}
    for job, specification in zip(jobs, plan['jobs']):
        assert (job['node'], job['arm']) == (specification['node'], specification['arm'])
        assert job['command'] == [str(out/'parent.pool') if x == POOL else x
                                  for x in specification['command']]
        raw = parse(out/job['log'])
        for key, value in raw.items():
            assert job[key] == value, (job['log'], key)
        final = re.findall(r'^FINAL best_lb (\S+) first_lb (\S+) gain (\S+)$',
                           (out/job['log']).read_text(), flags=re.MULTILINE)
        assert len(final) == 1 and tuple(map(float, final[0])) == tuple(
            job['final'][k] for k in ('best_lb', 'first_lb', 'gain'))
        assert math.isfinite(job['process_cpu_seconds']) and job['process_cpu_seconds'] > 0
        for i in range(len(job['rounds'])):
            prefix_cpu_interval(job, i)
    active = {}; started = set(); finished = set(); peak = 0
    for line in (out/'events.jsonl').read_text().splitlines():
        event = json.loads(line); key = (event['node'], event['arm'])
        if event['event'] == 'started':
            cpu = event['cpu']
            assert key not in started and cpu not in active and cpu in manifest['cores']
            assert by_id[key]['cpu'] == cpu
            if by_id[key]['arm'] in ('B', 'C'):
                assert ('n10', 'P') in finished
            active[cpu] = key; started.add(key); peak = max(peak, len(active))
        elif event['event'] == 'finished':
            assert key not in finished and event['exit_code'] == 0
            assert active.pop(by_id[key]['cpu']) == key
            assert event['process_cpu_seconds'] == by_id[key]['process_cpu_seconds']
            finished.add(key)
    assert not active and started == finished == set(by_id) and peak <= 8
    result = evaluate(manifest, plan)
    print(f'PASS: 19 single-attempt runs; pinned plan/source/input/pool; raw-log agreement; '
          f'parent dependency; distinct-core scheduling (peak {peak}); CPU accounting intervals.')
    print(f'Numerical gates: {result["decision"]}; median {result["median"]:.6f}; '
          f'd2 {result["d2_frozen_bound"]:.6f}; '
          f'validity questions {result["validity_questions"]}.')


if __name__ == '__main__':
    main()
