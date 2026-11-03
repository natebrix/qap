"""Check log interpretation and the predeclared gate contract before launch."""
from copy import deepcopy
from pathlib import Path
from run import ROOT, NODES, ARMS, parse_run, plan
from summarize import evaluate

assert len(plan()['jobs']) == 30
path = ROOT / 'experiments/wil50_feasibility/stronger_bounds/dual/nocap_6.log'
near = parse_run(path, True)
assert near['first_crossing_iteration'] == 1260 and near['outcome'] == 'fathomed'
assert near['quiet_bound'] is None and near['best_polyak_round'] == 3
unexpected = parse_run(path, False)
assert unexpected['outcome'] == 'validity_question' and unexpected['validity_question']
miss = parse_run(ROOT / 'experiments/wil50_feasibility/stronger_bounds/dual/tri_13.log', False)
assert miss['iterations'] == 5500 and miss['first_crossing_iteration'] is None
assert miss['best_polyak_round'] == 12 and miss['quiet_bound'] == 48755.4220

jobs = []
for n in NODES:
    for arm in ARMS:
        baseline = arm == 'B'
        value = 48815.0 if n['near'] else 48750.0
        jobs.append(dict(node=n['id'], arm=arm, state='completed', max_iterations=5500,
                         iterations=1200 if baseline else 1000,
                         first_crossing_iteration=(1200 if baseline else 1000) if n['near'] else None,
                         final=dict(best_lb=value, first_lb=48000.0, gain=value-48000),
                         outcome='fathomed' if n['near'] else 'short', validity_question=False,
                         process_cpu_seconds=100.0 if baseline else 70.0,
                         best_polyak_round=3, best_polyak_bound=value,
                         quiet_bound=None if n['near'] else value))
m=dict(jobs=jobs, source_commit='synthetic contract fixture', elapsed_wall_seconds=1)
s=evaluate(m)
assert s['shortened_winners'] == ['Y6','Y8s','Y12s']
assert all(abs(a['total_cpu_ratio']-0.7)<1e-12 for a in s['arms'])

for gate in (1, 2, 3, 'validity'):
    changed=deepcopy(m)
    row=next(j for j in changed['jobs'] if j['node']=='n14' and j['arm']=='Y6')
    if gate==1: row['final']['best_lb']-=3
    elif gate==2: row['first_crossing_iteration']=1200
    elif gate==3: row['process_cpu_seconds']=101
    else: row['validity_question']=True; row['outcome']='validity_question'
    result=evaluate(changed)
    assert 'Y6' not in result['shortened_winners']
    match=next(r for r in result['rows'] if r['node']=='n14' and r['arm']=='Y6')
    if gate in (1,2,'validity'): assert match['eligible_cpu_fraction'] is None

changed=deepcopy(m);changed['jobs'][0]['state']='failed'
try:
    evaluate(changed)
    raise AssertionError('Incomplete data accepted')
except ValueError:
    pass
print('PASS plan=30 historic_crossing=1260 quiet_parse=5500 validity_classification=pass gates_1_2_3=pass incomplete_rejection=pass')
