"""Apply the declared median-gain/d2 gates after all nineteen runs complete."""
import argparse
import json
import math
from pathlib import Path
import statistics
from run import CUTOFF, digest


def gate(gains, d2_bound):
    if len(gains)!=6 or not all(math.isfinite(x) for x in gains+[d2_bound]):
        raise ValueError('Six finite paired gains and a finite d2 bound required')
    median = statistics.median(gains)
    ambiguous = abs(median-40)<0.000002 or abs(d2_bound-48759.9)<0.000002
    return dict(median=median, minimum=min(gains), maximum=max(gains), sorted_gains=sorted(gains),
                median_pass=median>=40, d2_pass=d2_bound>=48759.9,
                decision='ROUNDING_AMBIGUOUS' if ambiguous else
                         'KILL' if median<40 or d2_bound<48759.9 else 'PASS')


def prefix_cpu_interval(job, index):
    # Each DNN call reports disjoint main-thread CPU to 0.1 s. All other
    # process CPU is nonnegative overhead whose placement is unknown. Charge
    # none/all of it before the boundary, with printing/accounting slack.
    rows=job['rounds'];total=job['process_cpu_seconds']
    before=sum(r['dnn_cpu'] for r in rows[:index+1])
    after=sum(r['dnn_cpu'] for r in rows[index+1:])
    lower=max(0.0,before-0.05*(index+1)-0.01)
    upper=min(total,total-after+0.05*(len(rows)-index-1)+0.01)
    if lower>upper+0.01:raise ValueError('Inconsistent round/whole-process CPU accounting')
    return dict(lower_seconds=lower,upper_seconds=max(lower,upper))


def first_target(job, target):
    if target is None: return None
    for i, row in enumerate(job['rounds']):
        crossed = row['status']=='fathomed' if target==CUTOFF else row['lb']>target+0.00005
        if crossed:
            return dict(stage=row['round'], iterations=row['cumulative_iterations'],
                        reported_bound=row['lb'],
                        process_cpu_interval=prefix_cpu_interval(job,i),
                        terminal_process_cpu_seconds=job['process_cpu_seconds'] if row['status']=='fathomed' else None)
    return None


def evaluate(manifest, definition):
    jobs = manifest['jobs']
    if len(jobs)!=19 or any(j['state']!='completed' or j['exit_code']!=0 for j in jobs):
        raise ValueError('All nineteen declared runs must finish successfully')
    lookup={(j['node'],j['arm']):j for j in jobs};assert len(lookup)==19
    parent=lookup['n10','P'];source=parent['saved_pool'];assert source is not None
    rows=[];deltas=[];diagnostics=[]
    for node in definition['nodes']:
        pair={arm:lookup[node['id'],arm] for arm in ('A','B','C')}
        for arm,j in pair.items():
            assert j['iterations']<=j['max_iterations']
            assert all(math.isfinite(x) for x in j['final'].values())
            assert len(j['rounds'])==len(j['stage_observations'])
            if arm!='A':
                p=j['pool'];assert p['kept']+p['dropped']==source['cuts']
                assert p['source_value']==source['source_value']
        assert pair['B']['pool']==pair['C']['pool']
        assert len(pair['A']['rounds'])==1 and pair['A']['rounds'][0]['round']=='round  0'
        assert len(pair['B']['rounds'])==1 and pair['B']['rounds'][0]['round']=='quiet round'
        delta=pair['B']['final']['best_lb']-pair['A']['final']['best_lb'];deltas.append(delta)
        difference=pair['B']['final']['best_lb']-source['source_value']
        if difference < -0.000001: diagnostics.append(dict(node=node['id'],below_parent_by=-difference))
        for arm,j in pair.items():
            rows.append(dict(node=node['id'],fixes=node['fixes'],arm=arm,bound=j['final']['best_lb'],
                frozen_minus_uncut=delta if arm=='B' else None,
                frozen_minus_parent=difference if arm=='B' else None,
                pool=j['pool'],fathomed=j['fathomed'],iterations=j['iterations'],
                first_fathom_iteration=j['first_fathom_iteration'],
                process_cpu_seconds=j['process_cpu_seconds'],
                historical_target=node['target'],first_reported_target=first_target(j,node['target']) if arm=='C' else None))
    # Every declared child fixes the same opposite-side index (location 5).
    # Its region remap varies, but cut membership/drop decisions must agree.
    assert len({(lookup[node['id'],'B']['pool']['kept'],
                 lookup[node['id'],'B']['pool']['dropped']) for node in definition['nodes']})==1
    d2=lookup['d2','B']['final']['best_lb'];result=gate(deltas,d2)
    result.update(d2_frozen_bound=d2,rows=rows,parent=parent['final'],parent_pool=source,
                  parent_cpu_seconds=parent['process_cpu_seconds'],parent_floor_diagnostics=diagnostics,
                  validity_questions=[dict(node=j['node'],arm=j['arm']) for j in jobs if j['node'] in ('n10','d2') and j['fathomed']],
                  elapsed_wall_seconds=manifest['elapsed_wall_seconds'],
                  total_process_cpu_seconds=sum(j['process_cpu_seconds'] for j in jobs))
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('directory',type=Path);args=ap.parse_args()
    manifest=json.loads((args.directory/'manifest.json').read_text())
    assert digest(args.directory/'parent.pool')==manifest['parent_pool_sha256']
    definition=json.loads((Path(__file__).parent/'plan.json').read_text())
    result=evaluate(manifest,definition)
    (args.directory/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
    print('| Node | Arm | Certified bound | B-A | B-parent | Kept/dropped | Fathomed | Iterations | Process CPU s |')
    print('|---|---|---:|---:|---:|---|---|---:|---:|')
    for r in result['rows']:
        delta='—' if r['frozen_minus_uncut'] is None else f"{r['frozen_minus_uncut']:+.6f}"
        floor='—' if r['frozen_minus_parent'] is None else f"{r['frozen_minus_parent']:+.6f}"
        count='—' if r['pool'] is None else f"{r['pool']['kept']}/{r['pool']['dropped']}"
        print(f"| {r['node']} | {r['arm']} | {r['bound']:.6f} | {delta} | {floor} | {count} | {r['fathomed']} | {r['iterations']} | {r['process_cpu_seconds']:.3f} |")


if __name__=='__main__':main()
