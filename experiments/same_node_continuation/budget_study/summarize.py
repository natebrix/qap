"""Apply the predeclared bound, fathom/iteration and whole-process CPU gates."""
import argparse
import json
import math
from pathlib import Path
from run import NODES, ARMS, CUTOFF


def evaluate(manifest):
    jobs = manifest['jobs']
    if len(jobs) != 30 or any(j['state'] != 'completed' for j in jobs):
        raise ValueError('All 30 runs must complete before declaring an all-node winner')
    lookup = {(j['node'], j['arm']): j for j in jobs}
    assert len(lookup) == 30
    rows = []
    for node in NODES:
        b = lookup[node['id'], 'B']
        for arm in ARMS:
            j = lookup[node['id'], arm]
            assert j['iterations'] <= j['max_iterations']
            assert all(math.isfinite(v) for v in j['final'].values())
            assert j['final']['first_lb'] == b['final']['first_lb'], 'Uncut controls differ'
            ratio = j['process_cpu_seconds']/b['process_cpu_seconds']
            diff = j['final']['best_lb']-b['final']['best_lb']
            valid = not (j['validity_question'] or b['validity_question'])
            ambiguous = abs(diff+2) <= 0.000002
            gate1 = valid and not ambiguous and diff >= -2
            gate2 = valid and j['outcome'] == b['outcome']
            if node['near']:
                gate2 = gate2 and b['outcome'] == 'fathomed' and j['first_crossing_iteration'] is not None and j['first_crossing_iteration'] < b['first_crossing_iteration']
            gate3 = j['process_cpu_seconds'] < b['process_cpu_seconds']
            row = dict(node=node['id'], arm=arm, best_polyak_round=j['best_polyak_round'],
                       best_polyak_bound=j['best_polyak_bound'], quiet_bound=j['quiet_bound'],
                       final_bound=j['final']['best_lb'], bound_difference=diff,
                       outcome=j['outcome'], validity_question=j['validity_question'],
                       first_crossing_iteration=j['first_crossing_iteration'],
                       iterations=j['iterations'], process_cpu_seconds=j['process_cpu_seconds'],
                       cpu_ratio=ratio, gate1=gate1, gate2=gate2, gate3=gate3,
                       rounding_ambiguous=ambiguous, reference=(arm=='B'),
                       eligible_cpu_fraction=ratio if gate1 and gate2 and arm!='B' else None)
            if arm == 'B':
                row.update(gate1=None, gate2=None, gate3=None)
            rows.append(row)
    aggregates = []
    total_b = sum(lookup[n['id'],'B']['process_cpu_seconds'] for n in NODES)
    for arm in ARMS:
        if arm == 'B':
            continue
        selected = [r for r in rows if r['arm']==arm]
        pass12 = [r for r in selected if r['gate1'] and r['gate2']]
        passall = all(r['gate1'] and r['gate2'] and r['gate3'] for r in selected)
        ratio = sum(r['process_cpu_seconds'] for r in selected)/total_b
        covered_b = sum(lookup[r['node'],'B']['process_cpu_seconds'] for r in pass12)
        aggregates.append(dict(arm=arm, passing_gate12_nodes=[r['node'] for r in pass12],
                               passing_all_nodes=[r['node'] for r in selected if r['gate1'] and r['gate2'] and r['gate3']],
                               passes_all_six=passall, total_cpu_ratio=ratio,
                               eligible_baseline_cpu_coverage=covered_b/total_b,
                               eligible_subset_cpu_fraction=sum(r['process_cpu_seconds'] for r in pass12)/covered_b if covered_b else None,
                               phase_cost_factor=(1/ratio if passall else None),
                               conditional_whole_run_factor=(1/(0.54+0.46*ratio) if passall else None)))
    return dict(source_commit=manifest['source_commit'], elapsed_wall_seconds=manifest['elapsed_wall_seconds'],
                rows=rows, arms=aggregates,
                shortened_winners=[a['arm'] for a in aggregates if a['arm']!='Y12' and a['passes_all_six']],
                validity_questions=[dict(node=r['node'], arm=r['arm']) for r in rows if r['validity_question']])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    data = evaluate(json.loads((args.directory/'manifest.json').read_text()))
    (args.directory/'summary.json').write_text(json.dumps(data, indent=2)+'\n')
    print('| Node | Arm | Best Polyak round | Quiet-only bound | Retained bound | Outcome | First crossing iteration | Total iterations | Process CPU s | CPU/B | Gates 1/2/3 | Eligible CPU/B |')
    print('|---|---|---:|---:|---:|---|---:|---:|---:|---:|---|---:|')
    for r in data['rows']:
        quiet='—' if r['quiet_bound'] is None else f"{r['quiet_bound']:.4f}"
        first='—' if r['first_crossing_iteration'] is None else str(r['first_crossing_iteration'])
        gates='reference' if r['reference'] else '/'.join('pass' if r[g] else 'FAIL' for g in ('gate1','gate2','gate3'))
        eligible='—' if r['eligible_cpu_fraction'] is None else f"{r['eligible_cpu_fraction']:.4f}"
        print(f"| {r['node']} | {r['arm']} | {r['best_polyak_round']} | {quiet} | {r['final_bound']:.6f} | {r['outcome']} | {first} | {r['iterations']} | {r['process_cpu_seconds']:.3f} | {r['cpu_ratio']:.4f} | {gates} | {eligible} |")
    print('\nShortened budgets passing all three gates on every node:', data['shortened_winners'])
    if data['validity_questions']:
        print('VALIDITY QUESTIONS:', data['validity_questions'])


if __name__ == '__main__':
    main()
