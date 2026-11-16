#!/usr/bin/env python3
"""pathtable.py: the on-path chain (runs/d*path/dnn.txt) as a markdown table; last line: the slack list for model.py."""
import glob, os, re
D = os.path.dirname(os.path.abspath(__file__)); best = 608215054; rows = []
for f in glob.glob(D + '/runs/d*path/dnn.txt'):
    for line in open(f):
        m = re.match(r'pid \d+ depth (\d+) fixes(.*?) \| first lb (\S+) .*? iters (\d+) cpu (\S+) .*?\| cut 1 best (\S+) rounds \d+ iters (\d+) cpu (\S+) status (\S+) cp \d active \d+ stalled (\d)', line)
        if m:
            d, fx, lb, it, cpu, b, pit, pcpu, st, stalled = m.groups()
            tag = f.split('/')[-2]; r, v = tag.split('_')[1], tag.split('_')[2][:-4]
            rows.append((int(d), f'{r}:{v}', float(lb), float(b), int(pit), 'stalled' if stalled == '1' else 'limit', float(cpu) + float(pcpu), int(it) + int(pit)))
rows.sort()
print('| depth | child (m) | first solve: certified, short | phase end: certified | short | % of best known | phase iterations (end) | CPU-s | s/it |')
print('|---|---|---|---|---|---|---|---|---|')
for d, fl, lb, b, pit, st, cpu, tot in rows:
    print(f'| {d} | `{fl}` ({60 - d}) | {lb:,.0f}, {best - lb:,.0f} | {b:,.0f} | {best - b:,.0f} | {100 * (best - b) / best:.3f}% | {pit} ({st}) | {cpu:,.0f} | {cpu / tot:.2f} |')
print('max certified on the path:', f'{max(r[3] for r in rows):,.0f}', '<=', f'{best:,}')
print('SLACK', ','.join(f'{(best - r[3]) / 1e6:.3f}' for r in rows))
