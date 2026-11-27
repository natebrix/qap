#!/usr/bin/env python3
"""summary.py RUNS NODES: probe.sh's firings beside sko42-261120b (NODES, its node trace): at each depth-1 node the line
rule 9 chose there and the subtree it grew, the line strong branching chooses, where rule 9's line ranks among the
candidates, and what the firing cost (child solves, under the node's pool or not, retried, still failed, CPU s)."""
import gzip, os, re, sys
runs, nodes = sys.argv[1], sys.argv[2]
prod, sub = {}, {}
with gzip.open(nodes, 'rt') as f:
    for ln in f:
        x = ln.split()
        if len(x) < 12 or x[0] != 'node' or x[1] == '0':
            continue
        k = x[2].split(',')[0]
        s = sub.setdefault(k, [0, 0.0])
        s[0] += 1
        s[1] += float(x[11])
        if x[1] == '1':
            prod[k] = x[9]
print('| node | rule 9 | its subtree: nodes, CPU h | strong branching | rule 9\'s line: rank of candidates, score / best | '
      'the runner-up / best | child solves | pooled | retried | failed | firing CPU s |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for d in sorted(os.listdir(runs), key=lambda d: int(d.split('_')[1]) if d.startswith('36_') else -1):
    if not d.startswith('36_'):
        continue
    k = '36:' + d.split('_')[1]
    t = os.path.join(runs, d, 'sb.txt')
    line = open(t).read().strip().splitlines()[-1] if os.path.exists(t) and os.path.getsize(t) else ''
    n, c = sub.get(k, [0, 0.0])
    if not line:
        p = os.path.join(runs, d, 'progress.txt')
        done = sum(1 for _ in open(p)) if os.path.exists(p) else 0
        print(f'| {k} | {prod.get(k, "?")} | {n:,}, {c / 3600:.1f} | (not finished: {done} child solves so far) | | | | | | | |')
        continue
    m = re.search(r'chose (\S+) score (\S+) open \d+ \| (\d+) solves( \(pooled\))?, (\d+) reused, (\d+) retried, (\d+) failed, '
                  r'\d+ iterations, ([\d.]+) CPU s', line)
    # the firing's whole cost from its progress file: a firing resumed after a stop counts only its new solves in
    # the trace line (the rest as reused); the progress file has every child solve, retried or not
    solves = cpu = retried = fathomed = 0
    for pl in open(os.path.join(runs, d, 'progress.txt')):
        f = pl.split()
        solves += 1
        cpu += float(f[f.index('cpu') + 1])
        retried += f[-2] == 'retry' and f[-1] != '0'
        fathomed += f[f.index('fathomed') + 1] == '1'
    cands = [c.split(':') for c in line.split('|')[-1].split()]
    scores = {c[0]: float(c[2]) for c in cands if len(c) == 4}
    order = sorted(scores, key=lambda l: scores[l])
    r9 = prod.get(k, '?')
    rank = f'{order.index(r9) + 1} of {len(order)}, {scores[r9] / max(scores[order[0]], 1e-12):.2f}x' if r9 in scores else 'not a candidate'
    second = f'{scores[order[1]] / max(scores[order[0]], 1e-12):.3f}x' if len(order) > 1 else '-'
    print(f'| {k} | {r9} | {n:,}, {c / 3600:.1f} | {m.group(1)}{" (the same)" if m.group(1) == r9 else ""} | {rank} | '
          f'{second} | {solves} ({fathomed} fathomed) | {"yes" if m.group(4) else "no"} | {retried} | {m.group(7)} | {cpu:,.0f} |')
