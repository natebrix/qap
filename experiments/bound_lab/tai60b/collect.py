#!/usr/bin/env python3
"""collect.py [PREFIX...]: one markdown row per node from runs/TAG/dnn.txt (QAP_DNN_TRACE) for the TAGs
starting with a PREFIX (default: all): fixes, first solve (shortfall, iterations, status), phase
(iterations, end, status), outcome, total iterations, CPU-s, s/iteration."""
import glob, os, re, sys
D = os.path.dirname(os.path.abspath(__file__))
pre = sys.argv[1:] or ['']
rows = []
for f in sorted(glob.glob(D + '/runs/*/dnn.txt')):
    tag = f.split('/')[-2]
    if not any(tag.startswith(p) for p in pre): continue
    for line in open(f):
        m = re.match(r'pid \d+ depth (\d+) fixes(.*?) \| first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+) warm (\d) pool (\d) kept (\d+) dropped (\d+)(?: pcp 1 handed \d+)? \| cut (\d)(.*?) \| cutoff (\S+) fathomed (\d)', line)
        if not m: continue
        depth, fixes, lb, st, it, cpu, warm, pool, kept, dropped, cut, crest, cutoff, fath = m.groups()
        lb, cutoff, it, cpu = float(lb), float(cutoff), int(it), float(cpu)
        pit = 0; pcpu = 0.0; pend = ''; pst = ''
        if cut == '1':
            c = re.search(r'best (\S+) rounds \d+ iters (\d+) cpu (\S+) status (\S+)(?: cp \d active \d+ stalled (\d) gaveup (\d))?', crest)
            best, pit, pcpu, pst = float(c.group(1)), int(c.group(2)), float(c.group(3)), c.group(4)
            if c.group(5) == '1': pst = 'stalled'
            if c.group(6) == '1': pst = 'gaveup'
            pend = f'{cutoff - best:,.0f}' if pst != 'fathomed' else 'fathomed'
        out = 'fathomed (first solve)' if fath == '1' and cut == '0' else ('fathomed (phase)' if fath == '1' else 'branched')
        tot = it + pit
        rows.append((tag, depth, fixes.strip().replace(' ', ','), f'{cutoff - lb:,.0f}', it, st, pit, pend, pst, out, tot, cpu + pcpu, (cpu + pcpu) / max(1, tot)))
print('| run | depth | fixes | first solve: short | its | status | phase its | phase end: short | phase status | outcome | total its | CPU-s | s/it |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    print('| ' + ' | '.join(str(x) if not isinstance(x, float) else f'{x:.2f}' if x < 10 else f'{x:,.0f}' for x in r) + ' |')
