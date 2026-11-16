#!/usr/bin/env python3
"""bbsum.py prefix arm...: per arm of a bnbdrv run traced with QAP_DNN_TRACE (logs/<prefix>_<arm>.dnn.txt):
DNN nodes, fathomed, give-ups, iterations (fathomed / not), DNN CPU."""
import sys, re
pre = sys.argv[1]
print('| arm | finished | DNN nodes | fathomed | not fathomed | iterations: fathomed nodes | not fathomed | total | DNN CPU s |')
print('|---|---|---|---|---|---|---|---|---|')
for arm in sys.argv[2:]:
    n = f = fi = gi = 0; cpu = 0.0
    for l in open('logs/%s_%s.dnn.txt' % (pre, arm)):
        m = re.search(r'status (\w+) iters (\d+) cpu ([\d.]+)', l)
        n += 1; i = int(m.group(2)); cpu += float(m.group(3))
        if m.group(1) == 'fathomed': f += 1; fi += i
        else: gi += i
    ex = [l.strip() for l in open('logs/%s_%s.log' % (pre, arm)) if l.startswith('exit')]
    print('| %s | %s | %d | %d | %d | %d | %d | %d | %.0f |' % (arm, 'yes' if ex == ['exit 0'] else 'no (%s)' % ex, n, f, n - f, fi, gi, fi + gi, cpu))
