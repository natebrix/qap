"""Summarize results_depth*.jsonl: counts, fathom rates, percentiles of QPB / lb / cpu / iterations."""
import json, glob, os
import numpy as np
D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn30/'
FATHOM = 6124 - 1.95
P = [0, 10, 25, 50, 75, 90, 100]
for fn in sorted(glob.glob(D + 'results_depth*.jsonl')):
    recs = [json.loads(l) for l in open(fn)]
    if not recs:
        continue
    print('==', os.path.basename(fn), 'nodes', len(recs), 'm', sorted(set(r['m'] for r in recs)))
    fath = [r for r in recs if r['lb'] > FATHOM]
    print('fathomed (certified lb > %.2f): %d / %d' % (FATHOM, len(fath), len(recs)))
    from collections import Counter
    print('status', Counter(r['status'] for r in recs))
    def pct(name, v):
        v = np.array(v, dtype=float)
        print('%-14s' % name, ' '.join('%8.1f' % x for x in np.percentile(v, P)), '  mean %.1f  sum %.1f' % (v.mean(), v.sum()))
    print('%-14s' % 'percentile', ' '.join('%8d' % x for x in P))
    pct('qpb', [r['qpb'] for r in recs if r['qpb'] is not None])
    pct('lb', [r['lb'] for r in recs])
    pct('lb-qpb', [r['lb'] - r['qpb'] for r in recs if r['qpb'] is not None])
    pct('6124-lb', [6124 - r['lb'] for r in recs])
    pct('obj(primal)', [r['obj'] for r in recs])
    pct('iters', [r['it'] for r in recs])
    pct('cpu s', [r['cpu'] for r in recs])
    pct('s/iter', [r['cpu'] / r['it'] for r in recs])
    if fath:
        pct('cpu s (fath)', [r['cpu'] for r in fath])
    nf = [r for r in recs if r['lb'] <= FATHOM]
    for r in nf:
        print('  NOT fathomed:', r['key'], 'qpb %.1f lb %.2f obj %.2f it %d cpu %.0f %s' % (r['qpb'], r['lb'], r['obj'], r['it'], r['cpu'], r['status']))
