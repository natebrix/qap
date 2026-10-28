"""Compute certified DNN lower bounds on nug20 open nodes.
usage: python3 batch.py depth2|depth3 [nsample] [maxit]
"""
import sys, json, time
import numpy as np
from collections import defaultdict
from qaputil import *
from dnn import dnn_admm

D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn/'
OPT = 2570.0
FATHOM = OPT - 1.95

n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug20.dat')
which = sys.argv[1]
nsample = int(sys.argv[2]) if len(sys.argv) > 2 else 0
maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 3000

if which == 'depth2':
    _, opened = read_state_nodes(D + 'n20_d1.state')     # the 98 depth-2 nodes as created
    _, kids = read_state_nodes(D + 'n20_d2.state')       # their surviving children (carry parent's QPB)
else:
    _, opened = read_state_nodes(D + 'n20_d2.state')     # the 1474 depth-3 nodes as created
    _, kids = read_state_nodes(D + 'n20_d3.state')
own_qpb = {}
for k in kids:
    key = parent_key(k['p'], k['hist'], k['depth'])
    own_qpb[key] = k['bound']

rows = []
for x in opened:
    key = tuple(sorted((i, l) for i, l in enumerate(x['p']) if l >= 0))
    rows.append(dict(key=key, p=x['p'], inherited=x['bound'], pred=x['pred'],
                     qpb=own_qpb.get(key), branched=key in own_qpb))
print('opened', len(rows), 'branched (QPB failed to fathom)', sum(r['branched'] for r in rows), flush=True)

if which == 'depth3':
    rng = np.random.default_rng(0)
    cand = [r for r in rows if r['branched']]
    idx = rng.choice(len(cand), nsample, replace=False)
    rows = [cand[i] for i in sorted(idx)]

import os
outfn = D + 'results_%s.jsonl' % which
done = set()
if os.path.exists(outfn):
    for l in open(outfn):
        try:
            done.add(tuple(tuple(t) for t in json.loads(l)['key']))
        except Exception:
            pass
print('already done:', len(done), flush=True)
out = open(outfn, 'a')
rng = np.random.default_rng(123)
for r in rows:
    if tuple(r['key']) in done:
        continue
    p = r['p']
    verify_reduction(a, b, c, p, rng)
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    res = dnn_admm(AFF, BLL, Cp, const, maxit=maxit, check_every=50, fathom_at=FATHOM,
                   time_limit=150, stall_window=6)
    rec = dict(key=[list(t) for t in r['key']], m=len(F), inherited=r['inherited'], pred=r['pred'],
               qpb=r['qpb'], branched=r['branched'], lb=res['lb'], lb_last=res['lb_last'], obj=res['obj'],
               it=res['it'], cpu=res['cpu'], status=res['status'],
               hist=[(h[0], round(h[1], 3), round(h[2], 3)) for h in res['hist']])
    out.write(json.dumps(rec) + '\n'); out.flush()
    print('%s m=%d qpb=%s lb=%.2f obj=%.2f it=%d cpu=%.1f %s fathom=%s' % (
        rec['key'], rec['m'], None if r['qpb'] is None else round(r['qpb'], 1), res['lb'], res['obj'],
        res['it'], res['cpu'], res['status'], res['lb'] > FATHOM), flush=True)
