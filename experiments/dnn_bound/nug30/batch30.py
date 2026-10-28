"""Certified DNN lower bounds on QPB-open nodes of nug30 (rule2 frontier).
usage: python3 batch30.py DEPTH NSAMPLE [MAXIT] [TIME_LIMIT_S] [SEED]
Nodes at depth d "as created" are read from n30_d{d-1}.state; their own QPB values come from the
children written in n30_d{d}.state (child.bound = parent's effective bound).  A uniform random
sample (seed) of the QPB-open nodes is processed in random order; results are appended per node to
results_depth{d}.jsonl and finished nodes are skipped on rerun.
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn')
from qaputil import *
from dnn import dnn_admm

D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn30/'
OPT = 6124.0
FATHOM = OPT - 1.95

depth = int(sys.argv[1])
nsample = int(sys.argv[2])
maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
tlim = float(sys.argv[4]) if len(sys.argv) > 4 else 900.0
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0

n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug30.dat')
_, opened = read_state_nodes(D + 'n30_d%d.state' % (depth - 1))
_, kids = read_state_nodes(D + 'n30_d%d.state' % depth)
own_qpb = {}
for k in kids:
    own_qpb[parent_key(k['p'], k['hist'], k['depth'])] = k['bound']
rows = []
for x in opened:
    key = tuple(sorted((i, l) for i, l in enumerate(x['p']) if l >= 0))
    rows.append(dict(key=key, p=x['p'], inherited=x['bound'], pred=x['pred'],
                     qpb=own_qpb.get(key), branched=key in own_qpb))
nopen = sum(r['branched'] for r in rows)
print('depth', depth, 'opened', len(rows), 'QPB-open (branched)', nopen, flush=True)
cand = [r for r in rows if r['branched']]
mode = sys.argv[6] if len(sys.argv) > 6 else 'random'
outsuffix = ''
if mode == 'low':        # the nsample QPB-open nodes with the LOWEST own QPB (tail where DNN misses sit)
    rows = sorted(cand, key=lambda r: r['qpb'])[:nsample]
    outsuffix = '_low'
elif mode in ('kids', 'kidslow'):  # random sample of the QPB-open children of the depth-(d-1) nodes DNN did NOT fathom
    miss = set()                     # ('kids': misses of the random sample; 'kidslow': misses of the lowest-QPB sample)
    for l in open(D + 'results_depth%d%s.jsonl' % (depth - 1, '' if mode == 'kids' else '_low')):
        rr = json.loads(l)
        if rr['lb'] <= FATHOM:
            miss.add(tuple(tuple(t) for t in rr['key']))
    cand = [r for r in cand if any(set(r['key']) >= set(mk) for mk in miss)]
    print('children of %d missed parents: %d' % (len(miss), len(cand)), flush=True)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(cand))[:min(nsample, len(cand))]
    rows = [cand[i] for i in idx]
    outsuffix = '_' + mode
else:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(cand))[:min(nsample, len(cand))]
    rows = [cand[i] for i in idx]

outfn = D + 'results_depth%d%s.jsonl' % (depth, outsuffix)
done = set()
if os.path.exists(outfn):
    for l in open(outfn):
        try:
            done.add(tuple(tuple(t) for t in json.loads(l)['key']))
        except Exception:
            pass
print('sample', len(rows), 'already done', len(done), flush=True)
out = open(outfn, 'a')
if mode == 'low':   # reuse nodes already computed in the random sample instead of recomputing them
    rfn = D + 'results_depth%d.jsonl' % depth
    if os.path.exists(rfn):
        for l in open(rfn):
            rr = json.loads(l)
            k = tuple(tuple(t) for t in rr['key'])
            if k in set(tuple(r['key']) for r in rows) and k not in done:
                out.write(l if l.endswith('\n') else l + '\n'); out.flush(); done.add(k)
                print('copied from random sample:', rr['key'], 'lb %.2f' % rr['lb'], flush=True)
vrng = np.random.default_rng(123)
for r in rows:
    if tuple(r['key']) in done:
        continue
    p = r['p']
    verify_reduction(a, b, c, p, vrng)
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    t0 = time.time()
    res = dnn_admm(AFF, BLL, Cp, const, maxit=maxit, check_every=50, fathom_at=FATHOM,
                   time_limit=tlim, stall_window=6)
    rec = dict(key=[list(t) for t in r['key']], m=len(F), inherited=r['inherited'], pred=r['pred'],
               qpb=r['qpb'], branched=r['branched'], lb=res['lb'], lb_last=res['lb_last'], obj=res['obj'],
               it=res['it'], cpu=res['cpu'], wall=time.time() - t0, status=res['status'],
               fathom=bool(res['lb'] > FATHOM),
               hist=[(h[0], round(h[1], 3), round(h[2], 3)) for h in res['hist']])
    out.write(json.dumps(rec) + '\n'); out.flush()
    print('%s m=%d qpb=%s lb=%.2f obj=%.2f it=%d cpu=%.1f %s fathom=%s' % (
        rec['key'], rec['m'], None if r['qpb'] is None else round(r['qpb'], 1), res['lb'], res['obj'],
        res['it'], res['cpu'], res['status'], rec['fathom']), flush=True)
