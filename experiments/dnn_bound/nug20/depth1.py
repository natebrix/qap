import json, numpy as np
from qaputil import *; from dnn import dnn_admm
D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn/'
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug20.dat')
_, opened = read_state_nodes(D + 'n20_d0.state'); _, kids = read_state_nodes(D + 'n20_d1.state')
own = {}
for k in kids: own[parent_key(k['p'], k['hist'], k['depth'])] = k['bound']
out = open(D + 'results_depth1.jsonl', 'a')
for x in opened:
    key = tuple(sorted((i, l) for i, l in enumerate(x['p']) if l >= 0))
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, x['p'])
    res = dnn_admm(AFF, BLL, Cp, const, maxit=4000, check_every=100, fathom_at=2568.05, time_limit=120)
    rec = dict(key=[list(t) for t in key], m=len(F), qpb=own.get(key), lb=res['lb'], obj=res['obj'], it=res['it'], cpu=res['cpu'], status=res['status'])
    out.write(json.dumps(rec) + '\n'); out.flush(); print(rec, flush=True)
