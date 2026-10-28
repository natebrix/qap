"""Run DNN ADMM to convergence (no early stop) on selected nodes to get actual DNN values."""
import sys, json
import numpy as np
from qaputil import *; from dnn import dnn_admm
D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn/'
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug20.dat')
fn = sys.argv[1]; recs = [json.loads(l) for l in open(D + fn)]
out = open(D + 'converged_' + fn, 'a')
for idx in map(int, sys.argv[2:]):
    r = recs[idx]; p = [-1] * n
    for i, l in r['key']: p[i] = l
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    res = dnn_admm(AFF, BLL, Cp, const, maxit=4000, check_every=100, tol=1e-5, time_limit=200)
    rec = dict(file=fn, idx=idx, key=r['key'], qpb=r['qpb'], lb_early=r['lb'], lb=res['lb'], obj=res['obj'], it=res['it'], cpu=res['cpu'], status=res['status'])
    out.write(json.dumps(rec) + '\n'); out.flush()
    print(rec, flush=True)
