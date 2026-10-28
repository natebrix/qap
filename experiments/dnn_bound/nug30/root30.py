import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn')
import numpy as np
from qaputil import *; from dnn import dnn_admm
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug30.dat')
res = dnn_admm(a, b, c, 0.0, maxit=6000, check_every=100, tol=1e-5, time_limit=600)
print(json.dumps(dict(lb=res['lb'], obj=res['obj'], it=res['it'], cpu=res['cpu'], status=res['status'], hist=[(h[0], round(h[1], 2), round(h[2], 2)) for h in res['hist']])))
