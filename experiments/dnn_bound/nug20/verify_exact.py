"""Exact check: write each chosen nug20 node's reduced instance (A_FF, B_LL, C') as a 3-matrix
QAPLIB file, solve it with bnbdrv, and compare optimum+const with the certified DNN lb.
usage: python3 verify_exact.py results_depth2.jsonl idx1 idx2 ...
"""
import sys, json, subprocess, re, time
import numpy as np
from qaputil import *
from dnn import dnn_admm

D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn/'
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug20.dat')
recs = [json.loads(l) for l in open(D + sys.argv[1])]
for idx in map(int, sys.argv[2:]):
    r = recs[idx]
    p = [-1] * n
    for i, l in r['key']:
        p[i] = l
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    assert np.allclose(Cp, np.round(Cp)) and float(const).is_integer()
    fn = D + 'node_%s_%d.dat' % (sys.argv[1].split('_')[1].split('.')[0], idx)
    with open(fn, 'w') as f:
        f.write('%d\n\n' % len(F))
        for M in (AFF, BLL, Cp):
            for row in M:
                f.write(' '.join('%d' % x for x in row) + '\n')
            f.write('\n')
    t0 = time.time()
    out = subprocess.run(['./bnbdrv', '-p', 'PARAM/default.par', fn], cwd='/home/user/qap/latest',
                         capture_output=True, text=True).stdout
    vals = [float(x) for x in re.findall(r'New objective value:\s*([-\d.]+)', out)]
    nodes = re.search(r'Nodes:\s*(\d+)', out)
    opt = min(vals) + const
    print('node %d key=%s m=%d const=%d subproblem optimum=%.0f  DNN lb=%.3f  valid=%s  QPB=%s  solver_nodes=%s  wall=%.1fs' % (
        idx, r['key'], len(F), const, opt, r['lb'], r['lb'] <= opt + 1e-6, r['qpb'], nodes.group(1) if nodes else '?', time.time() - t0), flush=True)
