"""Write the reduced 3-matrix QAPLIB instance (A_FF, B_LL, C') of chosen nug30 nodes for bnbdrv.
usage: python3 write_node.py results_depthD.jsonl idx1 idx2 ...   -> est/node_dD_idx.dat ; prints const
"""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn')
from qaputil import *
D = '/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnn30/'
n, a, b, c = read_qaplib('/home/user/qap/qaplib/nug30.dat')
recs = [json.loads(l) for l in open(D + sys.argv[1])]
tag = sys.argv[1].split('_')[1].split('.')[0]
for idx in map(int, sys.argv[2:]):
    r = recs[idx]
    p = [-1] * n
    for i, l in r['key']:
        p[i] = l
    F, L, AFF, BLL, Cp, const = reduce_subproblem(a, b, c, p)
    assert np.allclose(Cp, np.round(Cp)) and float(const).is_integer()
    fn = D + 'est/node_%s_%d.dat' % (tag, idx)
    with open(fn, 'w') as f:
        f.write('%d\n\n' % len(F))
        for M in (AFF, BLL, Cp):
            for row in M:
                f.write(' '.join('%d' % x for x in row) + '\n')
            f.write('\n')
    print('%s idx=%d key=%s m=%d const=%d incumbent_for_reduced=%d qpb=%.1f dnn_lb=%.1f' % (
        fn, idx, r['key'], len(F), const, 6124 - const, r['qpb'], r['lb']))
