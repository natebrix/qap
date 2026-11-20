"""childlp.py PARENT LINE K [EPS]: the AJ LP of K children of LINE at a tai60b truth node (truth.py), spread over the
range of their first DNN bound, by PDLP to EPS (default 1e-4). Prints, per child, its QPB, LP (valid dual bound),
DNN and outcome; appends to out/childlp.tsv. T1's 'LP of the child' predictor at m = 42."""
import sys, os, numpy as np, ajlp, truth
from scipy.stats import spearmanr
par, line, K = sys.argv[1], sys.argv[2], int(sys.argv[3]); eps = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4
nodes = truth.load(); ch = truth.parents(nodes)[par][line]
ch = sorted(ch, key=lambda t: t[1]['dnn'])
pick = [ch[int(round(i))] for i in np.linspace(0, len(ch) - 1, K)]
A, B, opt, perm = ajlp.inst('tai60b')
res = []
with open('out/childlp.tsv', 'a') as o:
    for fx, v in pick:
        Ar, Br, C, s, I, J = ajlp.reduce(A, B, par + ',' + fx)
        r = ajlp.pdlp_solve(ajlp.LP(Ar, Br, C), eps=eps, threads=int(os.environ.get('THREADS', 3)))
        row = (par.rsplit(',', 1)[-1], line, fx, v['qpb'], r['lagr'] + s, r['obj'] + s, v['dnn'], v['out'], r['wall'], r['cpu'])
        res.append(row); o.write('\t'.join(str(x) for x in row) + '\n'); o.flush()
        print(*row, flush=True)
q = [r[3] for r in res]; l = [r[4] for r in res]; d = [r[6] for r in res]
print('Spearman with child DNN: child LP', round(spearmanr(l, d)[0], 3), 'child QPB', round(spearmanr(q, d)[0], 3))
