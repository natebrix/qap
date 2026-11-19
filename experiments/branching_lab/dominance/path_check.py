import sys, numpy as np, collections
from qapdata import *
R = '/home/user/wt-cognac/infra/phase1/runs'
runs = {'tai60b': f'{R}/tai60b-261117/nodes.txt.gz', 'tai50b': f'{R}/tai50b-261116/nodes.txt.gz',
        'tai35b': '/home/user/wt-cognac/experiments/line_signals/runs/tai35b_base/nodes.txt.gz',
        'nug30': f'{R}/nug30-trace-a/nodes.txt.gz', 'tho40': f'{R}/tho40ab-32/nodes.txt.gz'}
for name, path in runs.items():
    val, p = read_sln(name); q = np.argsort(p)
    on = collections.Counter(); oninv = collections.Counter(); outc = {}
    for d, fx, oc, line, cpu, ln in read_trace(path):
        if all(p[r] == c for r, c in fx): on[d] += 1; outc[d] = oc
        if all(q[r] == c for r, c in fx): oninv[d] += 1
    print(name, 'p-consistent by depth', dict(sorted(on.items())), 'deepest outcome', outc.get(max(on) if on else 0), '| p^-1', dict(sorted(oninv.items())))
