"""truth.py: the bounded-children sets at tai60b nodes (linetest, lsla runs). parents() -> {parent fixes: {line:
[(child fix 'r:v', qpb, dnn, cut, cutoff, outcome)]}}, only lines whose children are all present."""
import glob, collections
L = '/home/user/qap/experiments/branching_lab'
FILES = (glob.glob(L + '/linetest/runs/s*/nodes.txt') + glob.glob(L + '/linetest/runs/arm17/*/nodes.txt') +
         glob.glob(L + '/linetest/runs/down/d18/nodes/*.txt') + glob.glob(L + '/lsla/runs/*/nodes.txt'))


def load(files=FILES):
    import gzip, os
    if os.environ.get('INST', 'tai60b') != 'tai60b':
        files = sorted(glob.glob(f'/home/user/wt-cognac/experiments/line_signals/runs/{os.environ["INST"]}_*/nodes.txt.gz'))
    nodes = {}
    for p in files:
        for ln in (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)):
            f = ln.split()
            if len(f) < 12 or f[0] != 'node': continue
            nodes[f[2]] = dict(depth=int(f[1]), qpb=float('nan') if f[3] == '-' else float(f[3]), dnn=None if f[4] == '-' else float(f[4]),
                               cut=None if f[5] == '-' else float(f[5]), cutoff=float(f[6]), out=f[7], line=f[9], cpu=float(f[11]))
    return nodes


def parents(nodes, n=None):
    import os
    n = n or int(''.join(ch for ch in os.environ.get('INST', 'tai60b') if ch.isdigit()))
    kids = collections.defaultdict(dict)
    for fx, v in nodes.items():
        if fx == '-': continue
        par, last = fx.rsplit(',', 1) if ',' in fx else ('-', fx); kids[par][last] = v
    res = {}
    for par, ch in kids.items():
        fixed = [] if par == '-' else [tuple(map(int, f.split(':'))) for f in par.split(',')]
        fr = {r for r, _ in fixed}; fc = {c for _, c in fixed}
        m = n - len(fixed); lines = collections.defaultdict(list)
        for last, v in ch.items():
            r, c = map(int, last.split(':'))
            lines[f'r{r}'].append((last, v)); lines[f'c{c}'].append((last, v))
        full = {ln: sorted(v) for ln, v in lines.items() if len(v) == m}
        if full: res[par] = full
    return res


if __name__ == '__main__':
    nd = load(); P = parents(nd)
    for par, ls in sorted(P.items(), key=lambda t: len(t[0])):
        pn = nd.get(par, {})
        print(0 if par == '-' else par.count(',') + 1, par[-40:], 'parent', pn.get('qpb'), pn.get('dnn'), pn.get('cut'), {ln: ''.join(sorted(o[1]['out'] for o in v)).count('B') for ln, v in ls.items()})
