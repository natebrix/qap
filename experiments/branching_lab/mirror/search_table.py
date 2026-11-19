"""The whole tai35b search with the automatic rule (runs/NAME) against the production trace (line_signals/runs/tai35b_base,
same settings, build 8df95737 via the coordinator): nodes and CPU in all and by subtree of the root's children, the
nodes the rule touched (pair, mu, outcome, and that node's outcome in production), the path check (every node on the
optimum's path, its outcome), the answer. usage: search_table.py NAME [BASE_TRACE]"""
import sys, gzip, re, os, collections
sys.path.insert(0, '/home/user/qap/experiments/branching_lab/dominance')
from qapdata import read_sln
H = os.path.dirname(os.path.abspath(__file__)); name = sys.argv[1]
base = sys.argv[2] if len(sys.argv) > 2 else '/home/user/wt-cognac/experiments/line_signals/runs/tai35b_base/nodes.txt.gz'
opt, p = read_sln('tai35b')
def load(path):
    op = gzip.open if path.endswith('.gz') else open
    out = {}
    for ln in op(path, 'rt'):
        t = ln.split()
        if t and t[0] == 'node': out[t[2]] = t
    return out
def norm(k):  # fixes in fixing order; key on the set for matching across runs
    return k if k == '-' else ','.join(sorted(k.split(','), key=lambda s: int(s.split(':')[0])))
R = load(f'{H}/runs/{name}/nodes.txt'); Bt = load(base)
rootn = load(f'{H}/runs/root/nodes.txt'); R.update(rootn)
Bn = {norm(k): v for k, v in Bt.items()}
def tot(T): return len(T), sum(float(t[11]) for t in T.values() if t[11] != '-')
print('| | nodes | CPU s |'); print('|---|---|---|')
for lab, T in [('production (coordinator trace)', Bt), (f'with the rule ({name}, root from runs/root)', R)]:
    n, c = tot(T); print(f'| {lab} | {n:,} | {c:,.0f} |')
print()
print('by subtree of the root child (nodes / CPU s): production -> with the rule')
for top in sorted({k.split(',')[0] for k in Bt if k != '-'} | {k.split(',')[0] for k in R if k != '-'}, key=lambda s: int(s.split(':')[1])):
    nb = [t for k, t in Bt.items() if k.split(',')[0] == top]; nr = [t for k, t in R.items() if k.split(',')[0] == top]
    cb = sum(float(t[11]) for t in nb); cr = sum(float(t[11]) for t in nr)
    if len(nb) > 1 or len(nr) > 1: print(f'  {top}: {len(nb)} / {cb:,.0f} -> {len(nr)} / {cr:,.0f}')
print()
err = open(f'{H}/runs/{name}/err.txt').read().splitlines()
auto = [l for l in err if l.startswith('mirror-auto')]
app = [l for l in err if l.startswith('mirror: ')]
print(f'mirror-auto lines {len(auto)}; applied at {len(app)} nodes:')
for l in app:
    fx = l.split(' | ')[0].split()[4:]; key = norm(','.join(fx))
    m = re.search(r'pair (\S+) mu (\S+)', l); f = re.search(r'fate (\S)', l).group(1)
    b = Bn.get(key); bo = b[7] if b else 'not in production'
    print(f'  {key}: pair {m.group(1)} mu {m.group(2)}  outcome {f} | production: {bo}')
chk = [l for l in err if l.startswith('mirror-check')]
worst = max([float(l.split()[-1]) for l in chk], default=0)
print(f'mirror-check lines {len(chk)}, worst difference {worst}')
print()
# path check
path = []
for k, t in R.items():
    if k == '-': continue
    if all(p[int(s.split(':')[0])] == int(s.split(':')[1]) for s in k.split(',')): path.append((int(t[1]), k, t[7], t[4], t[5]))
path.sort()
fath = [x for x in path if x[2] in 'DCQP']
print(f'path check: {len(path)} nodes on the optimum\'s path in the trace, depths {path[0][0]}..{path[-1][0]}; fathomed: {len(fath)}; outcomes ' + ''.join(x[2] for x in path))
log = open(f'{H}/runs/{name}/log.txt').read()
print('New objective values:', re.findall(r'New objective value: (\S+)', log))
st = open(f'{H}/runs/{name}/out.state').read() if os.path.exists(f'{H}/runs/{name}/out.state') else ''
print('open nodes in the final state:', len([l for l in st.splitlines() if l.startswith('node ')]))
