"""Q1, by kind: what share of tai60b-261117 lies in the subtrees OPENED (off-chain branching children and all below)
by given sets of branched nodes, against random sets of the same size at the same depths.  -> out/q1_sets.txt"""
import random, json
from collections import defaultdict
import pandas as pd
from common import *
nodes, tops = load('tai60b-261117')
obj, p = read_sln('tai60b'); pm = list(p); pm[43], pm[44] = p[44], p[43]
chain = set(id(x) for x in nodes.values() if all(p[f] == l for f, l in x.fixes) or all(pm[f] == l for f, l in x.fixes))
order = sorted(nodes.values(), key=lambda x: x.depth)
allN = len(nodes); allC = sum(x.cpu for x in nodes.values())
def opened(sel):
    sel = set(id(x) for x in sel); flag = {}; un = uc = 0
    for x in order:
        f = x.parent is not None and (flag.get(id(x.parent), False) or (id(x.parent) in sel and x.out == 'B' and id(x) not in chain))
        flag[id(x)] = f
        if f: un += 1; uc += x.cpu
    return un / allN, uc / allC
d = pd.read_csv('out/bn_tai60b-261117.tsv.gz', sep='\t')
byk = {r.key: r for r in d.itertuples()}
br = [x for x in nodes.values() if x.out == 'B']
brd = defaultdict(list)
for x in br: brd[x.depth].append(x)
rng = random.Random(2)
def base(s, reps=20):
    need = defaultdict(int)
    for x in s: need[x.depth] += 1
    a = c = 0
    for _ in range(reps):
        u, v = opened([y for dd, k in need.items() for y in rng.sample(brd[dd], k)]); a += u; c += v
    return a / reps, c / reps
def kind(r): return 'row' if r.side == 'r' else ('col-cl1' if r.ocl == 1 else 'col-other')
sets = {
 'b >= 4 (all)': [x for x in br if byk[x.key].b >= 4],
 'b >= 4, rule 9 levels (depths 0-6)': [x for x in br if byk[x.key].b >= 4 and byk[x.key].rule == 9],
 'b >= 4, rule 2 levels (depths 7+)': [x for x in br if byk[x.key].b >= 4 and byk[x.key].rule == 2],
 'b >= 4, rule 2, column into a cluster-1 location': [x for x in br if byk[x.key].b >= 4 and byk[x.key].rule == 2 and kind(byk[x.key]) == 'col-cl1'],
 'every rule-2 column into a cluster-1 location (any b)': [x for x in br if byk[x.key].rule == 2 and kind(byk[x.key]) == 'col-cl1'],
 'b >= 4, rule 2, any other line': [x for x in br if byk[x.key].b >= 4 and byk[x.key].rule == 2 and kind(byk[x.key]) != 'col-cl1'],
 'chain nodes (path and mirror)': [x for x in br if id(x) in chain],
}
with open('out/q1_sets.txt', 'w') as f:
    f.write('# tai60b-261117: subtrees opened by sets of branched nodes (q1_sets.py)\n')
    f.write(f'{allN} nodes, {allC/3600:.0f} CPU h. "opened" = the set\'s branching children off the chain (path, mirror) and everything below them, union.\n\n')
    f.write('| set | nodes in set | their b | opened: nodes | opened: CPU | random, same depths: nodes, CPU |\n|---|---|---|---|---|---|\n')
    for k, s in sets.items():
        u, c = opened(s); ru, rc = base(s)
        f.write(f'| {k} | {len(s)} | {sum(byk[x.key].b for x in s)} | {u:.1%} | {c:.1%} | {ru:.1%}, {rc:.1%} |\n')
print(open('out/q1_sets.txt').read())
