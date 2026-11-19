#!/usr/bin/env python3
"""score.py [ARM ...]: each arm's depth-8 children (runs/ARM/nodes.txt, a node's last line kept) against the
landscape prediction (out/predict.tsv) and production (tai60b-261117's trace, for c54). Prints, per arm, the
counts (D first solve, C cut phase, B branched), CPU, and child by child: outcome, bound, shortfall
(cutoff - final bound; negative = fathomed with that margin), x at production's slack and at the arm's own
(the arm's 10:39 cutoff - cut bound), under both search starts. Writes out/score_ARM.tsv."""
import sys, os, gzip, collections
H = os.path.dirname(os.path.abspath(__file__))
NODE = '44:45,43:58,1:33,6:24,7:19,0:32,10:39'
def trace(path):
    op = gzip.open if path.endswith('.gz') else open
    d = {}
    for ln in op(path, 'rt'):
        f = ln.split()
        if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
fv = lambda s: None if s == '-' else float(s)
def lb(f): return fv(f[5]) if f[5] != '-' else (fv(f[4]) if f[4] != '-' else fv(f[3]))
pred = collections.defaultdict(dict)
for l in open(H + '/out/predict.tsv').read().split('\n')[1:]:
    if l: s, ln, ch, e0, ek, x = l.split('\t'); pred[(s, ln)][ch] = float(ek)
PROD = '/home/user/wt-cognac/infra/phase1/runs/tai60b-261117/nodes.txt.gz'
prod = None
for arm in sys.argv[1:] or ['c54', 'c47', 'r31']:
    p = f'{H}/runs/{arm}/nodes.txt'
    if not os.path.exists(p): continue
    T = trace(p); top = T.get(NODE)
    line = arm[1:] if arm.startswith('s') else arm
    print(f'== arm {arm}')
    if top:
        s_arm = float(top[6]) - lb(top)
        print(f'10:39: dnn {top[4]} cut {top[5]} cutoff {top[6]} shortfall {s_arm:,.0f} line {top[9]} cpu {float(top[11]):.0f}')
    else:
        print('10:39 not bounded yet'); continue
    s_prod = 608215054.05 - 607809132.8147
    kids = [f for k, f in T.items() if int(f[1]) == 8 and k.startswith(NODE + ',')]
    if line == 'c54':
        prod = prod or trace(PROD)
    cnt = collections.Counter(f[7] for f in kids); cpu = sum(float(f[11]) for f in kids)
    print(f'children bounded {len(kids)}: D {cnt["D"]} C {cnt["C"]} B {cnt["B"]} other {len(kids) - sum(cnt[k] for k in "DCB")}; children CPU {cpu:,.0f} s; 10:39 + children {cpu + float(top[11]):,.0f} s')
    rows = []
    for f in kids:
        ch = f[2].rsplit(',', 1)[1]
        ex = pred[('p0', line)].get(ch); exi = pred[('ils', line)].get(ch)
        po = '-'
        if prod is not None and line == 'c54':
            g = prod.get(f[2]); po = '-' if g is None else f'{g[7]} {float(g[6]) - lb(g):,.0f}'
        rows.append((ex if ex is not None else 1e18, ch, f[7], lb(f), float(f[6]) - lb(f), ex, exi, float(f[11]), po))
    rows.sort()
    w = open(f'{H}/out/score_{arm}.tsv', 'w')
    w.write('child\toutcome\tlb\tshortfall\tex_p0\tx_prod\tx_arm\tex_ils\tx_ils_prod\tcpu\tproduction\n')
    print('child  out  shortfall    x(prod s)  x(arm s)  x_ils(prod s)  cpu  production(out, shortfall)')
    for _, ch, o, b, sf, ex, exi, c, po in rows:
        xp = ex / s_prod if ex is not None else float('nan'); xa = ex / s_arm if ex is not None else float('nan')
        xi = exi / s_prod if exi is not None else float('nan')
        print(f'{ch:6s} {o}  {sf:>12,.0f}  {xp:6.2f}  {xa:6.2f}  {xi:6.2f}  {c:7.0f}  {po}')
        w.write(f'{ch}\t{o}\t{b:.1f}\t{sf:.0f}\t{ex}\t{xp:.3f}\t{xa:.3f}\t{exi}\t{xi:.3f}\t{c:.1f}\t{po}\n')
    w.close()
