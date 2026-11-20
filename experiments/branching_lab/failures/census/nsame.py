"""A cluster count as predictor and rule (the coordinator's follow-up, Nov 20).  -> out/nsame*.txt, out/nsame_nodes.tsv.gz

n_same(line) at a node, from the instance's facility clusters (clusters.py) and the incumbent (= best known) p:
  column c_l: k = cluster of the incumbent's facility at l;  n_same = free facilities in cluster k
  row r_f:    k = cluster of f;                              n_same = free locations whose incumbent facility is in k
On tai60b's mirror chain the owner is the mirror's: the mirror exchanges facilities 43 and 44, which are the same
cluster, so every location has the same owner under both and the mirror needs no separate case. Off both chains
the incumbent's owner is used (the only owner known without the solver's marginals); results are split on/off chain.

Part 1: n_same of the taken line against b.  Part 2a: the taken line's place among every free line's n_same.
Part 2b: counterfactual pairs (out/q4_obs.tsv plus the lines measured in linetest/lsla READMEs).  Part 2c: what
"min n_same" picks at named nodes.  Tie-break for min n_same: the looser cluster first (larger tightness value,
clusters.py: members less interchangeable), then the line with more free children outside the owner cluster
(identical within a node, so in effect:) then the smallest index, rows and columns alternating by index -- no side
preference. The full tied set is printed.
"""
import gzip, json, itertools
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats
from common import *

import os
CL = json.load(open('out/clusters.json'))
# NSAME_MERGE=1: tai60b with clusters 1 (14-25) and 4 (48-51) merged, the five-cluster cut (MST edge 162 joins them;
# clusters.py). Outputs get the suffix _m14.
SUF = ''
if os.environ.get('NSAME_MERGE') == '1':
    c = CL['tai60b']['fac_clusters']
    CL['tai60b']['fac_clusters'] = [c[0], sorted(c[1] + c[4]), c[2], c[3], c[5]]
    SUF = '_m14'
RUNS_N = ['tai60b-261117', 'tai50b-261116', 'tai50b-r9', 'tai60b-c47', 'tai50b-ls']
BANDS = [-1, 2, 6, 10, 15, 22, 99]; BL = ['0-2', '3-6', '7-10', '11-15', '16-22', '23+']


class Inst:
    def __init__(self, inst):
        self.n, A, B = read_dat(inst)
        self.obj, self.p = read_sln(inst)
        self.q = [0] * self.n
        for f, l in enumerate(self.p):
            self.q[l] = f
        self.cls = CL[inst]['fac_clusters']
        self.fcl = {f: k for k, c in enumerate(self.cls) for f in c}
        m = np.array(A, dtype=float); mean = m.mean()
        self.tight = {k: (float(np.mean([np.abs(np.delete(m[i], [i, j]) - np.delete(m[j], [i, j])).mean() / mean
                                         for i, j in itertools.combinations(c, 2)])) if len(c) > 1 else 0.0)
                      for k, c in enumerate(self.cls)}
        pm = None
        if inst == 'tai60b':
            pm = list(self.p); pm[43], pm[44] = self.p[44], self.p[43]
        self.pm = pm

    def all_lines(self, fixes):
        ff = set(f for f, l in fixes); fl = set(l for f, l in fixes)
        freeF = [f for f in range(self.n) if f not in ff]; freeL = [l for l in range(self.n) if l not in fl]
        cntF = defaultdict(int); cntL = defaultdict(int)
        for f in freeF: cntF[self.fcl[f]] += 1
        for l in freeL: cntL[self.fcl[self.q[l]]] += 1
        out = {}
        for f in freeF:
            k = self.fcl[f]; out[f'r{f}'] = (cntL[k], k)
        for l in freeL:
            k = self.fcl[self.q[l]]; out[f'c{l}'] = (cntF[k], k)
        return out

    def pick(self, lines):
        """min n_same; ties: looser cluster first, then index (rows and columns interleaved)."""
        return min(lines, key=lambda ln: (lines[ln][0], -self.tight[lines[ln][1]], int(ln[1:]), ln[0]))

    def chain(self, fixes):
        if all(self.p[f] == l for f, l in fixes): return 'path'
        if self.pm and all(self.pm[f] == l for f, l in fixes): return 'mirror'
        return 'off'


def part1_2a():
    rows = []
    for run in RUNS_N:
        I = Inst(RUNS[run]['inst'])
        nodes, _ = load(run)
        md = min(y.depth for y in nodes.values())
        root_short = 2899357.0 if run == 'tai60b-c47' else [x for x in nodes.values() if x.depth == md][0].short
        for n in nodes.values():
            if n.out != 'B' or n.kids is None or len(n.kidsl) != n.kids:
                continue
            L = I.all_lines(n.fixes)
            if n.line not in L:
                continue
            b = sum(1 for c in n.kidsl if c.out == 'B')
            ns, k = L[n.line]
            vals = np.array([v[0] for v in L.values()])
            mn = vals.min()
            pk = I.pick(L)
            rows.append(dict(run=run, depth=n.depth, key=n.key, rule=rule_of(run, n), line=n.line, side=n.line[0], b=b,
                             nsame=ns, ocl=k, chain=I.chain(n.fixes), rel_short=(n.short or 0) / root_short,
                             nmin=mn, nmax=vals.max(), nlines=len(vals), rank_le=int((vals < ns).sum()),
                             n_at_min=int((vals == mn).sum()), rows_at_min=int(sum(1 for ln, v in L.items() if v[0] == mn and ln[0] == 'r')),
                             pick=pk, pick_ns=L[pk][0], pick_ocl=L[pk][1], n_le2=int((vals <= 2).sum()),
                             lsig=0))
    d = pd.DataFrame(rows)
    d['band'] = pd.cut(d.depth, BANDS, labels=BL)
    d.to_csv(f'out/nsame_nodes{SUF}.tsv.gz', sep='\t', index=False, compression='gzip')
    return d


def report1(d, w):
    w('## Part 1: n_same of the taken line against b')
    for run in RUNS_N:
        g = d[d.run == run]
        w(f'\n### {run}: {len(g)} branched nodes (on chain {int((g.chain != "off").sum())}, off {int((g.chain == "off").sum())})')
        rows = []
        for (bd, ch), h in g.groupby(['band', g.chain.where(g.chain == 'off', 'chain')], observed=True):
            if len(h) < 5:
                rows.append([bd, ch, len(h), np.nan, h.b.mean(), h.nsame.mean(), np.nan]); continue
            rho = stats.spearmanr(h.b, h.nsame).statistic if h.nsame.nunique() > 1 and h.b.nunique() > 1 else np.nan
            rs = stats.spearmanr(h.b, h.rel_short).statistic if h.b.nunique() > 1 else np.nan
            rows.append([bd, ch, len(h), rho, h.b.mean(), h.nsame.mean(), rs])
        w(pd.DataFrame(rows, columns=['band', 'chain/off', 'n', 'rho(b,n_same)', 'mean b', 'mean n_same', 'rho(b,rel_short)']).round(2).to_string(index=False))
        h = g.copy(); h['nb'] = pd.cut(h.nsame, [-1, 0, 2, 4, 6, 8, 10, 12, 15, 20, 60])
        t = h.groupby('nb', observed=True).b.agg(['size', 'mean', lambda x: x.quantile(0.9), 'max'])
        t.columns = ['n', 'mean b', 'p90 b', 'max b']
        w('b by n_same of the taken line:')
        w(t.round(2).to_string())
        # b / n_same by band and by shortfall quintile, nodes with n_same >= 3
        h = g[g.nsame >= 3].copy()
        if len(h) >= 10:
            h['ratio'] = h.b / h.nsame
            h['sq'] = pd.qcut(h.rel_short.rank(method='first'), 4, labels=['q1 (least short)', 'q2', 'q3', 'q4 (most short)'])
            t = h.groupby(['band', 'sq'], observed=True).ratio.agg(['size', 'mean']).unstack('sq')
            w('b / n_same (nodes with n_same >= 3), by depth band x shortfall quartile of the run (mean, n):')
            t2 = pd.DataFrame(index=t.index)
            for c in t.columns.levels[1]:
                if ('mean', c) in t.columns:
                    t2[c] = [f'{m:.2f} ({int(s)})' if s == s else '' for m, s in zip(t[('mean', c)], t[('size', c)])]
            w(t2.to_string())


def report2a(d, w):
    w('\n## Part 2a: where the taken line ranks among all free lines by n_same')
    rows = []
    for (run, rule), g in d.groupby(['run', 'rule']):
        ex = g.nsame - g.nmin
        rows.append([run, rule, len(g), (ex == 0).mean(), (ex >= 3).mean(), (ex >= 6).mean(), ex.mean(),
                     g.nmin.median(), (g.nmin <= 2).mean(), g.n_le2.median(),
                     (g.rows_at_min == g.n_at_min).mean(), (g.rows_at_min == 0).mean()])
    w(pd.DataFrame(rows, columns=['run', 'rule', 'nodes', 'taken = min', 'taken >= min+3', 'taken >= min+6',
                                  'mean excess', 'median min', 'min <= 2', 'median #lines <= 2',
                                  'min set all rows', 'min set all columns']).round(2).to_string(index=False))
    w('(excess = n_same(taken) - the node\'s minimum over all free lines. "min set all rows/columns": the lines '
      'attaining the minimum are all rows / all columns; otherwise both sides tie.)')
    w('\nb of the taken line by its excess over the node\'s minimum (all five runs):')
    h = d.copy(); h['ex'] = pd.cut(h.nsame - h.nmin, [-1, 0, 2, 5, 9, 60])
    t = h.groupby(['rule', 'ex'], observed=True).b.agg(['size', 'mean', lambda x: x.quantile(0.9), lambda x: (x >= 4).mean()])
    t.columns = ['n', 'mean b', 'p90 b', 'share b>=4']
    w(t.round(2).to_string())
    w('\nWhat min n_same would pick, by side and cluster (the pick\'s owner cluster), per run:')
    for run, g in d.groupby('run'):
        w(f'  {run}: picks rows {(g.pick.str[0] == "r").mean():.0%}; owner cluster of the pick '
          f'{g.pick_ocl.value_counts().sort_index().to_dict()}; pick n_same median {g.pick_ns.median():.0f}, '
          f'taken n_same median {g.nsame.median():.0f}; a line with n_same <= 2 exists at {(g.nmin <= 2).mean():.0%} of nodes')


MANUAL = [  # (inst, fixes, line, b, source)  lines bounded elsewhere, from the READMEs, strong parent states
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39', 'c59', 3, 'lsla sc59'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39', 'r11', 4, 'lsla sr11'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39', 'c47', 5, 'linetest sc47'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39', 'c54', 9, 'linetest sc54'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54', 'c59', 3, 'lsla at8'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54', 'r3', 4, 'lsla at8'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54', 'c2', 7, 'lsla at8'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'c59', 3, 'lsla at9'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'r8', 4, 'lsla at9'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'c53', 4, 'lsla at9'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'r3', 4, 'lsla at9'),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'c30', 12, 'lsla at9'),
]


def part2b(w):
    w('\n## Part 2b: counterfactual pairs: does the line with the smaller n_same have the smaller b?')
    obs = pd.read_csv('out/q4_obs.tsv', sep='\t')
    per = obs.groupby(['inst', 'node', 'line']).agg(b=('b', 'mean'), depth=('depth', 'first'), src=('src', lambda x: '+'.join(sorted(set(x))))).reset_index()
    cells = {}
    for r in per.itertuples():
        cells[(r.inst, r.node, r.line)] = (r.b, r.depth, r.src)
    for inst, fx, line, b, src in MANUAL:
        node = ','.join(sorted(fx.split(','), key=lambda s: tuple(map(int, s.split(':')))))
        node = ','.join(f'{a}:{c}' for a, c in sorted(tuple(map(int, s.split(':'))) for s in fx.split(',')))
        if (inst, node, line) not in cells:
            cells[(inst, node, line)] = (b, len(fx.split(',')), src)
    insts = {i: Inst(i) for i in set(k[0] for k in cells)}
    bynode = defaultdict(dict)
    for (inst, node, line), v in cells.items():
        bynode[(inst, node)][line] = v
    res = []; detail = []
    for (inst, node), L in bynode.items():
        if len(L) < 2:
            continue
        I = insts[inst]
        fixes = [] if node == '-' else [tuple(map(int, s.split(':'))) for s in node.split(',')]
        A = I.all_lines(fixes)
        for (l1, v1), (l2, v2) in itertools.combinations(L.items(), 2):
            n1, n2 = A[l1][0], A[l2][0]; b1, b2 = v1[0], v2[0]
            if n1 == n2:
                out = 'n_same tied'
            elif b1 == b2:
                out = 'b tied'
            else:
                out = 'won' if ((n1 < n2) == (b1 < b2)) else 'lost'
            res.append(dict(inst=inst, depth=v1[1], node=node, l1=l1, n1=n1, b1=b1, l2=l2, n2=n2, b2=b2, out=out,
                            gap=abs(b1 - b2)))
    R = pd.DataFrame(res)
    w(f'{len(R)} pairs at {R.groupby(["inst", "node"]).ngroups} nodes (q4_obs.tsv plus the README lines: tai60b 10:39, 14:54, 24:2)')
    w(R.out.value_counts().to_string())
    big = R[R.gap >= 3]
    w(f'\npairs with a b gap of 3 or more: {len(big)}')
    w(big.out.value_counts().to_string())
    w('\nby instance and depth (won / lost / b tied / n_same tied):')
    t = R.groupby(['inst', 'depth']).out.value_counts().unstack(fill_value=0)
    w(t.to_string())
    w('\nall pairs with a b gap of 3 or more:')
    w(big.sort_values(['inst', 'depth']).assign(node=lambda x: x.node.str[-28:])[['inst', 'depth', 'node', 'l1', 'n1', 'b1', 'l2', 'n2', 'b2', 'out']].to_string(index=False))
    R.to_csv(f'out/nsame_pairs{SUF}.tsv', sep='\t', index=False)


NAMED = [
    ('tai50b', '-', 'tai50b root', {'r31': 10, 'c26': 28}),
    ('tai60b', '44:45,43:58', 'tai60b depth 2, 43:58 (tai60b_lines)', {'r1': 6, 'r45': 10, 'r8': 12, 'c26': 20}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39', 'tai60b depth 7, 10:39', {'c59': 3, 'r11': 4, 'c47': 5, 'c54': 9}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54', 'tai60b depth 8, 14:54', {'c59': 3, 'r3': 4, 'c2': 7}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2', 'tai60b depth 9, 24:2', {'c59': 3, 'r8': 4, 'c53': 4, 'r3': 4, 'c30': 12}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:49', 'depth 17, 8:49', {'c30': 12, 'r51': 5, 'r2': 0}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:57', 'depth 17, 8:57', {'c30': 6, 'r51': 1, 'r2': 0}),
    ('tai60b', '44:45,43:58,1:33,6:24,7:19,0:32,10:39,47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6,22:31', 'depth 17, 22:31 (path)', {'c26': 7, 'r13': 2, 'c54': 3}),
]


def part2c(w):
    w('\n## Part 2c: what min n_same picks at named nodes (n_same of every measured line in brackets)')
    for inst, fx, lab, meas in NAMED:
        I = Inst(inst)
        fixes = [] if fx == '-' else [tuple(map(int, s.split(':'))) for s in fx.split(',')]
        A = I.all_lines(fixes)
        pk = I.pick(A); mn = A[pk][0]
        ties = sorted([ln for ln, v in A.items() if v[0] == mn], key=lambda s: (s[0], int(s[1:])))
        hist = defaultdict(int)
        for v in A.values(): hist[v[0]] += 1
        w(f'- {lab}: picks **{pk}** (n_same {mn}, owner cluster {A[pk][1]}); tied at the minimum: {" ".join(ties)}; '
          f'measured: ' + ', '.join(f'{ln} b={b} n_same={A[ln][0]}' for ln, b in meas.items()) +
          f'; n_same over all free lines: {dict(sorted(hist.items()))}')


if __name__ == '__main__':
    d = part1_2a()
    with open(f'out/nsame{SUF}.txt', 'w') as f:
        w = lambda s='': f.write(s + '\n')
        w('# A cluster count as predictor and rule (nsame.py)')
        report1(d, w); report2a(d, w); part2b(w); part2c(w)
    print(open(f'out/nsame{SUF}.txt').read())
