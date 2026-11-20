"""Q4 analysis of out/q4_obs.tsv (q4_counterfactual.py) -> out/q4.txt
1. noise: the same node under the same line in two sources (different parent states / builds): |delta b|
2. counterfactuals: the same node under different lines; for every pair of lines at a node, does a feature
   known before branching mark the worse line? Features: side (row/column); for tai60b, a column into a cluster-1
   location (the incumbent's facility at the location is one of 14-25); cand (free members of the line's cluster:
   column = the incumbent's cluster at the location, row = locations hosting the row facility's cluster);
   tight (that cluster's interchangeability, smaller = tighter); the node's shortfall is the same for both lines
   (up to the parent-state noise), so it cannot mark the line."""
import itertools
import numpy as np, pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
d = pd.read_csv('out/q4_obs.tsv', sep='\t')
out = open('out/q4.txt', 'w'); w = lambda s='': out.write(s + '\n')
w('# Q4: the same node under two or more lines (q4_counterfactual.py, q4_analyse.py)')
# 1. noise
rows = []
for (inst, node), g in d.groupby(['inst', 'node']):
    for line, h in g.groupby('line'):
        if len(h) >= 2:
            for (_, a), (_, b) in itertools.combinations(h.iterrows(), 2):
                rows.append([inst, int(a.depth), line, a.src, b.src, a.b, b.b, abs(a.b - b.b)])
nz = pd.DataFrame(rows, columns=['inst', 'depth', 'line', 'src1', 'src2', 'b1', 'b2', 'absdiff'])
w(f'\n## 1. the same line twice: {len(nz)} pairs at {nz.groupby(["inst"]).size().to_dict()} (instance: pairs)')
w(f'|delta b|: mean {nz.absdiff.mean():.2f}, share 0: {(nz.absdiff == 0).mean():.0%}, share >= 2: {(nz.absdiff >= 2).mean():.0%}, max {nz.absdiff.max()}')
w(nz[nz.absdiff >= 1].sort_values('absdiff', ascending=False).head(15).to_string(index=False))
# 2. counterfactuals
per = d.groupby(['inst', 'node', 'line']).agg(depth=('depth', 'first'), b=('b', 'mean'), nobs=('b', 'size'),
                                              srcs=('src', lambda x: '+'.join(sorted(set(x)))), side=('side', 'first'),
                                              ocl=('ocl', 'first'), cand=('cand', 'first'), tight=('tight', 'first'),
                                              dclus=('dclus', 'first'), ndom=('ndom', 'first'), path=('path', 'first'),
                                              ham=('ham', 'first'), short=('short', 'mean')).reset_index()
multi = per.groupby(['inst', 'node']).filter(lambda x: len(x) >= 2)
w(f'\n## 2. nodes under two or more lines: {multi.groupby(["inst", "node"]).ngroups} nodes, {len(multi)} (node, line) cells')
w(multi.sort_values(['inst', 'depth', 'node', 'b']).assign(node=lambda x: x.node.str[-34:])[
    ['inst', 'depth', 'node', 'line', 'b', 'srcs', 'side', 'ocl', 'cand', 'tight', 'dclus', 'ndom', 'path', 'ham']].round(3).to_string(index=False))
pairs = []
for (inst, node), g in multi.groupby(['inst', 'node']):
    for (_, a), (_, b) in itertools.combinations(g.iterrows(), 2):
        if a.b < b.b:
            a, b = b, a
        pairs.append(dict(inst=inst, depth=a.depth, node=node, worse=a.line, better=b.line, bw=a.b, bb=b.b, gap=a.b - b.b,
                          w_side=a.side, b_side=b.side, w_cl1=int(inst == 'tai60b' and a.side == 'c' and a.ocl == 1),
                          b_cl1=int(inst == 'tai60b' and b.side == 'c' and b.ocl == 1),
                          w_cand=a.cand, b_cand=b.cand, w_tight=a.tight, b_tight=b.tight,
                          w_dcl1=int(inst == 'tai60b' and a.side == 'c' and a.dclus == 1)))
P = pd.DataFrame(pairs)
w(f'\n## 3. pairs of lines at one node: {len(P)}; gap = b(worse) - b(better)')
for lab, Q in [('all pairs with gap >= 1', P[P.gap >= 1]), ('gap >= 3 (a bad failure: 3 or more extra branching children)', P[P.gap >= 3]),
               ('gap 1-2', P[(P.gap >= 1) & (P.gap < 3)])]:
    if not len(Q):
        continue
    w(f'\n{lab}: {len(Q)} pairs')
    w(f'  worse line is a column: {(Q.w_side == "c").mean():.0%}; better line is a column: {(Q.b_side == "c").mean():.0%}')
    w(f'  worse has more candidates (cand) than better: {(Q.w_cand > Q.b_cand).mean():.0%}, fewer: {(Q.w_cand < Q.b_cand).mean():.0%}, equal: {(Q.w_cand == Q.b_cand).mean():.0%}')
    w(f'  worse line\'s cluster tighter (smaller tight) than better\'s: {(Q.w_tight < Q.b_tight).mean():.0%}, looser: {(Q.w_tight > Q.b_tight).mean():.0%}')
    t6 = Q[Q.inst == 'tai60b']
    if len(t6):
        w(f'  tai60b ({len(t6)} pairs): worse is a cluster-1 column (incumbent at the location): {t6.w_cl1.mean():.0%}; better is: {t6.b_cl1.mean():.0%}; '
          f'worse line\'s branching children are mostly cluster 1 (column): {t6.w_dcl1.mean():.0%}')
w('\npairs with gap >= 3:')
w(P[P.gap >= 3].sort_values(['inst', 'depth']).assign(node=lambda x: x.node.str[-30:])[
    ['inst', 'depth', 'node', 'worse', 'bw', 'better', 'bb', 'w_cand', 'b_cand', 'w_tight', 'b_tight', 'w_cl1', 'w_dcl1']].round(3).to_string(index=False))
# the reverse test: when the feature fires, is the alternative better?
w('\n## 4. the reverse test: at nodes where one line is a tai60b cluster-1 column (by the incumbent) and another is not')
rows = []
for (inst, node), g in multi[multi.inst == 'tai60b'].groupby(['inst', 'node']):
    c1 = g[(g.side == 'c') & (g.ocl == 1)]; ot = g[~((g.side == 'c') & (g.ocl == 1))]
    if len(c1) and len(ot):
        rows.append([int(g.depth.iloc[0]), node[-30:], '/'.join(c1.line), c1.b.mean(), '/'.join(ot.line), ot.b.min(), ot.b.mean()])
R = pd.DataFrame(rows, columns=['depth', 'node', 'cl1 line(s)', 'b cl1', 'other lines', 'b other min', 'b other mean'])
w(R.round(2).to_string(index=False))
if len(R):
    w(f'cl1 column worse than the best other line at {(R["b cl1"] > R["b other min"]).mean():.0%} of {len(R)} nodes, '
      f'equal at {(R["b cl1"] == R["b other min"]).mean():.0%}; sum b: cl1 {R["b cl1"].sum():.0f}, best other {R["b other min"].sum():.0f}')
print(open('out/q4.txt').read())
