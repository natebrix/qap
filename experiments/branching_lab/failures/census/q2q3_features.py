"""Q2 (what high-b nodes have in common) and Q3 (sibling repetition), from out/bn_<run>.tsv.gz.
-> out/q2.txt, out/q3.txt.  Only branched nodes whose children are all in the trace.
'high' = b >= 4 (tai60b: the top 5.5% of branched nodes; the other runs' shares are printed)."""
import sys
import numpy as np, pandas as pd
from scipy import stats

pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200); pd.set_option('display.max_columns', 30)
RUNS = ['tai60b-261117', 'tai60b-c47', 'tai50b-261116', 'tai50b-r9', 'tai50b-ls', 'tai35b_base', 'tho40ab-32', 'nug30-trace-a']
BANDS = [-1, 2, 6, 10, 15, 22, 99]
HI = 4


def load(run):
    d = pd.read_csv(f'out/bn_{run}.tsv.gz', sep='\t')
    d = d[d.nkids == d.kidsf].copy()
    d['band'] = pd.cut(d.depth, BANDS, labels=['0-2', '3-6', '7-10', '11-15', '16-22', '23+'])
    d['hi'] = d.b >= HI
    return d


def bsh(row):
    if not isinstance(row.bshorts, str):
        return []
    return [float(x) for x in row.bshorts.split(';')]


def q2(out):
    w = lambda s='': out.write(s + '\n')
    w('# Q2: what the high-b nodes (b >= 4) have in common (q2q3_features.py)')
    for run in RUNS:
        d = load(run)
        w(f'\n## {run}: {len(d)} branched nodes, {d.hi.sum()} with b >= {HI} ({d.hi.mean():.1%}), '
          f'they hold {d[d.hi].b.sum() / max(d.b.sum(), 1):.0%} of all branching children')
        # side x rule x band
        t = d.groupby(['rule', 'band', 'side'], observed=True).agg(n=('b', 'size'), mean_b=('b', 'mean'), max_b=('b', 'max'),
                                                                   hi=('hi', 'sum')).reset_index()
        t['hi_rate'] = (t.hi / t.n).round(3); t['mean_b'] = t.mean_b.round(2)
        w('\nby rule, depth band and side (r = row: a facility\'s free locations; c = column: a location\'s free facilities)')
        w(t.to_string(index=False))
        # path / mirror / ham
        if d.path.notna().any():
            d['pos'] = np.where(d.path == 1, 'path', np.where(d.mirror == 1, 'mirror',
                                np.where(d.ham <= 2, 'ham1-2', np.where(d.ham <= 5, 'ham3-5', 'ham6+'))))
            t = d.groupby(['band', 'pos'], observed=True).agg(n=('b', 'size'), mean_b=('b', 'mean'), hi=('hi', 'sum')).reset_index()
            t['mean_b'] = t.mean_b.round(2)
            w('\nby depth band and position (path / mirror / number of fixes off the best known solution, or its mirror)')
            w(t.to_string(index=False))
        # shortfall and QPB gap: within each band, the high nodes' median percentile among the band's nodes
        w('\nshortfall (cutoff - bound) relative to the root\'s, and QPB\'s distance below the cutoff: per band, the '
          'median of the high nodes and of the rest; Spearman of b with each')
        rows = []
        for bd, g in d.groupby('band', observed=True):
            if len(g) < 5:
                continue
            r1 = stats.spearmanr(g.b, g.rel_short).statistic if g.rel_short.nunique() > 1 else np.nan
            r2 = stats.spearmanr(g.b, g.qpb_gap).statistic if g.qpb_gap.nunique() > 1 else np.nan
            rows.append([bd, len(g), g.hi.sum(), g[g.hi].rel_short.median(), g[~g.hi].rel_short.median(), r1,
                         g[g.hi].qpb_gap.median(), g[~g.hi].qpb_gap.median(), r2])
        w(pd.DataFrame(rows, columns=['band', 'n', 'hi', 'rel_short hi', 'rel_short rest', 'rho(b,rel_short)',
                                      'qpb_gap hi', 'qpb_gap rest', 'rho(b,qpb_gap)']).round(3).to_string(index=False))
        # clusters (tai-b only)
        if d.ocl.notna().any() and run.startswith('tai'):
            c = d[(d.side == 'c') & (d.b >= 2)]
            r = d[(d.side == 'r') & (d.b >= 2)]
            w(f'\nclusters: columns with b >= 2: {len(c)}; all branching children in ONE facility cluster: '
              f'{(c.nclus_b == 1).mean():.0%}; children in the commonest cluster / b: {(c.maxclus_b / c.b).mean():.2f}; '
              f'that cluster is the incumbent\'s at the location: {(c.dclus == c.ocl).mean():.0%}; b / free members of it: '
              f'{(c.b / c.fdom).mean():.2f} (b >= {HI}: {(c[c.hi].b / c[c.hi].fdom).mean():.2f})')
            if len(r):
                w(f'          rows with b >= 2: {len(r)}; branching children that are locations of the row facility\'s cluster '
                  f'(in the incumbent): {r.bhost.sum() / r.b.sum():.0%}; b / such free locations: {(r.b / r.hostfree).mean():.2f}')
            cc = d[d.side == 'c'].copy()
            cc['clus'] = cc.dclus.fillna(cc.ocl)
            t = cc.groupby(['rule', 'clus']).agg(n=('b', 'size'), mean_b=('b', 'mean'), max_b=('b', 'max'), hi=('hi', 'sum')).reset_index()
            t['mean_b'] = t.mean_b.round(2)
            w('columns by the cluster their branching children come from (the incumbent\'s cluster at the location when b = 0):')
            w(t.to_string(index=False))
            rr = d[d.side == 'r']
            t = rr.groupby(['rule', 'ocl']).agg(n=('b', 'size'), mean_b=('b', 'mean'), max_b=('b', 'max'), hi=('hi', 'sum')).reset_index()
            t['mean_b'] = t.mean_b.round(2)
            w('rows by the row facility\'s cluster:')
            w(t.to_string(index=False))
        # recurring lines
        t = d.groupby('line').agg(n=('b', 'size'), mean_b=('b', 'mean'), max_b=('b', 'max'), sum_b=('b', 'sum'),
                                  hi=('hi', 'sum')).sort_values('hi', ascending=False)
        t['mean_b'] = t.mean_b.round(2)
        w(f'\nlines with the most high-b nodes (of {d.line.nunique()} lines used):')
        w(t.head(12).to_string())
        # children: how the rest were fathomed; near misses against deep survivors
        rows = []
        for lab, g in [('b>=4', d[d.hi]), ('b 1-3', d[(d.b >= 1) & (d.b < HI)])]:
            fr = []
            for _, x in g.iterrows():
                s = bsh(x)
                if x.short and x.short > 0 and s:
                    fr += [v / x.short for v in s]
            fr = np.array(fr)
            fath = (g.nD + g.nC + g.nQ).sum()
            rows.append([lab, len(g), g.nC.sum() / max(fath, 1), g.nQ.sum() / max(fath, 1), g.near_fath.mean(),
                         (fr < 0.25).mean() if len(fr) else np.nan, (fr > 0.5).mean() if len(fr) else np.nan,
                         np.median(fr) if len(fr) else np.nan])
        w('\nchildren: share of the fathomed ones fathomed by the cut phase (C) and by QPB (Q); near-miss fathoms per node '
          '(D or C within 1% of the parent\'s shortfall above the cutoff); branching children\'s shortfall / the parent\'s: '
          'share < 0.25 (near misses), share > 0.5 (deep survivors), median')
        w(pd.DataFrame(rows, columns=['nodes', 'n', 'C share', 'Q share', 'near fathoms/node', 'surv <0.25', 'surv >0.5',
                                      'median']).round(3).to_string(index=False))


def q3(out):
    w = lambda s='': out.write(s + '\n')
    w('# Q3: sibling repetition (q2q3_features.py)')
    w('A sibling group = the branched children of one node (two or more, all children traced). "same line" = every '
      'sibling chose one line. ICC = intraclass correlation of b within groups (one-way ANOVA estimator).')
    for run in RUNS:
        d = load(run)
        g = d.groupby('parent_key')
        groups = [x for _, x in g if len(x) >= 2]
        if not groups:
            continue
        same = [x.line.nunique() == 1 for x in groups]
        # ICC(1)
        k = np.mean([len(x) for x in groups]); allb = pd.concat(groups).b
        msb = sum(len(x) * (x.b.mean() - allb.mean()) ** 2 for x in groups) / (len(groups) - 1)
        msw = sum(((x.b - x.b.mean()) ** 2).sum() for x in groups) / (len(allb) - len(groups))
        icc = (msb - msw) / (msb + (k - 1) * msw) if (msb + (k - 1) * msw) > 0 else np.nan
        # high-b mass in groups on a common line where two or more siblings are high
        hi_nodes = d[d.hi]
        in_rep = 0; in_same = 0
        for x in groups:
            if x.line.nunique() == 1:
                in_same += x[x.hi].b.sum()
            # siblings sharing this node's line with b >= 3
            for line, y in x.groupby('line'):
                if (y.b >= 3).sum() >= 2:
                    in_rep += y[y.hi].b.sum()
        tot_hi = hi_nodes.b.sum()
        no_sib = hi_nodes[~hi_nodes.parent_key.isin([x.parent_key.iloc[0] for x in groups])].b.sum()
        w(f'\n## {run}: {len(groups)} sibling groups (mean size {k:.1f}); all siblings on one line: {np.mean(same):.0%}; '
          f'ICC of b {icc:.2f}')
        w(f'   b of the high (b >= {HI}) nodes: {tot_hi}; in groups where every sibling took the same line: {in_same} '
          f'({in_same / max(tot_hi, 1):.0%}); in a same-line set of >= 2 siblings with b >= 3: {in_rep} ({in_rep / max(tot_hi, 1):.0%}); '
          f'high nodes with no branched sibling: b {no_sib} ({no_sib / max(tot_hi, 1):.0%})')
        # by band
        rows = []
        for bd in d.band.cat.categories:
            gg = [x for x in groups if x.band.iloc[0] == bd]
            if not gg:
                continue
            rows.append([bd, len(gg), np.mean([x.line.nunique() == 1 for x in gg]),
                         np.mean([len(x) for x in gg]),
                         stats.spearmanr(*zip(*[(a, b2) for x in gg if len(x) >= 2 for a, b2 in
                                                [(x.b.iloc[0], x.b.iloc[1])]])).statistic if len(gg) >= 5 else np.nan])
        w(pd.DataFrame(rows, columns=['band', 'groups', 'one line', 'mean size', 'rho(b of two siblings)']).round(2).to_string(index=False))
        # the largest groups of high-b siblings on one line
        rep = []
        for x in groups:
            for line, y in x.groupby('line'):
                if (y.b >= 3).sum() >= 2:
                    rep.append((y.b.sum(), len(y), line, int(y.depth.iloc[0]), x.parent_key.iloc[0][-40:], list(y.b)))
        rep.sort(reverse=True)
        w('   largest same-line sibling sets with two or more b >= 3 (sum b, siblings, line, depth, ...parent, b each):')
        for r in rep[:8]:
            w(f'   {r}')


if __name__ == '__main__':
    with open('out/q2.txt', 'w') as f:
        q2(f)
    with open('out/q3.txt', 'w') as f:
        q3(f)
    print(open('out/q2.txt').read()[:20000])
    print(open('out/q3.txt').read())
