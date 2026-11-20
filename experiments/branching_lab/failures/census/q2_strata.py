"""Q2, controlled: is the cluster-1 column effect in tai60b more than shortfall and depth?  -> out/q2_strata.txt
Within strata of (depth band, quintile of rel_short within the band), mean b and P(b >= 4) for three kinds of line:
columns whose branching children (or, at b = 0, the incumbent's facility at the location) are cluster 1 (facilities
14-25), other columns, rows. Then a Poisson regression of b on log(rel_short), depth band and line kind."""
import numpy as np, pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
out = open('out/q2_strata.txt', 'w')
w = lambda s='': out.write(s + '\n')
for run, rule in [('tai60b-261117', 2), ('tai60b-261117', 9)]:
    d = pd.read_csv(f'out/bn_{run}.tsv.gz', sep='\t'); d = d[(d.nkids == d.kidsf) & (d.rule == rule)].copy()
    d['band'] = pd.cut(d.depth, [-1, 2, 6, 10, 15, 22, 99], labels=['0-2', '3-6', '7-10', '11-15', '16-22', '23+'])
    d['clus'] = d.dclus.fillna(d.ocl)
    d['kind'] = np.where(d.side == 'r', 'row', np.where(d.clus == 1, 'col-cl1', 'col-other'))
    d['q'] = d.groupby('band', observed=True).rel_short.transform(lambda x: pd.qcut(x.rank(method='first'), 5, labels=False) if len(x) >= 5 else 0)
    t = d.groupby(['band', 'q', 'kind'], observed=True).b.agg(['size', 'mean', lambda x: (x >= 4).mean()]).unstack('kind')
    w(f'## {run}, rule {rule} levels: mean b (n) by depth band x shortfall quintile (0 = smallest) x line kind')
    t2 = pd.DataFrame(index=t.index)
    for k in ['col-cl1', 'col-other', 'row']:
        if ('size', k) in t.columns:
            t2[k] = [f'{m:.2f} ({int(n)})' if n == n else '' for m, n in zip(t[('mean', k)], t[('size', k)])]
    w(t2.to_string())
    try:
        import statsmodels.formula.api as smf
        import statsmodels.api as sm
        d['lrs'] = np.log(d.rel_short.clip(lower=1e-5))
        m = smf.glm('b ~ lrs + C(band) + C(kind, Treatment("row"))', data=d, family=sm.families.Poisson()).fit()
        w('\nPoisson regression of b (coefficients are log rate ratios):')
        w(m.summary2().tables[1][['Coef.', 'Std.Err.', 'P>|z|']].round(3).to_string())
    except Exception as e:
        w(f'(no regression: {e})')
    w()
print(open('out/q2_strata.txt').read())
