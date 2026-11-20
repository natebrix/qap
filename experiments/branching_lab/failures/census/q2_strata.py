"""Q2, controlled: is the cluster-1 column effect in tai60b more than shortfall and depth?  -> out/q2_strata.txt
Within strata of (depth band, quintile of rel_short within the band), mean b and P(b >= 4) for three kinds of line:
columns whose location holds a cluster-1 facility (14-25) in the incumbent (known before branching; it is also the
cluster of the branching children at every column with b >= 2 in this run), other columns, rows. Then a Poisson regression of b on log(rel_short), depth band and line kind."""
import numpy as np, pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
out = open('out/q2_strata.txt', 'w')
w = lambda s='': out.write(s + '\n')
for run, rule in [('tai60b-261117', 2), ('tai60b-261117', 9)]:
    d = pd.read_csv(f'out/bn_{run}.tsv.gz', sep='\t'); d = d[(d.nkids == d.kidsf) & (d.rule == rule)].copy()
    d['band'] = pd.cut(d.depth, [-1, 2, 6, 10, 15, 22, 99], labels=['0-2', '3-6', '7-10', '11-15', '16-22', '23+'])
    d['clus'] = d.ocl   # the incumbent's cluster at the location (column) / the facility's cluster (row): known before branching
    d['kind'] = np.where(d.side == 'r', 'row', np.where(d.clus == 1, 'col-cl1', 'col-other'))
    d['q'] = d.groupby('band', observed=True).rel_short.transform(lambda x: pd.qcut(x.rank(method='first'), 5, labels=False) if len(x) >= 5 else 0)
    t = d.groupby(['band', 'q', 'kind'], observed=True).b.agg(['size', 'mean', lambda x: (x >= 4).mean()]).unstack('kind')
    w(f'## {run}, rule {rule} levels: mean b (n) by depth band x shortfall quintile (0 = smallest) x line kind')
    t2 = pd.DataFrame(index=t.index)
    for k in ['col-cl1', 'col-other', 'row']:
        if ('size', k) in t.columns:
            t2[k] = [f'{m:.2f} ({int(n)})' if n == n else '' for m, n in zip(t[('mean', k)], t[('size', k)])]
    w(t2.to_string())
    # Poisson regression by IRLS (no statsmodels here)
    d = d[d.band.isin(d.band.value_counts()[lambda x: x >= 30].index)]
    X = pd.get_dummies(d[['band', 'kind']].astype(str), drop_first=False)
    X = X.drop(columns=[c for c in X.columns if c in ('kind_row',)] + [X.filter(like='band_').columns[0]])
    X.insert(0, 'log_rel_short', np.log(d.rel_short.clip(lower=1e-5)))
    X.insert(0, 'const', 1.0)
    A = X.values.astype(float); y = d.b.values.astype(float); beta = np.zeros(A.shape[1]); beta[0] = np.log(y.mean())
    for _ in range(50):
        mu = np.exp(A @ beta); W = mu; z = A @ beta + (y - mu) / mu
        beta_new = np.linalg.solve(A.T @ (A * W[:, None]), A.T @ (W * z))
        if np.max(np.abs(beta_new - beta)) < 1e-9: beta = beta_new; break
        beta = beta_new
    mu = np.exp(A @ beta); cov = np.linalg.inv(A.T @ (A * mu[:, None]))
    disp = np.sum((y - mu) ** 2 / mu) / (len(y) - len(beta))
    se = np.sqrt(np.diag(cov) * max(disp, 1))
    w(f'\nPoisson regression of b, quasi-Poisson s.e. (dispersion {disp:.2f}); rate ratio = exp(coef); baseline: rows, first band')
    w(pd.DataFrame({'coef': beta, 'se': se, 'rate ratio': np.exp(beta)}, index=X.columns).round(3).to_string())
    w()
print(open('out/q2_strata.txt').read())
