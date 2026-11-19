"""all_lines.py: tai60b's 43:58 (truth/tai60ball.tsv, every free line, ILS 100 kicks): each line's predicted
survivors (ex100 <= s and ex0 <= s) and model cost, for s = 1.76M (the line study's arms) and 1.11M (production);
where the study's four lines and the marginal score's choices rank. Writes out/all_lines.md."""
import json, numpy as np, pandas as pd, lsq
M = json.load(open(lsq.H + '/out/model.json'))
P = lambda x: 1 / (1 + np.exp(-(M['p0'] + M['p1'] * np.log(np.clip(x, 0.01, None)))))
d = pd.read_csv(lsq.H + '/truth/tai60ball.tsv', sep='\t')
out = ['# tai60b 43:58, every free line\n']
for s in (1.76e6, 1.113827e6):
    L = d.groupby('line').agg(cnt=('ex100', lambda v: int((v <= s).sum())), cnt0=('ex0', lambda v: int((v <= s).sum())),
                              mod=('ex100', lambda v: float((P(v.values / s) * np.exp(M['b'] * v.values / s)).sum())))
    L['rank_cnt'] = L.cnt.rank(method='min').astype(int); L['rank_mod'] = L['mod'].rank(method='min').astype(int)
    L = L.sort_values(['mod', 'cnt'])
    out.append(f'## s = {s:,.0f}: {len(L)} lines; the 10 best by the model, then the study\'s four\n')
    out.append('| line | cnt (ex100 <= s) | cnt0 (ex0 <= s) | model | rank by cnt | rank by model |'); out.append('|---|---|---|---|---|---|')
    for ln in list(L.index[:10]) + ['r1', 'r45', 'c26', 'r8']:
        r = L.loc[ln]; out.append(f'| {ln} | {r.cnt} | {r.cnt0} | {r["mod"]:.2f} | {r.rank_cnt} | {r.rank_mod} |')
    out.append(f'\ncnt over all lines: min {L.cnt.min()}, median {L.cnt.median():.0f}, max {L.cnt.max()}\n')
txt = '\n'.join(out); print(txt); open(lsq.H + '/out/all_lines.md', 'w').write(txt + '\n')
