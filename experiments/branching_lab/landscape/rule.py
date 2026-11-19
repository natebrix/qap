#!/usr/bin/env python3
"""rule.py RUN [RUN...]: the survive rule "a child branches iff excess < slack" against the production traces
(out/RUN.tsv from excess.py). A child of a branched node survives when its outcome is B (or E, the leaf).
Slacks (all known before the children are bounded):
  par     cutoff - lb(parent), the parent's shortfall after its cut phase
  pathd   cutoff - lb(on-path node at the child's depth)   (the bound lab's s(d))
  pathd1  cutoff - lb(on-path node at the parent's depth)
  gap     par + ex_k(parent): the parent's estimated gap to its own best completion
Excesses: ex_p0, ex_k, ex_r1, ex_r4, ex_r16, ex_rc (excess.py). Baselines, with the true number of survivors k of
each parent given (an oracle the rule does not get): the k children with the smallest excess, with the smallest
QPB bound, and (for reference only: it is the children's own first solve) with the smallest first DNN bound.
Writes out/rule_RUN.md and prints it."""
import sys, numpy as np, pandas as pd, lsq

def load(run):
    d = pd.read_csv(f'{lsq.H}/out/{run}.tsv', sep='\t', keep_default_na=False, dtype={'fixes': str, 'line': str, 'out': str})
    for c in ['qpb', 'dnn', 'cut', 'cutoff', 'nch', 'cpu'] + [c for c in d.columns if c.startswith('ex_')]:
        d[c] = pd.to_numeric(d[c].replace('-', np.nan))
    d['lb'] = d['cut'].fillna(d['dnn']).fillna(d['qpb'])
    d['parent'] = [','.join(f.split(',')[:-1]) if dp > 1 else ('-' if dp == 1 else None) for f, dp in zip(d.fixes, d.depth)]
    I = lsq.Inst(run.split('-')[0].split('_')[0])
    d['onpath'] = [f == '-' or all(I.p0[int(x.split(':')[0])] == int(x.split(':')[1]) for x in f.split(',')) for f in d.fixes]
    s = d[d.onpath & (d.out == 'B')].set_index('depth')
    sd = (s['cutoff'] - s['lb']).to_dict()
    P = d.set_index('fixes')
    kid = d[d.depth >= 1].copy()
    kid = kid[kid.parent.map(lambda p: p in P.index and P.loc[p, 'out'] == 'B')]
    kid['surv'] = kid.out.isin(['B', 'E'])
    kid['s_par'] = kid.parent.map(lambda p: P.loc[p, 'cutoff'] - P.loc[p, 'lb'])
    kid['s_pathd'] = kid.depth.map(lambda x: sd.get(x, np.nan))
    kid['s_pathd1'] = kid.depth.map(lambda x: sd.get(x - 1, np.nan))
    kid['s_gap'] = kid.s_par + kid.parent.map(lambda p: max(P.loc[p, 'ex_k'], 0))
    kid['run'] = run
    return kid, I

def prf(pred, truth):
    tp = int((pred & truth).sum()); fp = int((pred & ~truth).sum()); fn = int((~pred & truth).sum())
    pr = tp / (tp + fp) if tp + fp else float('nan'); rc = tp / (tp + fn) if tp + fn else float('nan')
    return tp, fp, fn, pr, rc

def oracle_k(kid, col):
    """per parent: the k = (true survivors) children with the smallest col predicted to survive"""
    pred = pd.Series(False, index=kid.index)
    for _, g in kid.groupby('parent'):
        k = int(g.surv.sum())
        if k: pred[g[col].fillna(np.inf).nsmallest(k).index] = True
    return pred

def fmt(x): return '-' if x != x else f'{x:.2f}'

if __name__ == '__main__':
    out = []
    for run in sys.argv[1:]:
        kid, I = load(run)
        out.append(f'## {run}: {len(kid)} children of {kid.parent.nunique()} branched nodes, {int(kid.surv.sum())} survive '
                   f'(B/E), cutoff {kid.cutoff.iloc[0]:,.0f}\n')
        out.append('All children; precision / recall of "survives iff excess <= slack" (tp, fp, fn):\n')
        exs = [c for c in ['ex_p0', 'ex_k', 'ex_r1', 'ex_r4', 'ex_r16', 'ex_rc'] if kid[c].notna().all()]
        out.append('| excess \\ slack | ' + ' | '.join(['par', 'pathd', 'pathd1', 'gap']) + ' |')
        out.append('|---|---|---|---|---|')
        for ex in exs:
            cells = []
            for s in ['s_par', 's_pathd', 's_pathd1', 's_gap']:
                ok = kid[s].notna()
                tp, fp, fn, pr, rc = prf((kid[ex] <= kid[s])[ok], kid.surv[ok])
                cells.append(f'{fmt(pr)} / {fmt(rc)} ({tp}, {fp}, {fn})')
            out.append(f'| {ex} | ' + ' | '.join(cells) + ' |')
        out.append('\nOracle k per parent (precision = recall when k is right): smallest excess vs smallest QPB vs smallest first DNN:\n')
        out.append('| ranking by | precision = recall |')
        out.append('|---|---|')
        for c in exs + ['qpb', 'dnn']:
            tp, fp, fn, pr, rc = prf(oracle_k(kid, c), kid.surv)
            out.append(f'| {c} | {fmt(pr)} ({tp} of {tp + fn}) |')
        # by depth, best excess (ex_k), slack par and pathd
        out.append('\nBy depth, ex_k (knows p0): rule with slack par / pathd (precision, recall), oracle-k excess and QPB:\n')
        out.append('| depth of child | children | survive | par: P / R | pathd: P / R | gap: P / R | oracle-k ex_k | oracle-k qpb |')
        out.append('|---|---|---|---|---|---|---|---|')
        bins = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 7), (8, 10), (11, 15), (16, 25), (26, 60)]
        ok_ex = oracle_k(kid, 'ex_k'); ok_q = oracle_k(kid, 'qpb')
        for a, b in bins:
            m = (kid.depth >= a) & (kid.depth <= b)
            if not m.any(): continue
            g = kid[m]; cells = []
            for s in ['s_par', 's_pathd', 's_gap']:
                ok = g[s].notna(); tp, fp, fn, pr, rc = prf((g.ex_k <= g[s])[ok], g.surv[ok]); cells.append(f'{fmt(pr)} / {fmt(rc)}')
            cells.append(fmt(prf(ok_ex[m], g.surv)[3])); cells.append(fmt(prf(ok_q[m], g.surv)[3]))
            out.append(f'| {a}-{b} | {len(g)} | {int(g.surv.sum())} | ' + ' | '.join(cells) + ' |')
        out.append('\nThreshold multiplier, ex_k <= alpha x slack: F1 (precision / recall):\n')
        al = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
        out.append('| slack | ' + ' | '.join(f'alpha {a}' for a in al) + ' |')
        out.append('|---|' + '---|' * len(al))
        for s in ['s_par', 's_pathd', 's_gap']:
            ok = kid[s].notna(); cells = []
            for a in al:
                tp, fp, fn, pr, rc = prf((kid.ex_k <= a * kid[s])[ok], kid.surv[ok])
                cells.append(f'{2 * tp / (2 * tp + fp + fn):.2f} ({fmt(pr)} / {fmt(rc)})')
            out.append(f'| {s[2:]} | ' + ' | '.join(cells) + ' |')
        kid.to_pickle(f'{lsq.H}/out/kid_{run}.pkl')
        out.append('')
    txt = '\n'.join(out); print(txt)
    open(f'{lsq.H}/out/rule_{"_".join(sys.argv[1:])}.md', 'w').write(txt)
