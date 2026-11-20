#!/usr/bin/env python3
"""failures.py: what rule 9 saw at the labelled nodes (README "Rule 9's failures"). Parses runs/dump/NAME/err.txt
(dumpnode.sh: QAP_LINE_SCORES, QAP_LINE_SCORES_LA, QAP_XBAR_DEBUG=2) and joins the outcomes: rule 9's line's branching
children (the descent, runs/down + runs/d9_*), rule 2's (runs/scan/NAME or runs/arm17/*_rule2), the signal's arms.
Per node: rule 9's candidates (top row, top column by sum x^2) and the lookahead values between them, the chosen
line, the alternative (rule 2's line, or the better arm) and its rank by marginal score and by lookahead; counts of
entries above 0.02 / 0.05 / 0.1 / 1/m on the chosen and the alternative line; features; and the "support" rule's
line (fewest entries above 1/m, ties by marginal score). Writes out/failures.tsv and out/failures.txt."""
import glob, os, re, collections, sys
H = os.path.dirname(os.path.abspath(__file__)); CUT = 608215054.05; ROOT = 605315696.0   # the root bound (tai60b.par)
t = open('/home/user/qap/QAPSOLN/tai60b.sln').read().split(); p0 = [int(x) - 1 for x in t[2:62]]
def trace(paths):
    d = {}
    for p in paths:
        if os.path.exists(p):
            for l in open(p):
                f = l.split()
                if len(f) >= 12 and f[0] == 'node': d[f[2]] = f
    return d
lb = lambda f: float(f[5] if f[5] != '-' else f[4] if f[4] != '-' else f[3])
DOWN = trace(['runs/sc47/nodes.txt'] + glob.glob(H + '/runs/d9_*/nodes.txt') + glob.glob(H + '/runs/down/d*/nodes/*.txt'))
def kids(T, k): return [f for kk, f in T.items() if kk.startswith(k + ',') and kk.count(',') == k.count(',') + 1]
def nb(T, k): return sum(f[7] == 'B' for f in kids(T, k))
ALT = {}   # node fixes -> {line: branching children} from every arm bounded at it
for d in glob.glob(H + '/runs/scan/*') + glob.glob(H + '/runs/arm17/*') + glob.glob(H + '/runs/alt/*'):
    T = trace([d + '/nodes.txt'])
    if not T or not os.path.exists(d + '/cmd.txt') or not any(l.startswith('exit') for l in open(d + '/cmd.txt')): continue
    dep = min(int(f[1]) for f in T.values()); f = [f for f in T.values() if int(f[1]) == dep][0]
    ALT.setdefault(f[2], {})[f[9]] = (nb(T, f[2]), os.path.basename(d))
def parse(path):
    R = dict(U=[], X=[], LArow=None, LAcol=None)
    for l in open(path):
        if l.startswith('ls-node:'):
            v = l.split(); R['depth'] = int(v[2]); R['m'] = int(v[4]); R['bound'] = float(v[6]); R['inc'] = float(v[8])
        elif l.startswith('ls-rows:'): R['rows'] = [int(x) for x in l.split()[1:]]
        elif l.startswith('ls-cols:'): R['cols'] = [int(x) for x in l.split()[1:]]
        elif l.startswith('ls-X '): R['X'].append([float(x) for x in l.split(':', 1)[1].split()])
        elif l.startswith('ls-LArow:'): R['LArow'] = [float(x) for x in l.split()[1:]]
        elif l.startswith('ls-LAcol:'): R['LAcol'] = [float(x) for x in l.split()[1:]]
        elif l.startswith('xbar9:') and 'chosen' not in R:
            m = re.match(r'xbar9: depth (\d+) (row|col) (\d+) \(row (\S+) col (\S+)\)', l)
            R['chosen'] = ('r' if m.group(2) == 'row' else 'c') + m.group(3); R['arb_row'] = float(m.group(4)); R['arb_col'] = float(m.group(5))
    return R
rows = []
for d in sorted(glob.glob(H + '/runs/dump/*')):
    name = os.path.basename(d)
    if not os.path.exists(d + '/done.txt'): continue
    R = parse(d + '/err.txt')
    if not R['X']: print(name, 'no dump', file=sys.stderr); continue
    m = R['m']; X = R['X']; rws = R['rows']; cls = R['cols']
    line_entries = {}
    for i, r in enumerate(rws): line_entries[f'r{r}'] = X[i]
    for j, c in enumerate(cls): line_entries[f'c{c}'] = [X[i][j] for i in range(m)]
    score = {ln: sum(v * v for v in e) for ln, e in line_entries.items()}
    order = sorted(score, key=lambda ln: -score[ln]); rank = {ln: i + 1 for i, ln in enumerate(order)}
    la = {}
    if R['LArow']:
        for i, r in enumerate(rws): la[f'r{r}'] = R['LArow'][i]
        for j, c in enumerate(cls): la[f'c{c}'] = R['LAcol'][j]
    laorder = sorted(la, key=lambda ln: -la[ln]); larank = {ln: i + 1 for i, ln in enumerate(laorder)}
    toprow = max((ln for ln in score if ln[0] == 'r'), key=lambda ln: score[ln]); topcol = max((ln for ln in score if ln[0] == 'c'), key=lambda ln: score[ln])
    T = trace([d + '/nodes.txt']); node = [f for f in T.values() if int(f[1]) == R['depth']][0]; k = node[2]
    ref = DOWN.get(k); r9line = ref[9] if ref else '?'; r9b = nb(DOWN, k) if ref else None
    cnt = lambda ln, th: sum(v > th for v in line_entries[ln])
    sup = sorted(score, key=lambda ln: (cnt(ln, 1.0 / m), -score[ln]))[0]
    alts = ALT.get(k, {}); best_alt = min(alts.items(), key=lambda kv: kv[1][0]) if alts else None
    on = all(p0[int(a)] == int(b) for a, b in (x.split(':') for x in k.split(',')))
    s1, s2 = score[order[0]], score[order[1]]
    rows.append(dict(name=name, depth=R['depth'], m=m, path='yes' if on else 'no', shortfall=CUT - lb(node), sf_root=(CUT - lb(node)) / (CUT - ROOT),
                     chosen=R['chosen'], r9line=r9line, r9b=r9b, toprow=toprow, topcol=topcol, arb_row=R['arb_row'], arb_col=R['arb_col'],
                     top=order[0], s1=s1, s2=s2, margin=s1 - s2, alts=alts, best_alt=best_alt, rank=rank, larank=larank,
                     cnt=cnt, sup=sup, score=score, order=order))
import pickle
out = []
hdr = ['node', 'depth', 'path', 'shortfall', 'sf/root slack', 'rule 9: line (B)', 'cands row/col (LA row, LA col)', 'top score, margin',
       'rule 9 line: entries >0.02/0.05/0.1/>1/m', 'alternative: line (B)', 'alt rank by x^2 / by LA', 'alt entries >0.02/0.05/0.1/>1/m', 'support line (B if known)']
out.append('| ' + ' | '.join(hdr) + ' |'); out.append('|' + '---|' * len(hdr))
for r in sorted(rows, key=lambda r: (r['depth'], r['name'])):
    c = r['cnt']; ln = r['chosen']; m = r['m']
    ba = r['best_alt']
    if ba: aln, (ab, _) = ba
    else: aln, ab = None, None
    ent = lambda L: f"{c(L, .02)}/{c(L, .05)}/{c(L, .1)}/{c(L, 1/m)}"
    alts = ', '.join(f"{l} {v[0]}" for l, v in sorted(r['alts'].items()))
    sb = r['alts'].get(r['sup'], (None,))[0]
    if r['sup'] == r['r9line']: sb = r['r9b']
    out.append(f"| {r['name']} | {r['depth']} | {r['path']} | {r['shortfall']/1e3:,.0f}k | {r['sf_root']:.3f} | {r['chosen']} ({r['r9b']}) | {r['toprow']}/{r['topcol']} ({r['arb_row']:.3g}, {r['arb_col']:.3g}) | "
               f"{r['s1']:.3f}, {r['margin']:.3f} | {ent(ln)} | {alts or '-'} | "
               + (f"{r['rank'][aln]} / {r['larank'].get(aln, '-')}" if aln else '-') + f" | {ent(aln) if aln else '-'} | {r['sup']}" + (f" ({sb})" if sb is not None else '') + ' |')
open(H + '/out/failures.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
pickle.dump([{k: v for k, v in r.items() if k != 'cnt'} for r in rows], open(H + '/out/failures.pkl', 'wb'))
