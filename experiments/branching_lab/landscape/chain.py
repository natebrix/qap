#!/usr/bin/env python3
"""chain.py RUN DMAX [KICKS]: at every chain node of depth 0..DMAX-1 of a production trace (a branched node holding
the best known solution p0 or, tai60b, its mirror), every child of every free line gets its excess: the best of the
descents from p0 and the mirror adjusted to the child's fixes, then an ILS of KICKS (default 30) kicks; ex0 and
exK are both kept. One descent per distinct child (a child r:v is in row r's line and column v's).
Writes chain/RUN.tsv: node, depth, node lb, cutoff, the line production took, line, child, ex0, exK, and the
child's production outcome where the trace has it (the taken line's children)."""
import sys, os, time, numpy as np, lsq

run, DMAX = sys.argv[1], int(sys.argv[2]); KICKS = int(sys.argv[3]) if len(sys.argv) > 3 else 30
inst = run.split('-')[0].split('_')[0]
I = lsq.Inst(inst)
path = f'{lsq.WT}/infra/phase1/runs/{run}/nodes.txt.gz'
if not os.path.exists(path): path = f'{lsq.WT}/experiments/line_signals/runs/{run}/nodes.txt.gz'
T = lsq.read_trace(path); byfix = {t['fixes']: t for t in T}
starts = [I.p0]
if inst == 'tai60b':
    mir = lsq.mirror(I)
    starts.append(mir)
rng = np.random.default_rng(5)
def excess(fix):
    c, p = min((I.descend(I.adjust(s, fix), fix) for s in starts), key=lambda x: x[0]); c0 = c
    free = np.array([i for i in range(I.n) if i not in fix])
    for _ in range(KICKS):
        q = p.copy()
        for _ in range(6):
            i, j = rng.choice(free, 2, replace=False); q[i], q[j] = q[j], q[i]
        cq, q = I.descend(q, fix)
        if cq < c: c, p = cq, q
    return c0 - I.best, c - I.best
chain = [t for t in T if t['out'] == 'B' and t['depth'] < DMAX and
         any(all(s[r] == v for r, v in t['fx']) for s in starts)]
if os.environ.get('ONLY') == 'mirror':      # the mirror chain alone (written to chain/RUN_mirror.tsv)
    chain = [t for t in chain if not all(I.p0[r] == v for r, v in t['fx'])]
chain.sort(key=lambda t: (t['depth'], t['fixes']))
os.makedirs(lsq.H + '/chain', exist_ok=True)
w = open(f'{lsq.H}/chain/{run}{"_mirror" if os.environ.get("ONLY") == "mirror" else ""}.tsv', 'w')
w.write('node\tdepth\tlb\tcutoff\ttaken\tline\tchild\tex0\texK\tout\n'); t0 = time.time()
for t in chain:
    nf = dict(t['fx']); lbv = lsq.lb(t); memo = {}
    fr = [r for r in range(I.n) if r not in nf]; fc = [c for c in range(I.n) if c not in nf.values()]
    for side, ks in (('r', fr), ('c', fc)):
        for k in ks:
            ch = [(k, v) for v in fc] if side == 'r' else [(r, k) for r in fr]
            for (r, v) in ch:
                if (r, v) not in memo:
                    fix = dict(nf); fix[r] = v; memo[(r, v)] = excess(fix)
                e0, ek = memo[(r, v)]
                cf = (t['fixes'] + ',' if t['fixes'] != '-' else '') + f'{r}:{v}'
                out = byfix[cf]['out'] if f'{side}{k}' == t['line'] and cf in byfix else '-'
                w.write(f"{t['fixes']}\t{t['depth']}\t{lbv:.1f}\t{t['cutoff']:.2f}\t{t['line']}\t{side}{k}\t{r}:{v}\t{e0:.0f}\t{ek:.0f}\t{out}\n")
    w.flush(); print(run, t['depth'], t['fixes'][-30:], len(memo), f'{time.time() - t0:.0f}s', flush=True)
w.close()
