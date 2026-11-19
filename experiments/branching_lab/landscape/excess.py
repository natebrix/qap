#!/usr/bin/env python3
"""excess.py RUN [RAND_MAXDEPTH] [norc]: every node of a production trace gets the excess (cost - best known) of the
best completion a 2-exchange descent finds for it, under several searches:
  ex_p0    one descent from the best known permutation p0, adjusted to the node's fixes (one exchange per fix)
  ex_par   one descent from the parent's best completion (min over the searches below that know p0), adjusted
  ex_mir   (tai60b) one descent from the mirror (p0 with the facilities at locations 43 and 44 exchanged)
  ex_k     min of ex_p0, ex_par, ex_mir: the "knows the optimum" search
  ex_r1, ex_r4, ex_r16   random completions of the free rows only, best of 1, 4, 16 descents (does not know p0);
           computed for nodes of depth <= RAND_MAXDEPTH (default 12) only
  ex_rc    the no-knowledge chain: one descent from the parent's ex_rc completion + 4 random (root: 16 random)
Writes out/RUN.tsv (one line per node). Instances: tai60b, tai50b, tai35b, tai40b by the run's name."""
import sys, os, time, numpy as np, lsq

run = sys.argv[1]
RMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 12
RC = not (len(sys.argv) > 3 and sys.argv[3] == 'norc')   # 'norc': skip ex_rc (tai60b: too slow for 331k nodes)
inst = run.split('-')[0].split('_')[0]
I = lsq.Inst(inst)
path = f'{lsq.WT}/infra/phase1/runs/{run}/nodes.txt.gz'
if not os.path.exists(path): path = f'{lsq.WT}/experiments/line_signals/runs/{run}/nodes.txt.gz'
T = sorted(lsq.read_trace(path), key=lambda t: t['depth'])
mir = None
if inst == 'tai60b':
    mir = lsq.mirror(I)
best_k, best_rc = {}, {}
os.makedirs(lsq.H + '/out', exist_ok=True)
w = open(f'{lsq.H}/out/{run}.tsv', 'w')
cols = ['depth', 'fixes', 'out', 'qpb', 'dnn', 'cut', 'cutoff', 'line', 'nch', 'cpu', 'ex_p0', 'ex_par', 'ex_mir', 'ex_k',
        'ex_r1', 'ex_r4', 'ex_r16', 'ex_rc']
w.write('\t'.join(cols) + '\n')
t0 = time.time(); seed = 1
for idx, t in enumerate(T):
    fix = dict(t['fx']); par = ','.join(t['fixes'].split(',')[:-1]) if t['depth'] > 1 else '-'
    e = {}
    c, p = I.descend(I.adjust(I.p0, fix), fix); e['ex_p0'] = c - I.best; bk = (c, p)
    if par in best_k:
        c, p = I.descend(I.adjust(best_k[par], fix), fix); e['ex_par'] = c - I.best
        if c < bk[0]: bk = (c, p)
    if mir is not None:
        c, p = I.descend(I.adjust(mir, fix), fix); e['ex_mir'] = c - I.best
        if c < bk[0]: bk = (c, p)
    e['ex_k'] = bk[0] - I.best; best_k[t['fixes']] = bk[1]
    if t['depth'] <= RMAX:
        cs = []
        for k in range(16):
            seed += 1
            rng = np.random.default_rng(seed)
            c, p = I.descend(I.adjust(rng.permutation(I.n).astype(np.int32), fix), fix); cs.append(c)
        e['ex_r1'] = cs[0] - I.best; e['ex_r4'] = min(cs[:4]) - I.best; e['ex_r16'] = min(cs) - I.best
    if RC:
        start = best_rc.get(par)
        seed += 1
        if start is None: c, p = I.multistart(I.adjust(np.random.default_rng(seed).permutation(I.n), fix), fix, 15, seed)
        else: c, p = I.multistart(I.adjust(start, fix), fix, 4, seed)
        e['ex_rc'] = c - I.best; best_rc[t['fixes']] = p
    row = [t['depth'], t['fixes'], t['out'], t['qpb'], t['dnn'], t['cut'], t['cutoff'], t['line'], t['nch'], t['cpu']] + \
          [e.get(k) for k in cols[10:]]
    w.write('\t'.join('-' if v is None else (f'{v:.1f}' if isinstance(v, float) else str(v)) for v in row) + '\n')
    if idx % 20000 == 0: print(run, idx, len(T), f'{time.time() - t0:.0f}s', flush=True)
w.close()
print(run, 'done', len(T), f'{time.time() - t0:.0f}s', 'moves', I.moves.value)
