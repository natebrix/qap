"""t2.py PARENT_FIXES [EPS] [SOLVER]: at a tai60b node with bounded children (truth.py), solve the AJ LP of the parent
(PDLP to EPS, default 1e-4; SOLVER 'ipm' for HiGHS IPM), and from its row dual u the valid bound of every child x_ij = 1
(ajlp.child_bounds: Lagrangian with the child's incompatible variables dropped). Then:
  T1: Spearman, over the bounded children, of {child-bound, child's own QPB} against the child's first DNN bound;
  T2: per line, the rule-2 analogue: count of children whose bound reaches the cutoff, sum of child bounds, sum of d;
      for the lines with truth, their branching-children counts.
Writes out/t2_<tag>.txt and out/t2_<tag>_lines.tsv (every free line)."""
import sys, os, time, numpy as np, ajlp, truth
from scipy.stats import spearmanr
par = sys.argv[1]; eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4; solver = sys.argv[3] if len(sys.argv) > 3 else 'pdlp'
tag = f'd{par.count(",") + 1}_{par.rsplit(",", 1)[-1].replace(":", "_")}_{solver}{eps:g}'
nodes = truth.load(); P = truth.parents(nodes)[par]; pn = nodes[par]
A, B, opt, perm = ajlp.inst('tai60b')
Ar, Br, C, s, I, J = ajlp.reduce(A, B, par)
lp = ajlp.LP(Ar, Br, C)
cache = f'/tmp/claude-0/lpdual_{tag}.npy'
if os.path.exists(cache):
    u = np.load(cache); r = dict(cpu=0, wall=0, iters=0, obj=float('nan'))
elif solver == 'ipm':
    r = lp.solve('ipm'); u = r['u']; np.save(cache, u)
else:
    r = ajlp.pdlp_solve(lp, eps=eps); u = r['u']; np.save(cache, u)
cb, base = lp.child_bounds(u); cb = cb + s; base += s
cutoff = pn['cutoff']; out = open(f'out/t2_{tag}.txt', 'w')
def say(*a):
    t = ' '.join(map(str, a)); print(t, flush=True); out.write(t + '\n')
say(f'parent {par} m={lp.m} QPB {pn["qpb"]:.0f} DNN {pn["dnn"]:.0f} cut {pn["cut"]:.0f} cutoff {cutoff:.0f}')
say(f'LP ({solver} eps {eps:g}): primal {r["obj"] + s:.0f} valid dual bound {base:.0f}; wall {r["wall"]:.0f} s cpu {r["cpu"]:.0f} s iters {r["iters"]}')
say(f'LP short of cutoff {cutoff - base:.0f} ({(cutoff - base) / cutoff * 100:.2f}%); QPB short {cutoff - pn["qpb"]:.0f}')
gain = cb - base
say(f'child gains over the LP (cb - LP), all {gain.size} children: median {np.median(gain):.0f} p90 {np.percentile(gain, 90):.0f} max {gain.max():.0f}; children with cb >= cutoff: {(cb >= cutoff).sum()}')
ri = {v: k for k, v in enumerate(I)}; ci = {v: k for k, v in enumerate(J)}
# T1 over all bounded children of all truth lines
seen = {}
for ln, ch in P.items():
    for fx, v in ch: seen[fx] = v
fx = list(seen); r_ = [ri[int(f.split(':')[0])] for f in fx]; c_ = [ci[int(f.split(':')[1])] for f in fx]
cbv = cb[r_, c_]; dnn = np.array([seen[f]['dnn'] if seen[f]['dnn'] is not None else np.nan for f in fx])
qp = np.array([seen[f]['qpb'] for f in fx]); br = np.array([seen[f]['out'] == 'B' for f in fx])
ok = ~np.isnan(dnn)
say(f'T1, {ok.sum()} bounded children: Spearman with the child DNN: LP child bound {spearmanr(cbv[ok], dnn[ok])[0]:.3f}, child QPB {spearmanr(qp[ok], dnn[ok])[0]:.3f}')
def auc(score, lab):   # P(score of a non-branching child > score of a branching child)
    a = score[~lab]; b = score[lab]
    if len(a) == 0 or len(b) == 0: return float('nan')
    return float(((a[:, None] > b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()) / (len(a) * len(b)))
say(f'T1 AUC (non-branching above branching): LP child bound {auc(cbv, br):.3f}, child QPB {auc(qp, br):.3f}, child DNN {auc(np.nan_to_num(dnn, nan=1e18), br):.3f}')
for ln, ch in P.items():
    rr = [ri[int(f.split(':')[0])] for f, _ in ch]; cc = [ci[int(f.split(':')[1])] for f, _ in ch]
    d = np.array([v['dnn'] for _, v in ch], dtype=float); q = np.array([v['qpb'] for _, v in ch])
    say(f'  line {ln}: B={sum(v["out"] == "B" for _, v in ch)}  Spearman(cb,DNN) {spearmanr(cb[rr, cc], d)[0]:.3f}  Spearman(QPB,DNN) {spearmanr(q, d)[0]:.3f}')
# T2: every free line
rows = []
for k in range(lp.m):
    for kind, vec in (('r', cb[k, :]), ('c', cb[:, k])):
        lab = f'{kind}{I[k] if kind == "r" else J[k]}'
        g = vec - base
        rows.append((lab, int((vec >= cutoff).sum()), g.sum(), np.sort(g)[:3].sum(), np.median(g), lab in P))
rows.sort(key=lambda t: (-t[1], -t[2]))
with open(f'out/t2_{tag}_lines.tsv', 'w') as f:
    f.write('rank\tline\tcount_fathomed\tsum_gain\tsum_3_smallest_gain\tmedian_gain\ttruth_B\n')
    for i, t in enumerate(rows):
        tb = sum(v['out'] == 'B' for _, v in P[t[0]]) if t[5] else ''
        f.write(f'{i + 1}\t{t[0]}\t{t[1]}\t{t[2]:.0f}\t{t[3]:.0f}\t{t[4]:.0f}\t{tb}\n')
say('T2, rule-2 analogue (count, then sum of gains): rank of the truth lines among', len(rows))
for i, t in enumerate(rows):
    if t[5]:
        say(f'  {t[0]}: rank {i + 1}  count {t[1]}  sum gain {t[2]:.0f}  3 smallest {t[3]:.0f}  B(truth) {sum(v["out"] == "B" for _, v in P[t[0]])}')
alt = sorted(rows, key=lambda t: -t[3])
say('  by the sum of the 3 smallest gains: ' + ', '.join(f'{t[0]} #{i + 1}' for i, t in enumerate(alt) if t[5]))
