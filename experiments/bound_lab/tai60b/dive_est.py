#!/usr/bin/env python3
"""dive_est.py: the tai60b estimate from the dives (runs/farD, runs/nearD: Knuth dives from the off-path children
at depth D of the on-path node at depth D-1, split by strata.py into a near and a far stratum), the measured
on-path chain (runs/dD_*path) and the root (runs/rootM). Each finished dive (one line of dives.txt) is an
unbiased estimate of the nodes and CPU-s of all subtrees rooted at its set's nodes; a set's estimate is the mean
over its dives. Iterations: the dive's CPU estimate times (iterations / CPU-s) of the dive's own nodes [e].
Total = root + sum over D of (on-path node + far_D + near_D). Interval: bootstrap over dives within sets,
5th-95th percentile. Also the model check: every dive node against the survive rule (strata.excess, s(d), f(d))."""
import sys, os, re, glob, random, statistics as st
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/../tai50b_size'); sys.path.insert(0, H)
import nodes as N
import strata as S
best = S.best
cpu = lambda d: d['fcpu'] + d['ccpu']; its = lambda d: d['fit'] + d['cit']
def dives_of(setdir):
    out = []
    for sd in sorted(glob.glob(setdir + '/s*')):
        if not os.path.exists(sd + '/dives.txt') or not os.path.exists(sd + '/dnn.txt'): continue
        est = [ln.split() for ln in open(sd + '/dives.txt')]
        ch, cur = [], []
        for d in N.load([sd + '/dnn.txt']):
            if cur and d['depth'] <= cur[-1]['depth']: ch.append(cur); cur = []
            cur.append(d)
        if cur: ch.append(cur)
        for e, c in zip(est, ch):        # finished dives only (a dive in flight has trace lines but no estimate)
            out.append(dict(nodes=float(e[1]), cpu=float(e[2]), deep=int(e[0]), chain=c,
                            its=float(e[2]) * sum(its(d) for d in c) / max(1e-9, sum(cpu(d) for d in c))))
    return out
def main():
    root = N.load([H + '/runs/rootM/dnn.txt'])[0]
    path = {}
    for f in glob.glob(H + '/runs/d*path/dnn.txt'):
        if 'under' in f: continue
        d = N.load([f])[0]; path[d['depth']] = d
    sets = {}
    for sd in sorted(glob.glob(H + '/runs/far*')) + sorted(glob.glob(H + '/runs/near*')):
        name = os.path.basename(sd); dv = dives_of(sd)
        if dv: sets[name] = dv
    size = {}
    for f in glob.glob(H + '/runs/strata/off*_*.state'):
        k = os.path.basename(f)[3:-6].split('_'); size[k[1] + k[0]] = sum(1 for ln in open(f) if ln.startswith('node '))
    levels = sorted(set(int(re.sub(r'\D', '', s)) for s in size))
    print('## dive sets\n\n| set | start nodes | dives | node bounds: first-solve fathom / phase fathom / branched | deepest dive | estimate: nodes | iterations | CPU-s |')
    print('|---|---|---|---|---|---|---|---|')
    key = lambda s: (int(re.sub(r'\D', '', s)), s)
    for s in sorted(sets, key=key):
        dv = sets[s]; nb = [d for x in dv for d in x['chain']]
        se = lambda k: (st.stdev([x[k] for x in dv]) / len(dv) ** .5) if len(dv) > 1 else float('nan')
        print(f"| `{s}` | {size.get(s, '?')} | {len(dv)} | {sum(d['out'] == 'D' for d in nb)} / {sum(d['out'] == 'C' for d in nb)} / {sum(d['out'] == 'B' for d in nb)} | "
              f"{max(x['deep'] for x in dv)} | {st.mean(x['nodes'] for x in dv):,.0f} ± {se('nodes'):,.0f} | {st.mean(x['its'] for x in dv):,.0f} | {st.mean(x['cpu'] for x in dv):,.0f} ± {se('cpu'):,.0f} |")
    def total(sets_):
        t = {'nodes': 1.0, 'its': float(its(root)), 'cpu': cpu(root)}; by = {0: dict(t)}
        for D in levels:
            b = by.setdefault(D, {'nodes': 0.0, 'its': 0.0, 'cpu': 0.0})
            if D in path: b['nodes'] += 1; b['its'] += its(path[D]); b['cpu'] += cpu(path[D])
            for k in ('far', 'near'):
                dv = sets_.get(f'{k}{D}')
                if dv:
                    for q in ('nodes', 'its', 'cpu'): b[q] += st.mean(x[q] for x in dv)
                elif size.get(f'{k}{D}', 0) and k == 'far':      # a far set without dives: its nodes at 45 iterations [e]
                    b['nodes'] += size[f'{k}{D}']
        tot = {q: sum(b[q] for b in by.values()) for q in ('nodes', 'its', 'cpu')}
        return tot, by
    tot, by = total(sets)
    rng = random.Random(1); reps = []
    for _ in range(3000):
        reps.append(total({s: [rng.choice(dv) for _ in dv] for s, dv in sets.items()})[0])
    def iv(q):
        v = sorted(r[q] for r in reps); return v[int(.05 * (len(v) - 1))], v[int(.95 * (len(v) - 1))]
    print(f"\n## estimate (E1: the sampling design alone)\n")
    for q, u in (('nodes', 'nodes'), ('its', 'ADMM iterations'), ('cpu', 'CPU-s')):
        lo, hi = iv(q)
        print(f"- {tot[q]:,.0f} {u} (bootstrap 5%-95%: {lo:,.0f} to {hi:,.0f})" + (f" = {tot[q] / 3600:,.0f} core-hours ({lo / 3600:,.0f} to {hi / 3600:,.0f})" if q == 'cpu' else ''))
    missing = [f'{k}{D}' for D in levels for k in ('far', 'near') if size.get(f'{k}{D}', 0) and f'{k}{D}' not in sets]
    print('- sets with start nodes and no finished dive (near: counted as 0; far: their nodes only):', ' '.join(missing) or 'none')
    print('\n## by level (the on-path node at depth D and the subtrees hanging off the path at depth D)\n\n| D | on-path node: iterations, CPU-s | near: start nodes, dives | far: start nodes, dives | nodes | iterations | CPU-s | share of CPU |\n|---|---|---|---|---|---|---|---|')
    for D in sorted(by):
        b = by[D]; p = path.get(D)
        print(f"| {D} | {('%d, %.0f' % (its(p), cpu(p))) if p else ('root: %d, %.0f' % (its(root), cpu(root)) if D == 0 else '-')} | {size.get(f'near{D}', 0)}, {len(sets.get(f'near{D}', []))} | {size.get(f'far{D}', 0)}, {len(sets.get(f'far{D}', []))} | {b['nodes']:,.0f} | {b['its']:,.0f} | {b['cpu']:,.0f} | {100 * b['cpu'] / tot['cpu']:.1f}% |")
    # the model check, node by node
    sl = {d: best - p['final'] for d, p in path.items()}; fl = {d: best - p['lb'] for d, p in path.items()}
    seen = {}
    for dv in sets.values():
        for x in dv:
            for d in x['chain']: seen[tuple(sorted(d['fixes']))] = d
    rows = []
    for fx, d in seen.items():
        fix = {int(a.split(':')[0]): int(a.split(':')[1]) for a in fx}
        if all(S.p0[r] == v for r, v in fix.items()): continue
        D = d['depth']; s_ = sl.get(D, sl[max(sl)]); ex = S.excess(fix)
        rows.append((ex / s_, d['out'], its(d), D, ','.join(fx), ex, (d['cutoff'] - d['final'])))
    print('\n## model check: distinct off-path dive nodes by excess / s(depth)\n\n| excess / s(depth) | nodes | branched | fathomed by the phase | fathomed by the first solve | mean iterations of the fathomed |\n|---|---|---|---|---|---|')
    for lo, hi in ((0, .5), (.5, 1), (1, 1.6), (1.6, 3), (3, 5), (5, 1e9)):
        R = [r for r in rows if lo <= r[0] < hi]; F = [r[2] for r in R if r[1] != 'B']
        print(f"| {lo} to {hi if hi < 1e9 else 'inf'} | {len(R)} | {sum(r[1] == 'B' for r in R)} | {sum(r[1] == 'C' for r in R)} | {sum(r[1] == 'D' for r in R)} | {st.mean(F) if F else float('nan'):.0f} |")
    ok = sum((r[0] < 1) == (r[1] == 'B') for r in rows)
    print(f"\nrule 'branches iff excess < s(depth)': right on {ok} of {len(rows)} nodes; wrong on:")
    for r in sorted(rows):
        if (r[0] < 1) != (r[1] == 'B'): print(f"  depth {r[3]} `{r[4]}` excess {r[5]:,.0f} = {r[0]:.2f} s(d): {'branched' if r[1] == 'B' else 'fathomed'} ({r[2]} iterations, end {r[6]:,.0f} short)")
    retries = sum(open(f).read().count('retry:') for f in glob.glob(H + '/runs/*/s*/err.txt'))
    print(f"\npool/cold retries in the dives (pool_retry.patch): {retries}")
if __name__ == '__main__': main()
