#!/usr/bin/env python3
"""est_by_level.py BOUND_LAB_TAI60B RESULTS: the tai60b estimate (the bound lab's E2) by tree level, beside the run's
per-level table, in the thesis's format (Table 4.1, "Actual vs. estimated performance").

BOUND_LAB_TAI60B is experiments/bound_lab/tai60b from claude/cool-babbage-ilnogj, with ../tai50b_size beside it
(git archive origin/claude/cool-babbage-ilnogj experiments/bound_lab/tai60b experiments/bound_lab/tai50b_size).
RESULTS is infra/phase1/runs/tai60b-261117/results.txt.

dive_est.py composes E2 by path depth D (the subtrees off the path at depth D). Each dive set's bnbdrv -e output
(est.out) prints, after every dive, the running mean of the dives so far by tree level; its block after the k-th
dive, times k, is the sum of the first k dives level by level. So a set's mean by level is the sum of those over
its seed directories over its dives, the same dives dive_est.py takes (the first min(dives, chains) of each); est.out
truncates each level's running mean of nodes to an integer, and the shortfall goes to the set's start level. The
root, the on-path nodes and E2's mirror copy go in as dive_est.py puts them: the mirror 44:58 at level 1 at the
on-path depth-1 node's cost, and everything from path depth 2 on twice."""
import sys, os, re, glob
L = os.path.abspath(sys.argv[1]); RES = sys.argv[2]
sys.path.insert(0, L + "/../tai50b_size")
import nodes as N
cpu = lambda d: d['fcpu'] + d['ccpu']

def blocks(path):
    """est.out's per-dive blocks: [{level: (nodes, time)}], one per finished dive, running means."""
    out, cur, inlev = [], None, False
    for ln in open(path):
        if 'Total Time:' in ln:
            cur = {}; out.append(cur); inlev = False
        elif ln.startswith('level') and 'nodes' in ln and '%fathom' in ln and '%c.elim' in ln and cur is not None:
            inlev = True
        elif inlev:
            f = ln.split()
            if len(f) == 5 and re.fullmatch(r'\d+', f[0]):
                cur[int(f[0])] = (float(f[1]), float(f[4]))
            else:
                inlev = False
    return out[1:]      # the first block is the estimator's start (no dive yet)

def set_by_level(setdir):
    """mean over the set's dives, by level, and the dive count; checked against dive_est.dives_of"""
    tot, n = {}, 0
    for sd in sorted(glob.glob(setdir + '/s*')):
        if not os.path.exists(sd + '/dives.txt') or not os.path.exists(sd + '/dnn.txt'): continue
        est = [ln.split() for ln in open(sd + '/dives.txt')]
        ch, cur = [], []
        for d in N.load([sd + '/dnn.txt']):
            if cur and d['depth'] <= cur[-1]['depth']: ch.append(cur); cur = []
            cur.append(d)
        if cur: ch.append(cur)
        k = min(len(est), len(ch))
        if not k: continue
        b = blocks(sd + '/est.out')
        assert len(b) >= k, (sd, len(b), k)
        blk = b[k - 1]
        # the block is the running mean of the first k dives: check it against dives.txt
        mn = sum(float(e[1]) for e in est[:k]) / k; mc = sum(float(e[2]) for e in est[:k]) / k
        sn = sum(v[0] for v in blk.values()); sc = sum(v[1] for v in blk.values())
        assert -0.01 <= mn - sn <= len(blk) and abs(sc - mc) <= 1 + 1e-3 * mc, (sd, sn, mn, sc, mc)
        # est.out prints the running mean of nodes truncated to an integer, level by level: the shortfall
        # against dives.txt (under one node a level) goes to the set's start level
        lo = min(blk); blk[lo] = (blk[lo][0] + mn - sn, blk[lo][1])
        for lv, (a, t) in blk.items():
            x = tot.setdefault(lv, [0.0, 0.0]); x[0] += k * a; x[1] += k * t
        n += k
    return ({lv: (a / n, t / n) for lv, (a, t) in tot.items()} if n else {}), n

def main():
    H = L
    root = N.load([H + '/runs/rootM/dnn.txt'])[0]
    path = {}
    for f in glob.glob(H + '/runs/d*path/dnn.txt'):
        if 'under' in f: continue
        d = N.load([f])[0]; path[d['depth']] = d
    size = {}       # start nodes by set, from dive_est.py's own output (the strata's state files are not committed)
    for ln in open(H + '/logs/dive_est.txt'):
        f = [x.strip() for x in ln.split('|')]
        if len(f) == 10 and re.fullmatch(r'\d+', f[1]):
            D = int(f[1]); size[f'near{D}'] = int(f[3].split(',')[0]); size[f'far{D}'] = int(f[4].split(',')[0])
    Ds = sorted(set(int(re.sub(r'\D', '', s)) for s in size if size[s]))
    sets = {}
    for sd in sorted(glob.glob(H + '/runs/far*')) + sorted(glob.glob(H + '/runs/near*')):
        m, n = set_by_level(sd)
        if n: sets[os.path.basename(sd)] = m
    est = {}                                    # level -> [nodes, CPU-s]
    def add(lv, a, t, w=1):
        x = est.setdefault(lv, [0.0, 0.0]); x[0] += w * a; x[1] += w * t
    add(0, 1, cpu(root))
    add(1, 1, cpu(path[1]))                     # E2's mirror 44:58, a copy of the on-path depth-1 child
    for D in Ds:
        w = 1 if D == 1 else 2                  # E2: everything from path depth 2 on twice (the mirror's subtree)
        if D in path: add(D, 1, cpu(path[D]), w)
        for k in ('far', 'near'):
            name = 'nearx1' if (D == 1 and k == 'near') else f'{k}{D}'
            if name in sets:
                for lv, (a, t) in sets[name].items(): add(lv, a, t, w)
            elif k == 'far' and size.get(name, 0):
                add(D, size[name], 0, w)        # a far set without a finished dive: its nodes only (as dive_est.py)
    # the run's per-level table
    act, inlev = {}, False
    for ln in open(RES):
        if ln.startswith('level') and '%fathom' in ln: inlev = True; continue
        if inlev:
            f = ln.split()
            if len(f) == 5 and re.fullmatch(r'\d+', f[0]): act[int(f[0])] = (int(f[1]), float(f[4]))
            elif act: break
    return act, est

HEAD = """tai60b: actual vs. estimated performance (the layout of thesis Table 4.1)

Actual: the run tai60b-261117 (Nov 17-18 2026) from the incumbent 608,215,055. The tree closed at
608,215,054, which proves that value optimal. Time: bnbdrv CPU seconds on AWS m7a cores.
Estimated: the estimate made before the run (the bound lab's E2, Nov 16). The path to the best known
solution was bounded node by node, the subtrees hanging off it were estimated by 267 Knuth dives from
the path's off-path children, and the depth-1 mirror 44:58 was counted as a copy of the path's subtree.
Time: CPU seconds on the 4-core cloud box the estimate ran on.
"""

if __name__ == '__main__':
    act, est = main()
    tn = sum(v[0] for v in est.values()); tt = sum(v[1] for v in est.values())
    an = sum(v[0] for v in act.values()); at = sum(v[1] for v in act.values())
    print(HEAD)
    print('               Actual                    Estimated')
    print('Level      Nodes          Time         Nodes          Time')
    for lv in range(max(max(act), max(est)) + 1):
        a = act.get(lv); e = est.get(lv)
        sa = f'{a[0]:10d} {a[1]:13.2f}' if a else f'{"-":>10} {"-":>13}'
        se = f'{e[0]:13.0f} {e[1]:13.2f}' if e and e[0] >= 0.5 else f'{"-":>13} {"-":>13}'
        print(f'{lv:5d} {sa} {se}')
    print(f'Total {an:10d} {at:13.2f} {tn:13.0f} {tt:13.2f}')
    sh = lambda d, lo, hi: 100 * sum(v[1] for l, v in d.items() if lo <= l <= hi) / sum(v[1] for v in d.values())
    print(f"""
- The run took {an / tn:.1f} times the nodes and {at / tt:.1f} times the time. The bound lab put the cloud box
  at about 0.70 of a c7a core's CPU per node (measured on tho40, not on tai60b; m7a and c7a are the same
  cores): in m7a seconds the estimate is about {tt / 0.70 / 1e6:.2f} million, and the run {at / (tt / 0.70):.1f} times that.
- Where the time went: levels 9-14 hold {sh(act, 9, 14):.0f}% of the run's time and levels 20-23 {sh(act, 20, 23):.0f}%; the
  estimate put {sh(est, 0, 8):.0f}% at levels 0-8, {sh(est, 9, 14):.0f}% at 9-14 and {sh(est, 20, 23):.0f}% at 20-23. Most of the run's time is
  inside the subtrees hanging off its path at depths 3 and 7, which reach levels 15 and 24 (68% of the
  CPU between them). A dive goes down until it picks a fathomed child, and most children are fathomed
  at every level: the dives from the path's children at depths 3 and 7 stopped by levels 6 and 8, and
  no dive went more than three levels below its start.
- The estimate also followed another tree: its path branched on row 12 at depth 2, where the run
  branched on row 1, and its dives ran the bound lab's build and parameters, not the run's.
- The run's results.txt gives a total time of 12,948,960.8 s; its levels sum to the figure above.
  The estimated levels are the dives' per-level output, composed as E2 was (dive_est.py); they sum to
  E2's 29,357 nodes and 1,003,936 CPU-s.""")
