"""Q1: how is the work distributed over branched nodes?  -> out/q1.txt

For every complete tree, every branched node v whose children are all in the trace:
  b(v)        children that branch
  excess(v)   b(v) - median b at v's depth (same run)
  direct(v)   nodes / CPU of bounding v's children (additive: every node but the root is exactly one node's child)
  under(v)    nodes / CPU strictly below v (nested: shares of a set of nodes are taken over the UNION of their subtrees)
Shares of the top 1/5/10% of branched nodes, ranked by b and by excess (ties broken by a hash of the node key,
never by the subtree), and the same within rule (9 / 2) and depth bands.
"""
import sys, statistics as st
from collections import defaultdict
from common import *

RUNS_Q1 = ['tai60b-261117', 'tai50b-261116', 'tai50b-r9', 'tai35b_base', 'tai40b_base', 'tho40ab-32', 'nug30-trace-a']
BANDS = [(0, 2), (3, 6), (7, 10), (11, 15), (16, 22), (23, 99)]


def band(d):
    for lo, hi in BANDS:
        if lo <= d <= hi:
            return f'{lo}-{hi}' if hi < 99 else f'{lo}+'


def analyse(run, out):
    nodes, tops = load(run)
    allN = len(nodes)
    allC = sum(n.cpu for n in nodes.values())
    br = [n for n in nodes.values() if n.out == 'B' and n.kids is not None and len(n.kidsl) == n.kids]
    for n in br:
        n.b = sum(1 for c in n.kidsl if c.out == 'B')
        n.dn = len(n.kidsl)
        n.dc = sum(c.cpu for c in n.kidsl)
    bydepth = defaultdict(list)
    for n in br:
        bydepth[n.depth].append(n.b)
    med = {d: st.median(v) for d, v in bydepth.items()}
    import zlib
    for n in br:
        n.ex = n.b - med[n.depth]
        n.tb = zlib.crc32(n.key.encode())   # ties broken by a hash of the node, not by its subtree
    order = sorted(nodes.values(), key=lambda x: x.depth)
    inst = RUNS[run]['inst']
    chain = set()
    try:
        obj, p = read_sln(inst)
        pm = None
        if inst == 'tai60b':
            pm = list(p); pm[43], pm[44] = p[44], p[43]
        for x in nodes.values():
            if all(p[f] == l for f, l in x.fixes) or (pm and all(pm[f] == l for f, l in x.fixes)):
                chain.add(id(x))
    except Exception:
        pass
    nchain = sum(1 for x in nodes.values() if id(x) in chain)

    def union_share(sel):
        """share of the tree in the subtrees OPENED by sel: each selected node's branching children that are
        not on the best known path (or tai60b's mirror), with everything below them (union, no double count)"""
        sel = set(id(x) for x in sel)
        flag = {}
        un = uc = 0
        for x in order:
            f = x.parent is not None and (flag.get(id(x.parent), False) or
                                          (id(x.parent) in sel and x.out == 'B' and id(x) not in chain))
            flag[id(x)] = f
            if f:
                un += 1; uc += x.cpu
        return un / allN, uc / allC

    import random
    rng = random.Random(1)
    brd = defaultdict(list)
    for n in br:
        brd[n.depth].append(n)

    def baseline(s, reps=20):
        """the same number of branched nodes per depth, drawn at random: mean opened shares"""
        need = defaultdict(int)
        for n in s:
            need[n.depth] += 1
        a = b2 = 0
        for _ in range(reps):
            pick = [y for d, k in need.items() for y in rng.sample(brd[d], k)]
            u, c = union_share(pick); a += u; b2 += c
        return a / reps, b2 / reps

    w = lambda s: out.write(s + '\n')
    w(f'\n## {run}: {allN} nodes, {allC / 3600:.1f} CPU h, {len(br)} branched nodes with all children traced; '
      f'sum b = {sum(n.b for n in br)}')
    # distribution of b
    bs = sorted((n.b for n in br), reverse=True)
    w(f'b: max {bs[0]}, top values {bs[:12]}, median {st.median(bs)}, mean {st.mean(bs):.2f}, b=0 at {sum(1 for x in bs if x == 0)}')
    w(f'chain (best known path{" and mirror" if inst == "tai60b" else ""}) nodes in the trace: {nchain}')
    w('| rank by | top | nodes | share of sum b | direct nodes | direct CPU | opened: nodes | opened: CPU | random at the same depths: nodes, CPU |')
    w('|---|---|---|---|---|---|---|---|---|')
    sb = sum(n.b for n in br); sdn = sum(n.dn for n in br); sdc = sum(n.dc for n in br)
    for key in ['b', 'excess']:
        f = (lambda n: (n.b, n.tb)) if key == 'b' else (lambda n: (n.ex, n.tb))
        rk = sorted(br, key=f, reverse=True)
        for pct in [1, 5, 10]:
            k = max(1, round(len(rk) * pct / 100))
            s = rk[:k]
            un, uc = union_share(s)
            rn, rc = baseline(s)
            thr = (s[-1].b if key == 'b' else s[-1].ex)
            ntie = sum(1 for n in br if (n.b if key == 'b' else n.ex) >= thr)
            w(f'| {key} | {pct}% | {k} (>= {thr:g}: {ntie}) | {sum(n.b for n in s) / sb:.1%} | {sum(n.dn for n in s) / allN:.1%} | '
              f'{sum(n.dc for n in s) / allC:.1%} | {un:.1%} | {uc:.1%} | {rn:.1%}, {rc:.1%} |')
    # by rule and depth band: direct shares, and within-band concentration of sum b
    w('\nby rule and depth band (direct cost = bounding the children of the band\'s branched nodes; "top 10%" = the '
      'band\'s top tenth by b, its share of the band\'s sum b and direct CPU)')
    w('| rule | depths | branched | sum b | mean b | var/mean b | max b | direct nodes | direct CPU | top 10%: k, share of b, of direct CPU | top 10% opened: nodes, CPU (random same depths) |')
    w('|---|---|---|---|---|---|---|---|---|---|---|')
    g = defaultdict(list)
    for n in br:
        g[(rule_of(run, n), band(n.depth))].append(n)
    for (r, bd) in sorted(g, key=lambda t: (int(t[1].split('-')[0].rstrip('+')), t[0])):
        s = g[(r, bd)]
        rk = sorted(s, key=lambda n: (n.b, n.tb), reverse=True)
        k = max(1, round(len(rk) / 10))
        tb = sum(n.b for n in s)
        mb = tb / len(s)
        vm = (st.pvariance([n.b for n in s]) / mb) if mb > 0 else float('nan')
        un, uc = union_share(rk[:k]); rn, rc = baseline(rk[:k], 10)
        w(f'| {r} | {bd} | {len(s)} | {tb} | {mb:.2f} | {vm:.2f} | {rk[0].b} | {sum(n.dn for n in s) / allN:.1%} | '
          f'{sum(n.dc for n in s) / allC:.1%} | {k}, {sum(n.b for n in rk[:k]) / max(tb, 1):.0%}, '
          f'{sum(n.dc for n in rk[:k]) / max(sum(n.dc for n in s), 1e-9):.0%} | {un:.1%}, {uc:.1%} ({rn:.1%}, {rc:.1%}) |')
    # the heaviest single branchings: b and what lies under
    return br, med


if __name__ == '__main__':
    with open('out/q1.txt', 'w') as out:
        out.write('# Q1: concentration of the work over branched nodes (q1_concentration.py)\n')
        for run in (sys.argv[1:] or RUNS_Q1):
            analyse(run, out)
    print(open('out/q1.txt').read())
