#!/usr/bin/env python3
"""Offline replay of dnn.C's give-up rule on the Nov 6 warm set's per-check
histories (w = 0.25 logs), for several floors (par.hopeless_after).

The rule (dnn.C 2161-2196): at a regular check with it >= floor, with win = 5
checks (50 iterations), d1 = best[h-1] - best[h-1-win], d2 = best[h-1-win] -
best[h-1-2 win], r = min(d1/d2, 0.95) (0.95 if d2 <= 0 < d1, 0 if both <= 0):
give up when best + 5 d1 r / (1 - r) < cutoff.
"""
import re, sys, glob, os

CUT = 48814.05
ZONE = 244.08

def parse(path):
    nodes = []
    cur = None
    for line in open(path):
        if line.startswith('node 0 ') or line.startswith('CHILDNODE'):
            cur = {'name': line.split()[1] if line.startswith('CHILDNODE') else 'parent', 'hist': []}
            nodes.append(cur)
        m = re.match(r'it\s+(\d+) lb (\S+) ', line)
        if m and cur is not None:
            cur['hist'].append((int(m.group(1)), float(m.group(2))))
        m = re.match(r'RESULT node=(\d+) .*lb=(\S+) .*status=(\S+) iters=(\d+) cpu=(\S+)', line)
        if m and cur is not None:
            cur.update(node=int(m.group(1)), lb=float(m.group(2)), status=m.group(3),
                       iters=int(m.group(4)), cpu=float(m.group(5)))
    return nodes

def replay(hist, floor, win=5):
    """Return (iteration, best) where the rule first fires, or None."""
    best = []
    b = float('-inf')
    for it, lb in hist:
        b = max(b, lb)
        best.append((it, b))
        h = len(best)
        if b > CUT:
            return ('fathom', it, b)
        if it >= floor and h > 2 * win:
            d1 = best[h - 1][1] - best[h - 1 - win][1]
            d2 = best[h - 1 - win][1] - best[h - 1 - 2 * win][1]
            r = min(d1 / d2, 0.95) if d2 > 0 else (0.95 if d1 > 0 else 0.0)
            if b + 5 * d1 * r / (1 - r) < CUT:
                return ('giveup', it, b)
    return ('end', best[-1][0], best[-1][1])

def main():
    d = sys.argv[1]
    floors = [110, 150, 200, 300]
    rows = []
    for path in sorted(glob.glob(os.path.join(d, 'w025_*.log'))):
        fam = os.path.basename(path)[5:-4]
        for nd in parse(path):
            if 'status' not in nd:
                continue
            rows.append((fam, nd))
    print('family child status iters lb300(short) | ' + ' | '.join('floor %d: it, short' % f for f in floors))
    tot = {f: 0 for f in floors}
    lost = {f: 0 for f in floors}
    for fam, nd in rows:
        if nd['name'] == 'parent':
            continue
        out = []
        for f in floors:
            kind, it, b = replay(nd['hist'], f)
            tot[f] += it
            if nd['status'] == 'fathomed' and kind != 'fathom':
                lost[f] += 1
            out.append('%s %d, %.1f' % (kind[0], it, CUT - b))
        print('%s %s %s %d %.1f | %s' % (fam, nd['name'], nd['status'], nd['iters'], CUT - nd['lb'], ' | '.join(out)))
    print('total iterations by floor:', tot, 'fathoms lost:', lost)
    # the cold parents as well (a cold call's history)
    print('\ncold parents: where the rule would fire by floor')
    for fam, nd in rows:
        if nd['name'] != 'parent':
            continue
        out = []
        for f in floors:
            kind, it, b = replay(nd['hist'], f)
            out.append('%s %d, %.1f' % (kind[0], it, CUT - b))
        print('%s parent %s %d %.1f | %s' % (fam, nd['status'], nd['iters'], CUT - nd['lb'], ' | '.join(out)))

main()
