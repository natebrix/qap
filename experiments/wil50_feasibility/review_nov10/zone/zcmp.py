#!/usr/bin/env python3
"""Compare two arms of the zone test on one parent: zcmp.py A/trace.txt B/trace.txt
(A: production zone, the parent branches without a phase; B: wide zone, the parent
runs a phase and its children inherit the pool). Children are matched by their fixes."""
import re, sys
CUT = 48814.05
def parse(fn):
    out = []
    for ln in open(fn):
        m = re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.]+) primal [\d.\-]+ status (\w+) iters (\d+) cpu ([\d.]+) warm (\d) pool (\d+)', ln)
        if not m: continue
        d = dict(depth=int(m.group(1)), fx=tuple(sorted(m.group(2).split())), lb=float(m.group(3)), st=m.group(4),
                 it=int(m.group(5)), cpu=float(m.group(6)), warm=int(m.group(7)), pool=int(m.group(8)), ph=None)
        mm = re.search(r'\| cut 1 best ([\d.]+) rounds (\d+) iters (\d+) cpu ([\d.]+) status (\w+)', ln)
        if mm: d['ph'] = (float(mm.group(1)), int(mm.group(3)), float(mm.group(4)), mm.group(5))
        out.append(d)
    return out
def summ(name, rows):
    root = rows[0]; kids = rows[1:]
    f = [k for k in kids if k['st'] == 'fathomed']; h = [k for k in kids if k['st'] != 'fathomed']
    sh = sorted(CUT - k['lb'] for k in h)
    cls = [sum(1 for s in sh if a <= s < b) for a, b in ((0, 50), (50, 130), (130, 244.08), (244.08, 1e9))]
    q = lambda p: sh[int(p * (len(sh) - 1))] if sh else float('nan')
    print('%s: parent first %.1f short after %d it (%.0f s)%s' % (name, CUT - root['lb'], root['it'], root['cpu'],
          '; phase to %.1f short in %d it (%.0f s, %s)' % (CUT - root['ph'][0], root['ph'][1], root['ph'][2], root['ph'][3]) if root['ph'] else '; no phase'))
    print('   children %d: fathomed %d (%.2f), first-solve iterations %d, CPU %.0f s (fathomed: %d it; misses: %d it)' % (
          len(kids), len(f), len(f) / max(1, len(kids)), sum(k['it'] for k in kids), sum(k['cpu'] for k in kids),
          sum(k['it'] for k in f), sum(k['it'] for k in h)))
    print('   misses %d: short q25/med/q75/max %.0f/%.0f/%.0f/%.0f; by class 0-50: %d, 50-130: %d, 130-244: %d, 244+: %d' % (
          len(h), q(.25), q(.5), q(.75), sh[-1] if sh else float('nan'), *cls))
    return root, kids
A = parse(sys.argv[1]); B = parse(sys.argv[2])
ra, ka = summ('A (no phase)', A); rb, kb = summ('B (phase)   ', B)
ma = {k['fx']: k for k in ka}; mb = {k['fx']: k for k in kb}
both = [fx for fx in ma if fx in mb]
print('matched children: %d' % len(both))
if both:
    fa = sum(ma[x]['st'] == 'fathomed' for x in both); fb = sum(mb[x]['st'] == 'fathomed' for x in both)
    ia = sum(ma[x]['it'] for x in both); ib = sum(mb[x]['it'] for x in both)
    ca = sum(ma[x]['cpu'] for x in both); cb = sum(mb[x]['cpu'] for x in both)
    print('   fathomed %d / %d; iterations %d / %d (%.2fx); CPU %.0f / %.0f s (%.2fx)' % (fa, fb, ia, ib, ia / max(1, ib), ca, cb, ca / max(1e-9, cb)))
    gain = sorted((mb[x]['lb'] - ma[x]['lb']) for x in both if ma[x]['st'] != 'fathomed' and mb[x]['st'] != 'fathomed')
    if gain: print('   both miss (%d): B - A first-solve bound q25/med/q75 %.0f/%.0f/%.0f' % (len(gain), gain[len(gain)//4], gain[len(gain)//2], gain[(3*len(gain))//4]))
    flips = sum(1 for x in both if ma[x]['st'] != 'fathomed' and mb[x]['st'] == 'fathomed'); back = sum(1 for x in both if ma[x]['st'] == 'fathomed' and mb[x]['st'] != 'fathomed')
    print('   miss in A, fathomed in B: %d; fathomed in A, miss in B: %d' % (flips, back))
