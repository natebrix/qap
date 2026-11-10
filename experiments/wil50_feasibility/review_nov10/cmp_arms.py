#!/usr/bin/env python3
"""cmp_arms.py BASE=trace.txt ARM=trace.txt ...: each arm against the first, by kind of work
and on the leaves both arms fathom in the first solve (matched by fixes)."""
import re, sys
CUT = 48814.05
def parse(fn):
    out = {}
    for ln in open(fn):
        m = re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.]+) primal [\d.\-]+ status (\w+) iters (\d+) cpu ([\d.]+)', ln)
        if not m: continue
        d = dict(depth=int(m.group(1)), lb=float(m.group(3)), st=m.group(4), it=int(m.group(5)), cpu=float(m.group(6)), ph=None)
        mm = re.search(r'\| cut 1 best ([\d.]+) rounds (\d+) iters (\d+) cpu ([\d.]+) status (\w+)', ln)
        if mm: d['ph'] = (float(mm.group(1)), int(mm.group(3)), float(mm.group(4)), mm.group(5))
        d['fath'] = ' fathomed 1' in ln
        out[tuple(sorted(m.group(2).split()))] = d
    return out
arms = [(a.split('=')[0], parse(a.split('=')[1])) for a in sys.argv[1:]]
bname, base = arms[0]
def tot(t):
    leaves = [d for d in t.values() if d['st'] == 'fathomed']; miss = [d for d in t.values() if d['st'] != 'fathomed']; ph = [d for d in t.values() if d['ph']]
    return dict(n=len(t), fath=sum(d['fath'] for d in t.values()), leaf_it=sum(d['it'] for d in leaves), leaf_cpu=sum(d['cpu'] for d in leaves), nleaf=len(leaves),
                miss_it=sum(d['it'] for d in miss), miss_cpu=sum(d['cpu'] for d in miss), nmiss=len(miss),
                ph_it=sum(d['ph'][1] for d in ph), ph_cpu=sum(d['ph'][2] for d in ph), nph=len(ph), phf=sum(d['ph'][3] == 'fathomed' for d in ph))
bt = tot(base)
print('arm | nodes / fathomed | leaves: n, it, CPU | give-ups: n, it, CPU | phases: n (fathomed), it, CPU | total it, CPU | against %s: it, CPU' % bname)
for name, t in arms:
    x = tot(t); ti = x['leaf_it'] + x['miss_it'] + x['ph_it']; tc = x['leaf_cpu'] + x['miss_cpu'] + x['ph_cpu']
    bi = bt['leaf_it'] + bt['miss_it'] + bt['ph_it']; bc = bt['leaf_cpu'] + bt['miss_cpu'] + bt['ph_cpu']
    print('%s | %d / %d | %d, %d, %.0f | %d, %d, %.0f | %d (%d), %d, %.0f | %d, %.0f | %.3fx, %.3fx' % (name, x['n'], x['fath'], x['nleaf'], x['leaf_it'], x['leaf_cpu'],
          x['nmiss'], x['miss_it'], x['miss_cpu'], x['nph'], x['phf'], x['ph_it'], x['ph_cpu'], ti, tc, bi / ti, bc / tc))
    rootfx = set(min(t, key=len))
    for fx, d in sorted(t.items(), key=lambda kv: kv[1]['depth']):
        if d['ph']: print('     phase at depth %d (%s): first %.1f short, phase to %.1f short in %d it, %.0f s, %s' % (d['depth'], ' '.join(f for f in fx if f not in rootfx) or 'root', CUT - d['lb'], CUT - d['ph'][0], d['ph'][1], d['ph'][2], d['ph'][3]))
print()
for name, t in arms[1:]:
    both = [fx for fx in base if fx in t and base[fx]['st'] == 'fathomed' and t[fx]['st'] == 'fathomed']
    slow = [fx for fx in both if base[fx]['it'] >= 60]
    fast = [fx for fx in both if base[fx]['it'] < 60]
    s = lambda arm, keys, f: sum(arm[k][f] for k in keys)
    print('%s against %s on %d leaves in both: iterations %d / %d (%.3fx), CPU %.0f / %.0f (%.3fx); slow leaves (%d, >= 60 it in %s): %d / %d; the others (%d): %d / %d; nodes only in %s: %d, only in %s: %d' % (
        name, bname, len(both), s(base, both, 'it'), s(t, both, 'it'), s(base, both, 'it') / max(1, s(t, both, 'it')), s(base, both, 'cpu'), s(t, both, 'cpu'), s(base, both, 'cpu') / max(1e-9, s(t, both, 'cpu')),
        len(slow), bname, s(base, slow, 'it'), s(t, slow, 'it'), len(fast), s(base, fast, 'it'), s(t, fast, 'it'), bname, len([k for k in base if k not in t]), name, len([k for k in t if k not in base])))
