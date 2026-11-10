#!/usr/bin/env python3
"""subtree.py A=trace.txt B=trace.txt: the root's children one by one, each with its whole
subtree (own first solve and phase plus every descendant), in both arms. A child's subtree is
complete when every node that branched under it has all its children in the trace."""
import re, sys
CUT = 48814.05
def load(fn):
    N = {}
    for ln in open(fn):
        m = re.match(r'pid \d+ depth (\d+) fixes (.*?) \| first lb ([\d.]+) primal [\d.\-]+ status (\w+) iters (\d+) cpu ([\d.]+) warm (\d) pool (\d+)', ln)
        if not m: continue
        fx = frozenset(m.group(2).split())
        ph = re.search(r'\| cut 1 best ([\d.]+) rounds (\d+) iters (\d+) cpu ([\d.]+) status (\w+)', ln)
        N[fx] = dict(d=int(m.group(1)), st=m.group(4), it=int(m.group(5)), cpu=float(m.group(6)), short=CUT - float(m.group(3)), fath=' fathomed 1' in ln,
                     ph=(CUT - float(ph.group(1)), int(ph.group(3)), float(ph.group(4)), ph.group(5)) if ph else None, kids=[])
    root = min(N, key=len)
    for fx, n in N.items():
        if fx == root: continue
        for f in fx:
            p = fx - {f}
            if p in N: N[p]['kids'].append(fx); n['new'] = f; break
    return N, root
def tot(N, fx):
    n = N[fx]
    it = n['it'] + (n['ph'][1] if n['ph'] else 0); cpu = n['cpu'] + (n['ph'][2] if n['ph'] else 0); cnt = 1; ok = True
    if not n['fath']:
        ok = len(n['kids']) == 50 - n['d']
        for k in n['kids']:
            a, b, c, o = tot(N, k); it += a; cpu += b; cnt += c; ok = ok and o
    return it, cpu, cnt, ok
arms = [(a.split('=')[0], *load(a.split('=')[1])) for a in sys.argv[1:]]
(an, A, ra), (bn, B, rb) = arms
for name, N, r in arms:
    n = N[r]
    print('%s root: first solve %.1f short, %d it%s; children in the trace %d of %d; nodes %d' % (name, n['short'], n['it'],
          ', phase to %.1f short in %d it (%s)' % (n['ph'][0], n['ph'][1], n['ph'][3]) if n['ph'] else ', no phase', len(n['kids']), 50 - n['d'], len(N)))
ka = {N['new']: fx for fx in A[ra]['kids'] for N in [A[fx]]}; kb = {N['new']: fx for fx in B[rb]['kids'] for N in [B[fx]]}
rows = []
for c in ka:
    if c not in kb: continue
    ta = tot(A, ka[c]); tb = tot(B, kb[c])
    if ta[3] and tb[3]: rows.append((c, A[ka[c]], B[kb[c]], ta, tb))
def kind(n): return 'first' if n['st'] == 'fathomed' else ('phase' if n['fath'] else 'branch')
print('children complete in both arms: %d' % len(rows))
for nm, idx in ((an, 3), (bn, 4)):
    ks = [r[idx - 2] for r in rows]
    print('  %s: fathom in the first solve %d, in their own phase %d, branch %d; own cost %d it (first %d + phase %d); subtrees %d it, %.0f s, %d nodes' % (
        nm, sum(kind(k) == 'first' for k in ks), sum(kind(k) == 'phase' for k in ks), sum(kind(k) == 'branch' for k in ks),
        sum(k['it'] + (k['ph'][1] if k['ph'] else 0) for k in ks), sum(k['it'] for k in ks), sum(k['ph'][1] for k in ks if k['ph']),
        sum(r[idx][0] for r in rows), sum(r[idx][1] for r in rows), sum(r[idx][2] for r in rows)))
sa = sum(r[3][0] for r in rows); sb = sum(r[4][0] for r in rows); ca = sum(r[3][1] for r in rows); cb = sum(r[4][1] for r in rows)
php = B[rb]['ph'][1] if B[rb]['ph'] else 0; phc = B[rb]['ph'][2] if B[rb]['ph'] else 0
print('  subtrees of these children: %s / %s = %.3fx in iterations, %.3fx in CPU; the parent phase of %s is %d it, %.0f s; scaled to all %d children it is %.3fx, %.3fx net of that phase' % (
      an, bn, sa / sb, ca / cb, bn, php, phc, 50 - A[ra]['d'], (sa * (50 - A[ra]['d']) / len(rows)) / (sb * (50 - A[ra]['d']) / len(rows) + php), (ca * (50 - A[ra]['d']) / len(rows)) / (cb * (50 - A[ra]['d']) / len(rows) + phc)))
if '-v' in sys.argv or True:
    for c, na, nb, ta, tb in sorted(rows, key=lambda r: -r[3][0]):
        def d(n, t): return '%-6s %5.1f short%s; subtree %5d it %4d nodes' % (kind(n), n['short'], (' -> %5.1f in %4d' % (n['ph'][0], n['ph'][1])) if n['ph'] else '                ', t[0], t[2])
        print('  %-6s | %s | %s' % (c, d(na, ta), d(nb, tb)))
