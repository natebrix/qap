#!/usr/bin/env python3
"""ident_cmp.py [NAME...] : the switch-off identity. For each NAME the node traces (QAP_NODE_TRACE) and the DNN traces
(QAP_DNN_TRACE) of ident/NAME_base and ident/NAME_new are compared line by line, in order, after removing what is
CPU time: the node trace's last field, and the DNN trace's pid and every `cpu X` pair. Everything else (fixes, the
QPB, DNN and cut bounds to four decimals, the cutoff, the outcome, the rule, the line, the children; iterations,
statuses, pool sizes) must be equal. Prints one line per NAME with the line counts, a digest of the stripped node
trace, and the first difference if there is one. Exit 1 on any difference or unfinished run."""
import hashlib, os, re, sys
H = os.path.dirname(os.path.abspath(__file__))

def nodes(path):
    for l in open(path):
        yield l.rsplit(None, 1)[0]

def dnn(path):
    for l in open(path):
        yield re.sub(r'^pid \d+ ', '', re.sub(r'\b(cp)?cpu [-\d.e+]+', 'cpu _', l.rstrip('\n')))

def cmp(a, b, strip):
    n = 0; h = hashlib.sha256(); first = None
    ia, ib = strip(a), strip(b)
    while True:
        x, y = next(ia, None), next(ib, None)
        if x is None and y is None: break
        n += 1
        if x != y and first is None: first = (n, x, y)
        if x is not None: h.update((x + '\n').encode())
    return n, h.hexdigest()[:16], first

names = sys.argv[1:] or sorted({d.rsplit('_', 1)[0] for d in os.listdir(H + '/ident') if d.endswith(('_base', '_new'))})
bad = 0
for nm in names:
    A, B = (H + '/ident/%s_%s' % (nm, b) for b in ('base', 'new'))
    ex = []
    for d in (A, B):
        t = open(d + '/cmd.txt').read() if os.path.exists(d + '/cmd.txt') else ''
        ex.append(t.split('exit ')[-1].split()[0] if 'exit ' in t else 'running')
    if ex != ['0', '0']:
        print('%-16s NOT FINISHED (exit base %s, new %s)' % (nm, ex[0], ex[1])); bad = 1; continue
    out = []
    for f, strip in (('nodes.txt', nodes), ('dnn.txt', dnn)):
        if not os.path.exists(A + '/' + f) and not os.path.exists(B + '/' + f): out.append('%s: none' % f); continue
        n, dig, first = cmp(A + '/' + f, B + '/' + f, strip)
        if first:
            bad = 1; out.append('%s: DIFFER at line %d of %d\n    base: %s\n    new:  %s' % (f, first[0], n, first[1], first[2]))
        else:
            out.append('%s: %d lines identical (sha256 %s)' % (f, n, dig))
    tot = [next((l.strip() for l in open(d + '/log.txt') if l.strip().startswith('Nodes:')), '?') for d in (A, B)]
    print('%-16s %s; %s; statistics base "%s" new "%s"' % (nm, out[0], out[1], tot[0], tot[1]))
sys.exit(bad)
