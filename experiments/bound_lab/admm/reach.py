#!/usr/bin/env python3
"""reach.py base.log var.log [plateau] [k1,k2,...]: for the baseline's certified running best at
iterations k, the first check at which the variant reaches it (linear interpolation between checks),
and the ratio k / that. With a plateau, also the baseline's remaining gap to it at k."""
import sys, re
def load(f):
    d = []
    for l in open(f):
        m = re.match(r'it\s+(\d+) lb ([-\d.e+]+)', l)
        if m and ' extra' not in l[-8:]: d.append((int(m.group(1)), float(m.group(2))))
    # the first run in the file only (dnntest -K logs hold several)
    out = []
    for it, v in d:
        if out and it <= out[-1][0]: break
        out.append((it, v))
    return out
b, v = load(sys.argv[1]), load(sys.argv[2])
plat = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != '-' else None
ks = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [200, 400, 600, 1000, 1500]
bd = dict(b)
print('| baseline at k | value | gap to plateau | variant reaches it at | ratio | variant at k (delta) |')
print('|---|---|---|---|---|---|')
for k in ks:
    if k not in bd: continue
    V = bd[k]; at = None
    for (i0, v0), (i1, v1) in zip([(0, -1e300)] + v, v):
        if v1 >= V:
            at = i1 if v0 == -1e300 or v1 == v0 else i0 + (i1 - i0) * (V - v0) / (v1 - v0)
            break
    vd = dict(v)
    print('| %d | %.0f | %s | %s | %s | %s |' % (k, V, '%.0f' % (plat - V) if plat else '-',
          '%.0f' % at if at else '> %d' % v[-1][0], '%.2fx' % (k / at) if at else '-',
          '%.0f (%+.0f)' % (vd[k], vd[k] - V) if k in vd else '-'))
