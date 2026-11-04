#!/usr/bin/env python3
"""Warm vs cold siblings of the same parent in a QAP_DNN_TRACE file: in the
search a node is cold when its parent's DNN state did not survive a task or
unit boundary (state files carry the pool, not Y and Z). Siblings of the same
parent with the same pool status differ only in that.
usage: sibling_warm_cold.py trace.txt [trace.txt ...]
"""
import collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_breakdown import parse


def s(x, st):
    y = [r['it'] for r in x if r['status'] == st]
    return (f"{len(y):3d} @ {sum(y) / len(y):5.0f}" if y else '      -'), y


for name in sys.argv[1:]:
    rows = parse(name)
    g = collections.defaultdict(list)
    for r in rows:
        f = set(r['fixes'].split())
        par = None
        for p in rows:
            if p['depth'] == r['depth'] - 1 and set(p['fixes'].split()) <= f:
                par = p['fixes']
        g[(r['depth'], par, r['pool'])].append(r)
    print(name)
    tw, tc = [], []
    for k, v in sorted(g.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        w = [r for r in v if r['warm']]
        c = [r for r in v if not r['warm']]
        if w and c:
            sw, yw = s(w, 'fathomed'); sc, yc = s(c, 'fathomed')
            hw, _ = s(w, 'hopeless'); hc, _ = s(c, 'hopeless')
            print(f"  d{k[0]} pool {k[2]} parent [{k[1]}]: fathoms warm {sw}  cold {sc} | misses warm {hw}  cold {hc}")
            if yw and yc:
                tw.append(sum(yw) / len(yw)); tc.append(sum(yc) / len(yc))
    if tw:
        print(f"  mean over sibling groups: warm fathom {sum(tw)/len(tw):.0f} it, cold {sum(tc)/len(tc):.0f} it")
