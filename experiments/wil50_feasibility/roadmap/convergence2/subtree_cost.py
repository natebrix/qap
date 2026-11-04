#!/usr/bin/env python3
"""DNN CPU (first solve + cut phase) and iterations of the subtree under each
given node of a QAP_DNN_TRACE file (the node itself included).
usage: subtree_cost.py trace.txt "25:49 30:48 39:27 49:45" [...]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_breakdown import parse

rows = parse(sys.argv[1])
tot = sum(r['cpu'] + r.get('ccpu', 0) for r in rows)
for spec in sys.argv[2:]:
    s = set(spec.split())
    sub = [r for r in rows if s <= set(r['fixes'].split())]
    cpu = sum(r['cpu'] + r.get('ccpu', 0) for r in sub)
    it = sum(r['it'] + r.get('cit', 0) for r in sub)
    ph = sum(1 for r in sub if r['cut'])
    print(f"[{spec}]: {len(sub)} nodes, {ph} phases, {it} iterations, {cpu:.0f} CPU-s "
          f"({100 * cpu / tot:.0f}% of the trace's {tot:.0f})")
