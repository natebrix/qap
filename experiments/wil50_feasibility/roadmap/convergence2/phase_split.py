#!/usr/bin/env python3
"""Split the search's cut phases (QAP_DNN_TRACE) into Polyak rounds and the
last solve; a last solve longer than one round (300) is the quiet round.
usage: phase_split.py trace.txt [round_iters=300]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_breakdown import parse

rows = parse(sys.argv[1])
ri = int(sys.argv[2]) if len(sys.argv) > 2 else 300
tot = quiet = uncut = 0
for r in rows:
    if not r['cut']:
        continue
    last = r['cit'] - ri * (r['crounds'] - 1)
    q = last if last > ri else 0
    tot += r['cit']; quiet += q; uncut += r['it']
    print(f"d{r['depth']} {r['cutoff'] - r['lb']:6.1f} short after {r['it']:4d} uncut it | rounds {r['crounds']} "
          f"phase it {r['cit']:4d} (last solve {last:4d}{', quiet' if q else ''}) -> {r['cstatus']}")
print(f"phase iterations {tot}, of which in quiet rounds {quiet} ({100 * quiet / tot:.0f}%); "
      f"uncut iterations before these phases {uncut}")
