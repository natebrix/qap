#!/usr/bin/env python3
"""Replay dnn.C's hopeless rule on a -v trace: at every check from iteration
`start`, the projection best + 5 d1 r / (1 - r) (d1, d2: gains of the running
best over the last two windows of 5 checks, r = min(d1/d2, 0.95)) against the
cutoff. Prints the closest call (projection - cutoff) and the smallest
give-up threshold H that never fires before the node fathoms.

  python3 replay.py cutoff trace.log [start]
"""
import re, sys
cut = float(sys.argv[1]); f = sys.argv[2]; start = int(sys.argv[3]) if len(sys.argv) > 3 else 200
pts = [(int(m.group(1)), float(m.group(2))) for m in
       re.finditer(r"^it\s+(\d+) lb (\S+)", open(f).read(), re.M)]
best = []; b = -float("inf")
for it, lb in pts:
    b = max(b, lb); best.append((it, b))
w = 5
proj = {}
for h in range(2 * w + 1, len(best) + 1):
    it, bb = best[h - 1]
    d1 = bb - best[h - 1 - w][1]; d2 = best[h - 1 - w][1] - best[h - 1 - 2 * w][1]
    r = min(d1 / d2, 0.95) if d2 > 0 else (0.95 if d1 > 0 else 0.0)
    proj[it] = bb + 5 * d1 * r / (1 - r) if r < 1 else float("inf")
fath = next((it for it, bb in best if bb > cut), None)
fires = [it for it, p in proj.items() if p < cut and (fath is None or it < fath)]
safe_H = (max(fires) + 10) if fires else min(proj)
close = min(((p - cut, it) for it, p in proj.items() if it >= start and (fath is None or it < fath)), default=None)
print(f"{f}: fathoms at {fath}; closest call from {start}: {close[0]:.1f} at it {close[1]}" if close else
      f"{f}: fathoms at {fath}; no checks between {start} and the fathom")
print(f"  rule would fire at iterations {fires[:8]}{'...' if len(fires) > 8 else ''}; smallest safe H = {safe_H}")
