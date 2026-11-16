"""The wc sweep on single phases (the root of a subtree, cold, -d 3): the gap to the cutoff by phase
iteration, and the iteration at which each wc first reaches the gaps the base (wc=0.7) phase had at
300, 600 and its end.   python3 wcsweep.py"""
import os
from tracelib import load
H = os.path.dirname(os.path.abspath(__file__))
for title, runs in [("wil50 a220 root (254.0 short after the first solve)", [("0.7", "a220_base"), ("1.4", "a220r_wc14"), ("2.8", "a220r_wc28"), ("5.6", "a220r_wc56")]),
                    ("sko49 s232 root (245.5 short)", [("0.7", "s232_base"), ("1.4", "s232r_wc14"), ("2.8", "s232r_wc28"), ("5.6", "s232r_wc56")])]:
    hs = {}
    for wc, r in runs:
        p = os.path.join(H, "runs", r, "dnn.txt")
        if os.path.exists(p) and load(p): e = load(p)[0]; hs[wc] = [(t, e["cutoff"] - b) for t, b in e["hist"]]
    b = dict(hs["0.7"]); end = max(b)
    targets = [b[300], b[600], b[end]]
    print("\n%s\n| wc | gap at 100 | 200 | 300 | 500 | 700 | 900 | it to base's gap at 300 (%.1f) | at 600 (%.1f) | at its end, %d (%.1f) |" % (title, targets[0], targets[1], end, targets[2]))
    print("|---|---|---|---|---|---|---|---|---|---|")
    for wc, h in hs.items():
        g = dict(h)
        reach = [next((t for t, v in h if v <= tg + 1e-9), None) for tg in targets]
        print("| %s | %s | %s |" % (wc, " | ".join("%.1f" % g[t] if t in g else "-" for t in (100, 200, 300, 500, 700, 900)),
              " | ".join(("%d (%.2fx)" % (r, t0 / r)) if r else "not reached" for r, t0 in zip(reach, (300, 600, end)))))
