#!/usr/bin/env python3
"""L = mu * sum_k G_k over the local7 cuts dumped by local7_separate.py (each
<G_k, Y> <= beta_k), and the constant mu * sum beta_k to pass as -Ladd L.bin:CONST.
  python3 make_ladd7.py --json PREFIX.local7.json --mu MU --out L.bin [--top K]"""
import argparse, json
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--json', required=True); ap.add_argument('--mu', type=float, required=True)
ap.add_argument('--out', required=True); ap.add_argument('--top', type=int, default=0)
a = ap.parse_args()
d = json.load(open(a.json)); cuts = d['dumped'][:a.top] if a.top > 0 else d['dumped']
L = None; const = 0.0
for c in cuts:
    G = np.fromfile(c['file']); L = a.mu * G if L is None else L + a.mu * G; const += a.mu * c['beta']
L.tofile(a.out)
print(f"{len(cuts)} cuts, violations {cuts[-1]['violation']:.4f}..{cuts[0]['violation']:.4f}, mu {a.mu}, |L| {np.abs(L).sum():.3f}, "
      f"const {const:.6f} -> use: -Ladd {a.out}:{const:.9g}")
