#!/usr/bin/env python3
"""Assemble the external objective term for dnntest/dnncuts -Ladd from
triplet_separate.py's cuts: L = mu * sum_k G_k over the selected cuts, where each
<G_k, Y> <= 0 is a valid triple-consistency inequality (right-hand side 0, so
no shift constant). Writes raw doubles, (m^2+1)^2 row-major.

  python3 make_ladd.py --prefix P --mu MU --out L.bin [--top K] [--min-viol V]
Reads P.triplet_*.G.bin and the matching .json (rounded_cut_violation)."""
import argparse, glob, json, os
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--prefix', required=True); ap.add_argument('--mu', type=float, required=True)
ap.add_argument('--out', required=True); ap.add_argument('--top', type=int, default=0)
ap.add_argument('--min-viol', type=float, default=1e-3)
a = ap.parse_args()
cuts = []
for j in sorted(glob.glob(a.prefix + '.triplet_*.json')):
    d = json.load(open(j))
    g = j[:-5] + '.G.bin'
    if d.get('success') and d['rounded_cut_violation'] > a.min_viol and os.path.exists(g):
        cuts.append((d['rounded_cut_violation'], g, d['facilities_local']))
cuts.sort(reverse=True)
if a.top > 0: cuts = cuts[:a.top]
if not cuts: raise SystemExit("no violated cuts")
L = None
for v, g, fac in cuts:
    G = np.fromfile(g)
    L = a.mu * G if L is None else L + a.mu * G
L.tofile(a.out)
n1 = int(round(len(L) ** 0.5))
print(f"{len(cuts)} cuts, violations {cuts[-1][0]:.4f}..{cuts[0][0]:.4f} (sum {sum(c[0] for c in cuts):.3f}), "
      f"mu {a.mu}, |L| {np.abs(L).sum():.3f}, nonzeros {np.count_nonzero(L)}, n1 {n1} -> {a.out}")
