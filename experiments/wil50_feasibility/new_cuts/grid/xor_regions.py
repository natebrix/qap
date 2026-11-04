#!/usr/bin/env python3
"""Triangles on regions built from a row and a column threshold: R&C (corner
rectangle), R xor C (quadrant checkerboard), and the 12-super-cell partition's
2-colourings. python3 xor_regions.py SPEC YFILE"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from analyze_Y import tri_viol, region_pW, summarize
spec, yf = sys.argv[1], sys.argv[2]
nd = gl.Node(spec); Y = nd.load_Y(yf); Y4, _ = nd.blocks(Y)
rr = nd.loc_r; cc = nd.loc_c
print('region |S| max n_viol n>0.01 n>0.05 sum_pos')
for t in (1, 2, 3, 4):
    for s in (2, 3, 4, 5, 6, 7, 8):
        R = rr < t; C = cc < s
        for nm, ind in (('and', R & C), ('xor', R ^ C)):
            S = np.where(ind)[0]
            if not 0 < len(S) < nd.m: continue
            p, W = region_pW(nd, Y4, S); V, _ = tri_viol(p, W); sm = summarize(V)
            print(f'row<{t} {nm} col<{s}: {len(S):3d} {sm["max"]:.4f} {sm["n_viol"]:6d} {sm["n>0.01"]:6d} {sm["n>0.05"]:5d} {sm["sum_pos"]:8.3f}', flush=True)
