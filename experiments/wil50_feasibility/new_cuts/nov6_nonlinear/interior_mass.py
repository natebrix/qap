#!/usr/bin/env python3
"""Violation mass of the box inequalities on the INTERIOR lines only (rows 1-3,
columns 1-8) at a dumped iterate, for README.md section 3.
  python3 interior_mass.py SPEC Y.bin"""
import sys

import numpy as np

from nlc_analyze import Node, ROWS, COLS
from line_boxes import moments, separate_boxes

nd = Node(sys.argv[1]).load(sys.argv[2]); m = nd.m
vals = []
for r in range(1, ROWS - 1):
    cells = [c for c in range(m) if nd.r[c] == r]; p, W = moments(nd, cells)
    vals += list(separate_boxes(p, W, m, 40, 24, 7, seed=hash(f'row={r}') & 0xffff).values())
for c_ in range(1, COLS - 1):
    cells = [c for c in range(m) if nd.c[c] == c_]; p, W = moments(nd, cells)
    vals += list(separate_boxes(p, W, m, 40, 24, 7, seed=hash(f'col={c_}') & 0xffff).values())
vals = sorted(vals, reverse=True)
print(f'interior lines: violated {len(vals)}  top-200 sum {sum(vals[:200]):.3f}  max {vals[0]:.4f}'
      f'  median of top 200 {np.median(vals[:200]):.4f}')
