#!/usr/bin/env python3
"""Modified wil50 instances for the objective-splitting test (Lagrangian
decomposition along a threshold subset K): A' = s*D - s*alpha*D_K, B = flow.
The DNN bound of a node of A' divided by s, plus alpha * OPT_K(node), is a
valid bound for the node whenever OPT_K(node) is a valid lower bound of the
K-part (here OPT_K comes from SA: a measurement of potential, not a proof).
Also D_row-only, D_col-only and D_K-only instances.
  python3 make_inst.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from coupled_split import thr_matrix
HERE = os.path.dirname(os.path.abspath(__file__))
n, A, B = gl.read_qaplib(gl.ROOT + '/qaplib/wil50.dat')
Dr, Dc = gl.grid_dist()
mid5 = sum(thr_matrix(nm) for nm in ['row<2', 'row<3', 'col<4', 'col<5', 'col<6'])
gl.write_qaplib(HERE + '/inst/wil50_row.dat', Dr, B)
gl.write_qaplib(HERE + '/inst/wil50_col.dat', Dc, B)
gl.write_qaplib(HERE + '/inst/wil50_mid5only.dat', mid5, B)
for s, k in ((4, 1), (2, 1), (4, 3)):   # alpha = k/s
    gl.write_qaplib(HERE + f'/inst/wil50_mid5_a{k}of{s}.dat', s * A - k * mid5, B)
print('ok')
