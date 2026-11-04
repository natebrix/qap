#!/usr/bin/env python3
"""For threshold subsets K: the K-part of a node's iterate against the SA
optimum of the K-part at the ROOT (no fixings: a single global problem whose
exact optimum would give a cut valid at every node) and at the node.
python3 root_vs_node.py SPEC YFILE"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridlib as gl
from coupled_split import thr_matrix, part_value
spec, yf = sys.argv[1], sys.argv[2]
nd = gl.Node(spec); Y = nd.load_Y(yf); Y4, _ = nd.blocks(Y); dg = np.einsum('aiai->ai', Y4)
root = gl.Node('-')
sets = [['row<3', 'col<5'], ['row<2', 'col<4'], ['row<2', 'row<3', 'col<5'], ['row<3', 'col<4', 'col<5', 'col<6'],
        ['row<2', 'row<3', 'col<4', 'col<5', 'col<6'], ['row<2', 'row<3', 'col<3', 'col<4', 'col<5', 'col<6', 'col<7'],
        ['row<1', 'row<2', 'row<3', 'row<4', 'col<3', 'col<4', 'col<5', 'col<6', 'col<7']]
print('K | iterate | node SA | root SA | iterate-node iterate-root')
for K in sets:
    D = sum(thr_matrix(nm) for nm in K)
    vy, nodeopt = part_value(nd, Y4, dg, D)
    _, rootopt = part_value(root, None, None, D)
    print(f'{"+".join(K):52s} | {vy:9.1f} | {nodeopt:7.0f} | {rootopt:7.0f} | {vy - nodeopt:7.1f} {vy - rootopt:7.1f}', flush=True)
