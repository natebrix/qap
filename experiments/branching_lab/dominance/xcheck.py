# cross-check onefixed_fast against Inst.onefixed (validated by brute force in validate.py) on random nodes of real instances
import numpy as np
from qapdata import read_dat
from dom import Inst, onefixed_fast
rng = np.random.default_rng(1); bad = 0; N = 0
for nm in ['tai35b', 'nug30', 'tho40', 'tai12b', 'chr12a']:
    inst = Inst(*read_dat(nm)); n = inst.n
    for t in range(40):
        d = int(rng.integers(1, n)); rows = rng.choice(n, d, replace=False); cols = rng.choice(n, d, replace=False)
        F = {int(r): int(c) for r, c in zip(rows, cols)}
        v1, _ = inst.onefixed(F); v2, _ = onefixed_fast(inst, F); N += 1
        if (v1 is None) != (v2 is None) or (v1 is not None and abs(v1 - v2) > 0.5): bad += 1; print(nm, d, v1, v2)
print('xcheck nodes', N, 'mismatches', bad)
