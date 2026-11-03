"""Regenerate the archived wil50 node-1, 1,000-triangle, mu=1 objective."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

root = Path(__file__).resolve().parents[2]
source = root / 'experiments/wil50_feasibility/stronger_bounds/astra'
sys.path.insert(0, str(source))
from region_cuts import cut_matrix

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
a = parser.parse_args()
spec = json.loads((source / 'tri1000_spec.json').read_text())
L = np.zeros((47*47+1, 47*47+1))
h = 0
for cut in spec['cuts']:
    G, rhs = cut_matrix(cut, 47)
    L += G
    h += rhs
assert len(spec['cuts']) == 1000 and h == 261
L.tofile(a.output)
print(json.dumps({'file': str(a.output), 'shift_add': -h,
                  'sha256': hashlib.sha256(a.output.read_bytes()).hexdigest()}))
