"""Verify the reported grids and exact threshold decompositions (no SDP solve)."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def validate(directory):
    results = []
    for name, height, width, axis in [
        ('nug12', 3, 4, 1), ('nug20', 4, 5, 1),
        ('nug24', 4, 6, 1), ('nug25', 5, 5, 1),
        ('nug30', 5, 6, 1), ('wil50', 5, 10, 0),
    ]:
        raw = (directory / (name + '.dat')).read_bytes()
        values = np.array([int(v) for v in raw.split()], dtype=np.int64)
        n = int(values[0])
        assert n == height * width and len(values) == 1 + 2*n*n, name
        matrices = values[1:].reshape(2, n, n)
        row, col = np.arange(n)//width, np.arange(n)%width
        distance = abs(row[:, None]-row) + abs(col[:, None]-col)
        masks = [row < k for k in range(1, height)]
        masks += [col < k for k in range(1, width)]
        decomposition = sum((s[:, None] != s).astype(np.int64) for s in masks)
        assert np.array_equal(matrices[axis], distance), name
        assert np.array_equal(distance, decomposition), name
        results.append(dict(instance=name, n=n, grid=[height, width],
                            distance_matrix='AB'[axis], thresholds=len(masks),
                            regional_m2_dimension=1+n+n*(n-1)//2,
                            dnn_dimension=1+(n-1)**2,
                            sha256=hashlib.sha256(raw).hexdigest()))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qaplib', type=Path,
                        default=Path(__file__).resolve().parents[2]/'qaplib')
    args = parser.parse_args()
    print(json.dumps(validate(args.qaplib), indent=2))
