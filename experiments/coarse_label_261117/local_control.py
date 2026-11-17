#!/usr/bin/env python3
"""Screen candidates against the exact binary moment hull for every label union."""

import argparse
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

from separate import local_vertices, write_pool


def local_model(t, caps):
    k = len(caps)
    v, pairs = local_vertices(t, caps)
    center = v.mean(axis=0)
    _, sv, vt = np.linalg.svd(v - center, full_matrices=len(v) < len(center))
    orthogonal = vt[np.count_nonzero(sv > 1e-9):]
    all_binary = np.array(list(itertools.product(range(2), repeat=t)))
    union_specs = []
    total = len(center)
    for mask in range(1, 1 << (k - 1)):
        groups = [g for g in range(k) if mask & (1 << g)]
        left = sum(caps[g] for g in groups)
        binary = all_binary[(all_binary.sum(axis=1) <= left) & (t - all_binary.sum(axis=1) <= sum(caps) - left)]
        union_specs.append((groups, binary, slice(total, total + len(binary))))
        total += len(binary)

    equations, rhs = [], []
    for q in orthogonal:
        equations.append(np.r_[q, np.zeros(total - len(center))])
        rhs.append(q @ center)
    for groups, binary, indices in union_specs:
        row = np.zeros(total)
        row[indices] = 1
        equations.append(row)
        rhs.append(1)
        for i in range(t):
            row = np.zeros(total)
            row[i * k + np.array(groups)] = 1
            row[indices] = -binary[:, i]
            equations.append(row)
            rhs.append(0)
        for q, (i, j) in enumerate(pairs):
            row = np.zeros(total)
            for g, h in itertools.product(groups, repeat=2):
                row[t * k + q * k ** 2 + k * g + h] = 1
            row[indices] = -binary[:, i] * binary[:, j]
            equations.append(row)
            rhs.append(0)
    return center, pairs, csr_matrix(equations), rhs


def moment(p, k, t, pairs):
    mat = np.zeros((t * k + 1, t * k + 1))
    mat[0, 0] = 1
    mat[1:, 0] = mat[0, 1:] = p[:t * k]
    mat[1:, 1:] = np.diag(p[:t * k])
    for idx, (i, j) in enumerate(pairs):
        ii, jj = slice(1 + i * k, 1 + (i + 1) * k), slice(1 + j * k, 1 + (j + 1) * k)
        block = p[t * k + idx * k ** 2:t * k + (idx + 1) * k ** 2].reshape(k, k)
        mat[ii, jj], mat[jj, ii] = block, block.T
    return mat


def screen(cut, caps, model):
    center, pairs, a_eq, b_eq = model
    t, k, d = len(cut['members']), len(caps), len(center)
    coef = np.array(cut['coef'], dtype=float) / 1024
    rhs = cut['rhs'] / 1024
    objective = -np.r_[coef, np.zeros(a_eq.shape[1] - d)]
    result = linprog(objective, A_eq=a_eq, b_eq=b_eq, bounds=(0, 1), method='highs-ds')
    if not result.success:
        raise RuntimeError(result.message)
    upper = -result.fun - rhs
    # Mixing with the mean of feasible coarse vertices keeps every linear
    # constraint and can recover a PSD witness in the center's affine face.
    base = moment(center, k, t, pairs)
    e, u = np.linalg.eigh(base)
    basis = u[:, e > 1e-9] / np.sqrt(e[e > 1e-9])
    extreme = result.x[:d]
    relative = np.linalg.eigvalsh(basis.T @ moment(extreme, k, t, pairs) @ basis)[0]
    alpha = max(0., (1e-7 - relative) / (1 - relative)) if relative < 1e-7 else 0.
    witness = (1 - alpha) * extreme + alpha * center
    mineig = np.linalg.eigvalsh(moment(witness, k, t, pairs))[0]
    witness_v = coef @ witness - rhs
    return {'binary_union_lp_violation': upper, 'mixed_psd_violation': float(witness_v),
            'mixed_psd_min_eigenvalue': float(mineig),
            'binary_union_lp_max_eq_error': float(np.max(np.abs(result.eqlin.residual)))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('pool_json', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--filter-pool', type=Path, help='write directions that survive the union LP screen')
    args = p.parse_args()
    data = json.loads(args.pool_json.read_text())
    start = time.process_time()
    models, results = {}, []
    for rank, cut in enumerate(data['cuts']):
        t = len(cut['members'])
        if t not in models:
            models[t] = local_model(t, data['capacity'])
        results.append({'rank': rank, **screen(cut, data['capacity'], models[t])})
    def witnessed(result):
        return result['mixed_psd_violation'] > 1e-8 and result['mixed_psd_min_eigenvalue'] >= -1e-9

    kept = [r for r in results if r['binary_union_lp_violation'] > 1e-8]
    # Keep separation order inside each class, using local PSD witnesses first.
    kept.sort(key=lambda r: not witnessed(r))
    if args.filter_pool:
        side = data['side']
        fixes = [] if data['fixes'] == '-' else [tuple(map(int, p.split(':'))) for p in data['fixes'].split(',')]
        fixed = {p[side] for p in fixes}
        labels = np.array([g for a, g in enumerate(data['original_labels']) if a not in fixed])
        write_pool(args.filter_pool, side, labels, [data['cuts'][r['rank']] for r in kept])
    summary = {
        'cases': len(results),
        'locally_redundant_in_binary_lp': len(results) - len(kept),
        'positive_psd_witnesses': sum(witnessed(r) for r in results),
        'retained_ranks': [r['rank'] for r in kept],
        'cpu_s': time.process_time() - start,
    }
    args.output.write_text(json.dumps({**summary, 'cuts': results}, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
