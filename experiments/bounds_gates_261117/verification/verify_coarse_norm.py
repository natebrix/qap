#!/usr/bin/env python3
"""Minimize lifted normal norms modulo local affine identities in a coarse pool."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np


def model(order, capacity):
    k = len(capacity)
    labels = np.array([
        config for config in itertools.product(range(k), repeat=order)
        if np.all(np.bincount(config, minlength=k) <= capacity)
    ])
    pairs = list(itertools.combinations(range(order), 2))
    pair_size = k * k
    vertices = np.zeros((len(labels), order * k + len(pairs) * pair_size), dtype=np.int64)
    rows = np.arange(len(labels))
    for i in range(order):
        vertices[rows, i * k + labels[:, i]] = 1
    for q, (i, j) in enumerate(pairs):
        offset = order * k + q * pair_size
        vertices[rows, offset + k * labels[:, i] + labels[:, j]] = 1

    differences = vertices - vertices[0]
    has_short_matrix = differences.shape[0] < differences.shape[1]
    _, singular, right = np.linalg.svd(differences, full_matrices=has_short_matrix)
    rank = np.count_nonzero(singular > 1e-8)
    affine_null = right[rank:].T

    # Off-diagonal lifted coefficients are shared by two symmetric entries.
    pair_weight = (np.outer(capacity, capacity) - np.diag(capacity)) / 2
    weights = np.r_[np.tile(capacity, order), np.tile(pair_weight.ravel(), len(pairs))]
    if np.min(weights) <= 0:
        raise ValueError("This audit requires a nonzero lifted weight for every feature")
    root_weight = np.sqrt(weights)
    weighted_null_inverse = np.linalg.pinv(root_weight[:, None] * affine_null)
    return vertices, affine_null, weights, root_weight, weighted_null_inverse


def compare(cut, basis):
    vertices, affine_null, weights, root_weight, weighted_null_inverse = basis
    initial = np.asarray(cut["coef"], dtype=np.int64)
    initial_rhs = int(np.max(vertices @ initial))
    if initial_rhs != cut["rhs"]:
        raise ValueError("Saved RHS differs from the enumerated maximum")
    minimum = initial - affine_null @ (weighted_null_inverse @ (root_weight * initial))
    translated_rhs = initial_rhs - vertices[0] @ (initial - minimum)
    difference = (vertices @ minimum - translated_rhs) - (vertices @ initial - initial_rhs)
    translation_error = float(np.max(np.abs(difference)))
    orthogonality_error = float(np.max(np.abs(affine_null.T @ (weights * minimum))))
    if translation_error > 1e-9 or orthogonality_error > 1e-8:
        raise ValueError("Minimum norm or affine equivalence check failed")

    # The loader accepts dyadic numerators in [-1024, 1024].
    multiplier = min(1.0, 1024 / np.max(np.abs(minimum)))
    roundtrip = np.rint(multiplier * minimum).astype(np.int64)
    rounded_rhs = int(np.max(vertices @ roundtrip))
    initial_residual = vertices @ initial - initial_rhs
    rounded_residual = (vertices @ roundtrip - rounded_rhs) / multiplier
    original_norm = np.linalg.norm(root_weight * initial)
    norm = np.linalg.norm(root_weight * minimum)
    rounded_norm = np.linalg.norm(root_weight * roundtrip) / multiplier
    if abs(original_norm / 1024 - cut["lift_norm"]) > 1e-10:
        raise ValueError("Saved lifted norm differs from the diagonal-weight formula")
    return {
        "order": len(cut["members"]),
        "original_over_minimum_norm": float(original_norm / norm),
        "original_over_dyadic_norm": float(original_norm / rounded_norm),
        "dyadic_scale": float(multiplier),
        "dyadic_rhs_numerator": rounded_rhs,
        "max_unrounded_residual_error_numerator_units": translation_error,
        "max_rounded_residual_error": float(np.max(np.abs(rounded_residual - initial_residual)) / 1024),
    }


def quantiles(values):
    keys = ("min", "q1", "median", "q3", "max")
    return dict(zip(keys, np.quantile(values, [0, 0.25, 0.5, 0.75, 1]).tolist()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pool_json", type=Path, help="coarse separator's JSON output")
    args = parser.parse_args()
    raw = args.pool_json.read_bytes()
    source = json.loads(raw)
    capacity = np.asarray(source["capacity"])
    models = {}
    comparisons = []
    for cut in source["cuts"]:
        order = len(cut["members"])
        if order not in models:
            models[order] = model(order, capacity)
        comparisons.append(compare(cut, models[order]))
    report = {
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "capacity": capacity.tolist(),
        "order_counts": {order: sum(c["order"] == order for c in comparisons) for order in models},
        "original_over_minimum_norm": quantiles([c["original_over_minimum_norm"] for c in comparisons]),
        "original_over_dyadic_norm": quantiles([c["original_over_dyadic_norm"] for c in comparisons]),
        "max_unrounded_residual_error_numerator_units": max(
            c["max_unrounded_residual_error_numerator_units"] for c in comparisons
        ),
        "max_rounded_residual_error": max(c["max_rounded_residual_error"] for c in comparisons),
        "rows": comparisons,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
