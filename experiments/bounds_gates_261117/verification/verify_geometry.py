#!/usr/bin/env python3
"""Check the CP normal and overlap formulas against explicit lifted matrices."""

import itertools

import numpy as np


def check_side(side, rng):
    m = 5
    n = 1 + m * m
    location, facility = np.exp(rng.normal(size=(2, m)))
    regions = [{0, 1, 2}, {1, 2, 4}, {0, 1}, {3}, {0, 1, 2, 3, 4}]
    perms = list(itertools.permutations(range(m)))
    checks = 0

    def lift(i, a):
        return 1 + (i * m + a if side == 1 else a * m + i)

    d = np.empty(n)
    d[0] = 0.63
    for i, a in itertools.product(range(m), repeat=2):
        d[lift(i, a)] = facility[i] * location[a]
    weights = 1 / np.outer(d, d) ** 2

    def p_normal(i, region):
        normal = np.zeros((n, n))
        for a in region:
            normal[lift(i, a), lift(i, a)] = 1
        return normal

    def w_normal(i, j, region):
        normal = np.zeros((n, n))
        for a, b in itertools.permutations(region, 2):
            r, c = lift(i, a), lift(j, b)
            normal[r, c] = normal[c, r] = 0.5
        return normal

    for s, t in itertools.product(regions, repeat=2):
        intersection = s & t
        u2 = sum(location[a] ** -2 for a in intersection)
        u4 = sum(location[a] ** -4 for a in intersection)
        for i, j in itertools.product(range(m), repeat=2):
            if i == j:
                actual = np.sum(p_normal(i, s) * p_normal(i, t) * weights)
                expected = facility[i] ** -4 * u4
            else:
                actual = np.sum(w_normal(i, j, s) * w_normal(i, j, t) * weights)
                expected = 0.5 * (facility[i] * facility[j]) ** -2 * (u2 * u2 - u4)
            assert np.isclose(actual, expected), (side, s, t, i, j, actual, expected)
            checks += 1

    assignment_supports = [[0] + [lift(i, a) for i, a in enumerate(p)] for p in perms]
    for s in regions:
        ps = [p_normal(i, s) for i in range(m)]
        cardinality = sum(ps)
        rows = [sum(w_normal(i, j, s) for j in range(m) if j != i) - (len(s) - 1) * ps[i]
                for i in range(m)]
        for support in assignment_supports:
            submatrix = np.ix_(support, support)
            assert cardinality[submatrix].sum() == len(s)
            assert all(row[submatrix].sum() == 0 for row in rows)
            checks += 1

    disjoint = [{0}, {1, 2}, {3, 4}]
    members = [0, 2, 3, 4]
    singletons = rng.integers(-1024, 1025, size=(4, 3))
    pairs = rng.integers(-1024, 1025, size=(4, 4, 3, 3))
    normal = np.zeros((n, n))
    norm_squared = 0
    for u, i in enumerate(members):
        for g, s in enumerate(disjoint):
            c = int(singletons[u, g]) / 1024
            norm_squared += c * c * len(s)
            for a in s:
                normal[lift(i, a), lift(i, a)] += c
            for v in range(u + 1, len(members)):
                for h, t in enumerate(disjoint):
                    w = int(pairs[u, v, g, h]) / 1024
                    norm_squared += 0.5 * w * w * (len(s) * len(t) - len(s & t))
                    for a, b in itertools.product(s, t):
                        if a != b:
                            r, c = lift(i, a), lift(members[v], b)
                            normal[r, c] += 0.5 * w
                            normal[c, r] += 0.5 * w
    assert norm_squared == np.sum(normal * normal)

    scores = []
    for labels in itertools.product(range(len(disjoint)), repeat=len(members)):
        if any(labels.count(k) > len(s) for k, s in enumerate(disjoint)):
            continue
        score = sum(singletons[u, g] for u, g in enumerate(labels))
        score += sum(pairs[u, v, labels[u], labels[v]]
                     for u, v in itertools.combinations(range(len(members)), 2))
        scores.append(score)
    rhs = max(scores) / 1024
    values = [normal[np.ix_(support, support)].sum() for support in assignment_supports]
    assert rhs == max(values)
    return checks + len(values)


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print(sum(check_side(side, rng) for side in (0, 1)), "lifted reference checks passed")
