#!/usr/bin/env python3
"""Check coarse cuts and CP certificates against exhaustive small-node solutions."""

import argparse
import itertools
import math
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile


def synthetic_pool(m, side, rng):
    k = 3
    labels = [a % k for a in range(m)]
    cuts = []
    for i, j in ((0, 1), (2, 3), (4, 5)):
        for sign in (1, -1):
            cuts.append(([i, j], [0] * (2 * k) + [sign * 1024] * k ** 2))
            for g in range(k):
                single = [0] * (2 * k)
                single[g] = -sign * 1024
                pair = [sign * 1024 if a == g else 0 for a, _ in itertools.product(range(k), repeat=2)]
                cuts.append(([i, j], single + pair))
    for _ in range(8):
        subset = sorted(rng.sample(range(m), 4))
        cuts.append((subset, [rng.randint(-1024, 1024) for _ in range(4 * k + 6 * k ** 2)]))

    def evaluate(features, coefficients):
        value = sum(coefficients[i * k + a] for i, a in enumerate(features))
        start = len(features) * k
        for i, j in itertools.combinations(range(len(features)), 2):
            value += coefficients[start + features[i] * k + features[j]]
            start += k ** 2
        return value

    perms = list(itertools.permutations(range(m)))
    capacities = [labels.count(g) for g in range(k)]
    maxima = []
    for subset, coef in cuts:
        choices = [s for s in itertools.product(range(k), repeat=len(subset))
                   if all(s.count(g) <= capacities[g] for g in range(k))]
        maximum = max(evaluate(choice, coef) for choice in choices)
        realized = max(evaluate([labels[p[i]] for i in subset], coef) for p in perms)
        assert maximum == realized
        maxima.append(maximum)
    assert min(maxima) < 0 < max(maxima)
    return (f"QAP_COARSE_POOL 1 {side} {m} {k} {len(cuts)}\n" + " ".join(map(str, labels)) + "\n"
            + "\n".join(" ".join(map(str, [len(t), *t, *c])) for t, c in cuts) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", default="experiments/dnn_cuts/dnncuts")
    args = parser.parse_args()
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    worst = -math.inf
    lazy = checked = cases = 0
    n = 7
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        recipe = work / "test.par"
        recipe.write_text("[settings]\ncut.cp_lazycert = 1\ndnn.beta_sched = w=0.36,wc=2.88\n")
        for sample in range(4):
            rng = random.Random(10519 + sample)
            a = [[rng.randint(-20, 31) for _ in range(n)] for _ in range(n)]
            b = [[rng.randint(-9, 43) for _ in range(n)] for _ in range(n)]
            instance = work / f"sample{sample}.dat"
            instance.write_text(f"{n}\n" + "\n".join(" ".join(map(str, row)) for row in a + b) + "\n")
            fixes = [(0, 4)] if sample % 2 else []
            m = n - len(fixes)
            feasible = (p for p in itertools.permutations(range(n)) if all(p[i] == j for i, j in fixes))
            optimum = min(sum(a[i][j] * b[p[i]][p[j]] for i in range(n) for j in range(n)) for p in feasible)
            for side in (0, 1):
                pool = work / "pool"
                pool.write_text(synthetic_pool(m, side, rng))
                options = ["-p", str(recipe), "-prefix", "-side", str(side), "-I", "25", "-k", "20", "-warmY",
                           "-cpsep", "10:10:100:20", "-cpiters", "120", "-cpquiet", "12", "-c", str(optimum + 1000),
                           "-coarse", str(pool)]
                if fixes:
                    options += ["-x", ",".join(f"{i}:{j}" for i, j in fixes)]
                call = subprocess.run([args.solver, *options, str(instance)], capture_output=True, text=True, env=env)
                assert call.returncode == 0, (call.returncode, call.stdout, call.stderr)
                output = call.stdout
                activity = re.findall(r"coarse \d+ active (\d+)", output)
                assert activity and max(map(int, activity)) > 0, output
                pattern = (r"^(?:cp (?:it\s+\d+ cert|solve\s+lb|quiet round\s+lb)"
                           r"|lazycert exit (?:last check|stored best): cert) (\S+)")
                finite = [v for s in re.findall(pattern, output, re.M) if math.isfinite(v := float(s))]
                assert "cp quiet round" in output and len(finite) >= 3 and max(finite) <= optimum + 1e-8, output
                worst = max(worst, max(finite) - optimum)
                lazy += output.count(", lazy)")
                checked += len(finite)
                cases += 1
    assert lazy > 0
    print(f"{cases} active cases, {checked} finite certificates, {lazy} lazy checks, maximum certificate minus exact optimum {worst:.9g}")


if __name__ == "__main__":
    main()
