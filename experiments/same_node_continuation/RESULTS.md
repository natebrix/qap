# Same-node continuation: findings and decision

Based on cognac `7392c10e673707061f98b44c30eea1bd1c5109e2`. Experiment date: 2 November 2026, Los Angeles.

**Result:** the formal full-spectrum continuation test passes. Restoring Y as well as Z also improves the certified bound in all three short, equal-iteration large-instance probes. This supports a full cut-phase A/B; it does not establish a search speedup or justify a 4–12x forecast.

## What the formal test establishes

With a fixed objective, node, indexing, parameters and full spectral projection, the unaccelerated ADMM/PRSM recurrence is a deterministic map of (Y,Z). Loading both matrices reproduces the interrupted recurrence; loading only Z does not preserve its input state. See the [contract and reproduction instructions](README.md).

The implemented regression compares 45 uninterrupted iterations with 10+15+20, using aligned certification checkpoints. Across 72 combinations of instance, full-spectrum precision and ADMM/PRSM:

- Maximum final Y/Z entry difference: **0**.
- Maximum difference between the uninterrupted best certificate and the maximum of the chunk certificates: **0**.
- Z-only negative controls differ in **72/72** cases; largest state difference **1.22938177428831**.

This is deterministic regression evidence under the documented assumptions, not a claim of bitwise reproducibility across every compiler and LAPACK implementation. In particular, certificate Cholesky hints reset between calls.

Additional validation passed:

| Check | Result |
|---|---:|
| Certificates checked against exact enumerated optima | 648 |
| Changed-objective multiplier-rescaling comparisons | 72 |
| Actual changing-cut-pool cases, including quiet rounds | 108 |
| Malformed/misapplied starts rejected | 7 |
| Existing brute-force solver regression | 4,536/4,536 |
| Existing checkpoint/resume regression | 5/5 |
| Existing parallel and seven-instance benchmark checks | Passed |
| Existing DNN quick validation | Passed |

The 72 rescaling comparisons are included in the 648 certificate checks. The 108 pool checks are additional. The exact-optimum tests use 18 seeded problems of order 3–8, signed matrix entries, nonzero linear costs and shifts. No certificate formula, rounding margin or child-transfer policy changed.

## Large-instance paired probes

Each comparison starts from the same 200-iteration state and spends another 100 iterations. The reported split-arm bound is the maximum certificate across the first and second calls. All three use tracked eigenspaces, PRSM, one BLAS thread, exact inner LAPs, degree 3 and tight-lambda 2. Early stops are disabled.

| Node / fixed objective | Z-only restart | Y/Z restart | Improvement |
|---|---:|---:|---:|
| wil50 node 13, uncut | 48,526.171632 | 48,572.799702 | **+46.628070** |
| wil50 node 1, archived 1,000 triangles at multiplier 1 | 48,372.824570 | 48,404.428553 | **+31.603982** |
| sko49 depth-4 node 13, uncut | 23,211.860279 | 23,229.222595 | **+17.362316** |

These are improvements in the returned certified lower bound, not primal-objective improvements or projection distances. The fixed-cut case holds its L_add and constant shift fixed throughout all arms; it does not test the full current adaptive cut pool. The two uncut cases are also early-iteration probes, not converged quiet-round comparisons.

Explicit zero-based assignments identify the nodes independently of historical numbering:

- wil50 node 13: `21:28,8:41,28:19`.
- wil50 node 1: `34:40,28:22,14:29`.
- sko49 depth-4 node 13: `7:48,26:41,2:12,44:30`.

A separate cold 300-iteration reference gives:

| Node | Uninterrupted bound | Y/Z split minus uninterrupted | Max Y difference | Max Z difference |
|---|---:|---:|---:|---:|
| wil50 node 13 | 48,572.727578 | +0.072124 | 0.000152003 | 0.00182034 |
| wil50 node 1, fixed cuts | 48,404.322986 | +0.105567 | 0.000143157 | 0.0014442 |
| sko49 depth-4 node 13 | 23,229.094543 | +0.128052 | 0.00011812 | 0.0025998 |

The tracked eigensolver retains more than Y and Z: its basis, previous projected matrix, spectral estimates, refresh counters and expansion seed affect subsequent iterations. This patch resets that workspace. Consequently Y/Z restoration is a warm initialization in tracked mode, not exact continuation of the whole implementation. The uninterrupted arms diagnose this distinction; their 300-iteration states need not match the split arms.

CPU values are retained in the logs, but probes and validation ran concurrently. No timing speedup is claimed. Three hand-selected nodes also cannot establish a population average, tail behavior or tree-size reduction.

## Recommendation for Fable

Use the new `-warmY` option for a paired A/B of the full current cut protocol on the same representative wil50 and sko49 nodes, including near-converged quiet-round starts. Keep the node, cutoff, pool settings and first-iteration budget identical. Record the best certified bound, total phase CPU and additional fathoms. Because Y changes cut separation and multiplier updates, this full-phase comparison measures the complete resulting policy; the frozen-objective probes above isolate the restart choice more directly.

After that, compare complete subtrees and include cut cost before claiming a core-second reduction. The present +17 to +47 results are a reason to run that experiment, not quantities that can be added to existing late-stage cut gains.

A separate follow-up can persist the tracked spectral state and rerun the split-versus-uninterrupted identity test. This should preserve the basis and associated refresh/spectral state together. Restoring a basis alone is not a complete checkpoint. Acceleration history, certification schedule, stopping counters and the running best need an explicit policy as well if the aim is a fully resumable solver session.

The experiment remains **off by default**. `make continuationcheck` runs the formal suite; the README lists the large probes. This research branch has not been merged into cognac and has not completed the repository's production performance-benchmark gate.

## Reproduction record

- [Formal test output](results/formal_checks.log)
- [wil50 node 13](results/wil50_n13_subspace.log)
- [wil50 node 1, fixed cuts](results/wil50_n1_fixed_cuts_subspace.log)
- [sko49 node](results/sko49_d4_13_subspace.log)
- [Existing regression output](results/regression.log)
- [DNN quick validation](results/dnn_quick_validation.log)
- [CLI wiring smoke test](results/warmY_cli_smoke.log), using fixed unit multipliers; this is not a performance comparison.
- [Environment and build notes](results/manifest.json)

The dense fixed-cut objective is regenerated from an already versioned specification, with its expected SHA-256 in the README. The final malformed-start overflow guard was added after the large probes were launched; it changes no arithmetic for their valid finite starts. The final formal suite and CLI smoke test were run on the final implementation.
