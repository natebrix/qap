# Conditioned DNN: measured results

Experiment off cognac `6293cb7`, 3 November 2026. All bounds below are certified lower bounds. Child-policy CPU totals include the initial parent solve; the unsplit control is charged for its full cold solve. No production search code changed.

## Finding

The fixed-budget conditioning policy is **not supported for promotion on this wil50 pilot**. Its short-probe pass gained zero over the initial parent. After the additional child work, reuse gives 48767.419682 at 3963.447 solver CPU seconds, compared with 48770.604275 at 1832.520 seconds for the unsplit control.

The cache successfully reuses child computations. Its savings against discarding probes should be distinguished from any advantage over direct branching: without early stopping, the per-child iteration budgets are 150 for reuse, 200 for discard, and 150 for direct branching. The direct comparison below is the relevant one for judging extra search efficiency.

This is a **finite-budget policy result for one wil50 node and one facility**. Exact conditioned DNN dominates the exact parent relaxation; finite-run certified bounds can be below those optima. The candidate retains the 200-iteration parent bound, while the separate control spends more work on the parent itself. The experiment does not rule out a stronger parent warm start, adaptive refinement of the weakest children, another facility, or conditioning with inherited cuts.

## Completed cases

### nug20 root, facility 0

Complete split: **20 children**. Cutoff: 2568.05. Parent/probe/additional iteration caps: 100/50/100.

| Policy | Certified parent bound | Total solver CPU (s) | Children fathomed | Child iterations |
|---|---:|---:|---:|---:|
| Initial parent | 2494.698459 | 1.435 | — | — |
| Unsplit parent, probe-matched CPU | 2506.328407 | 10.966 | — | — |
| Conditioned: short probes | 2528.987245 | 10.937 | 0/20 | 1000 |
| Probes + cached Y/Z continuation | 2540.955010 | 24.637 | 4/20 | 2680 |
| Probes + recomputation from parent Z | 2540.967487 | 34.022 | 4/20 | 3680 |
| Direct branching, no probes | 2540.967487 | 24.520 | 4/20 | 2680 |

The short conditioned bound is **+22.658838** relative to the CPU-matched unsplit parent. The latter ran 1170 iterations (time_limit); its stopping checkpoint exceeded the target by 0.030 CPU seconds.

Reuse saved **27.58%** of total solver CPU versus discarding the probes. Its CPU ratio to direct branching was **1.0048**. The final minimum bound differed from direct branching by **-0.012477**. Maximum absolute per-child bound difference: 0.029975. Per-child differences include different early-stop iterations and overshoot after the cutoff; they are not certificate errors.

Locations with different reuse/direct pruning decisions: **none**.

Weakest locations: short probes [1]; reuse [17]; direct [2]. Retained child Y/Z payload: 0.0419 GB (0.0391 GiB). Summed probe snapshot write time: 0.095 s.

[Summary JSON](results/nug20/summary.json), [per-child CSV](results/nug20/children.csv), [manifest](results/nug20/manifest.json).

### wil50 node 9, facility 37

Complete split: **47 children**. Cutoff: 48814.05. Parent/probe/additional iteration caps: 200/50/100.

| Policy | Certified parent bound | Total solver CPU (s) | Children fathomed | Child iterations |
|---|---:|---:|---:|---:|
| Initial parent | 48653.907938 | 139.490 | — | — |
| Unsplit parent, probe-matched CPU | 48770.604275 | 1832.520 | — | — |
| Conditioned: short probes | 48653.907938 | 1827.518 | 5/47 | 2330 |
| Probes + cached Y/Z continuation | 48767.419682 | 3963.447 | 32/47 | 5110 |
| Probes + recomputation from parent Z | 48767.221341 | 5393.183 | 32/47 | 7190 |
| Direct branching, no probes | 48767.221341 | 3896.681 | 32/47 | 5090 |

The short conditioned bound is **-116.696337** relative to the CPU-matched unsplit parent. The latter ran 2100 iterations (time_limit); its stopping checkpoint exceeded the target by 5.001 CPU seconds.

Reuse saved **26.51%** of total solver CPU versus discarding the probes. Its CPU ratio to direct branching was **1.0171**. The final minimum bound differed from direct branching by **+0.198341**. Maximum absolute per-child bound difference: 21.998147. Per-child differences include different early-stop iterations and overshoot after the cutoff; they are not certificate errors.

Locations with different reuse/direct pruning decisions: **none**.

Weakest locations: short probes [3]; reuse [35]; direct [35]. Retained child Y/Z payload: 3.3702 GB (3.1388 GiB). Summed probe snapshot write time: 3.395 s.

[Summary JSON](results/wil50_node9/summary.json), [per-child CSV](results/wil50_node9/children.csv), [manifest](results/wil50_node9/manifest.json).

## Relation to the existing cut policy

The archived [shorter cut-phase study](../same_node_continuation/budget_study/RESULTS.md) already reports fathoming this same wil50 node 9 with Y12s at 720 total iterations and retained bound 48816.286927. Its 1050.166 whole-process CPU seconds come from a different run/build environment and are not a matched timing control here. Consequently a gain over the uncut parent is not evidence of a gain over cognac's cut-strengthened policy. This experiment should not displace that policy on the strength of an uncut comparison.

## Interpretation boundaries

- The complete disjunction is essential: these minima include every free location, without symmetry reduction or cherry-picking children.
- Reuse versus discard measures avoidable repeated work. Direct branching is the relevant baseline for a search-speed claim. No complete subtree search was run.
- This compares uncut DNN. It does not establish superiority over the existing 12-round cut-strengthened DNN policy or measure inherited cut pools.
- The saved state includes Y/Z and the best certified bound, not the tracked eigensolver workspace. Thus fixed-total-iteration results can differ from direct solves; no exact continuation claim is made for tracked mode.
- Three one-thread child workers ran concurrently; the parent control overlapped the second phase. CPU ratios are diagnostic on this machine, not controlled sequential wall-clock speedups. Disk I/O, cache hashing, scheduling, and file reads are excluded from solver CPU; raw process wall times are retained.
- Child states are kept on disk, one loaded per active worker. An in-memory implementation retaining every child would need the reported aggregate payload. A production cache would need an explicit memory/eviction policy.

## Validation and provenance

- Ten fully enumerated signed/asymmetric QAPs of orders 3–7: 50 children, all inherited/probe/reuse/fresh bounds safe; disjunction coverage exact; 50 wrong-node snapshot rejections. Full-spectrum state and bound continuation errors were both zero.
- Existing continuation checks: 72 identity cases, 648 safe certificates, 72 rescaling cases, 108 changing-cut cases, seven malformed states; passed.
- `make test`: 4,536 brute-force runs, five resume cases, parallel tests, and seven exact regression matches; passed. Logs are in `results/`.
- The post-run [audit](results/validation.json) checks all raw child logs, commands/budgets, complete coverage, cache-bound inheritance, input/binary hashes, summary hashes, and every locally retained child-state hash.
- Measurements used local commit `0ae156a31ae9c2f0af53abe4e2f59c83aef34dd2`. GitHub publication through the connector produced commit `b36b6b6cbe71e56f919189e8f4a25b58b00fcf9c` with the **identical tree** `994a1903d6f002be32c0923c3e3cf9f693077085`. Reporting files were added afterward; the measured executable was not changed.
- Intel Xeon Platinum 8573C; one OpenBLAS thread per solver; SciPy LP64 OpenBLAS, compiler `-O3 -march=native -ffp-contract=off`. AMX was not enabled for these runs. Full command lines and executable/input SHA-256 hashes are retained in the manifests.
