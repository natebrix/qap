# Full cut-phase continuation A/B

Protocol committed before launch at `fda28a6922b2c7661142264a93676c2b5948e31c` on `research/same-node-continuation`. Solver code remains based on cognac `7392c10`; cognac was checked again at `d88b3e0`, whose intervening changes are a search experiment script, parameters and partial log, with no DNN implementation change.

## Predeclared comparison

Two nodes, each run with the existing Z-only initialization and the opt-in Y/Z initialization. Within each pair the command differs only by `-warmY`. The cut pool starts empty. Four processes run concurrently, each pinned to a distinct physical core with one BLAS thread.

| Node | Zero-based assignments | Grid | Incumbent / cutoff |
|---|---|---|---|
| wil50 node 13 | `21:28,8:41,28:19` | 5 x 10, side 0 | 48,816 / 48,814.05 |
| sko49 depth-4 node 13 | `7:48,26:41,2:12,44:30` | 7 x 7, side 0 | 23,386 / 23,384.05 |

Protocol: 400 initial uncut iterations, 12 cut rounds of 300, up to 1,000 new regional triangle cuts per round, uncapped Polyak target best + 50, theta 1, and up to 1,500 quiet iterations from the best round's multipliers. Cut-round stalling is disabled; valid early fathoming remains enabled. Maximum total is 5,500 DNN iterations per arm.

Use the existing dnncuts eigensolver settings (tracked mode, degree 4, PRSM, default tight-lambda 0 and inner-LAP auto policy). These differ from the earlier short diagnostic probes (degree 3, tight-lambda 2, exact inner LAP, checkpoints every five); the full-phase comparison matches the established experiment protocol instead.

Primary outcomes: final best certified lower bound and fathoming. Secondary outcomes: total iterations and whole-process CPU, including cut separation and objective assembly. Timing is conditional on this four-process workload: shared-memory contention, core differences, and different finishing times limit interpretation. No serial production benchmark or whole-search speedup is inferred.

The uncut first bounds must match within each pair. Cut selection and multiplier trajectories after the first restart may differ; that is part of the policy being tested. Both arms retain the strongest certificate, and the quiet phase uses the corresponding best-round state.

## Running and reproducing

Build `dnncuts` with the same LAPACK settings as the formal checks, then:

```bash
python3 experiments/same_node_continuation/run_cut_phase_ab.py --out /tmp/qap-cut-phase-ab
python3 experiments/same_node_continuation/summarize_cut_phase_ab.py /tmp/qap-cut-phase-ab
```

The runner refuses an existing output directory; use `--out` for a fresh repetition. The wall deadline is 3.75 hours, including all four arms, to stay within the user's four-hour authorization. Every arm retains its log even if interrupted. The manifest records the exact source commit, binary SHA-256, commands, core assignments, elapsed wall time, whole-process CPU and maximum RSS.

## Results — completed

**The full A/B completed in 7,109.526 wall seconds (1 hour 58 minutes 30 seconds).** All four arms finished successfully and used 5,500 iterations. Both pairs had identical initial bounds at the logged precision. The only solver-option difference within a pair was `-warmY`.

Y/Z restoration improved the final certified bound by **7.939138 on wil50** and **0.864434 on sko49**. **Neither node was fathomed by either arm.** These small final gains supersede the larger short-probe gains when judging the current full cut-phase policy.

| Node | Z-only final bound | Y/Z final bound | Improvement | Fathomed | Remaining gap to cutoff Z / YZ |
|---|---:|---:|---:|---|---:|
| wil50 node 13 | 48,754.619389 | 48,762.558527 | **+7.939138** | Neither | 59.430611 / 51.491473 |
| sko49 depth-4 node 13 | 23,379.294892 | 23,380.159326 | **+0.864434** | Neither | 4.755108 / 3.890674 |

The quiet solve recovered most of the baseline's deficit:

| Node | Best Z-only before quiet | Best Y/Z before quiet | Pre-quiet lead | Final lead |
|---|---:|---:|---:|---:|
| wil50 node 13 | 48,708.8084 | 48,757.0704 | +48.2620 | +7.939138 |
| sko49 depth-4 node 13 | 23,365.3361 | 23,377.8950 | +12.5589 | +0.864434 |

Per-round values have four printed decimal places; final values have six. The lead shrank from +48.2620 to +7.939138 on wil50 and from +12.5589 to +0.864434 on sko49. This is consistent with much of the early advantage being faster inner convergence. The cut pools and multiplier paths also differ, so this is not an exact decomposition of the effect into convergence and pool quality.

Whole-process CPU, including separation and assembly:

| Node | Z-only CPU seconds | Y/Z CPU seconds | YZ/Z ratio | Wall minutes Z / YZ |
|---|---:|---:|---:|---:|
| wil50 node 13 | 7,024.041 | 6,875.861 | 0.9789 | 118.49 / 115.96 |
| sko49 depth-4 node 13 | 5,893.691 | 5,868.306 | 0.9957 | 99.46 / 98.96 |

Observed CPU reductions are 2.11% on wil50 and 0.43% on sko49. This was one concurrent run per arm, with different cores and changing shared-memory contention as arms finished. Those small differences do not establish a repeatable speedup. Total CPU across the four arms was 25,661.899 seconds (7.13 core-hours).

## Assessment and next decision

The experiment supports a modest finite-budget improvement on these two nodes. It does not demonstrate an extra fathom, a whole-search reduction, or a contribution approaching the required 4–12x. Keep the option experimental pending a broader, representative comparison.

There is a useful budget-tuning observation: on wil50, Y/Z already certified **48,757.0704 after 4,000 iterations**, before the quiet solve. That exceeds the Z-only arm's eventual 5,500-iteration result of **48,754.619389**. Thus a shorter Y/Z phase could match the current full baseline quality on this particular node. This is a retrospective prefix comparison, not a validated stopping rule or an exact time-to-target measurement; the logs do not record every quiet checkpoint. Sko49 did not match its baseline's final bound by the end of the cut rounds.

The next focused test would use a predeclared shorter Y/Z budget against the current baseline on held-out near-miss and clear-miss nodes, retaining the bound, fathom and total-cost gates. A subtree A/B would then price any benefit in search work. No further solver experiments were launched after the four declared arms completed.

This standalone protocol always runs the quiet solve when no fathom has occurred. Production search can additionally gate the quiet phase on its primal/cutoff policy and can enter from different warm states; these results do not directly measure that complete search policy.

## Retained evidence

- [Summary JSON](results/full_cut_phase_ab/summary.json)
- [Manifest, exact commands, binary hash and process statistics](results/full_cut_phase_ab/manifest.json)
- [wil50 Z-only log](results/full_cut_phase_ab/wil50_n13_Z_only.log) and [Y/Z log](results/full_cut_phase_ab/wil50_n13_YZ.log)
- [sko49 Z-only log](results/full_cut_phase_ab/sko49_d4_13_Z_only.log) and [Y/Z log](results/full_cut_phase_ab/sko49_d4_13_YZ.log)
- [Runner](run_cut_phase_ab.py) and [comparison validator](summarize_cut_phase_ab.py)

The solver and runner used for the measurements are pinned at `fda28a6`. After completion, the runner received a macOS-compatible path that leaves CPU placement to the OS when affinity APIs are absent, plus a platform-correct RSS conversion. This compatibility path was checked without rerunning the solver and has not had a full macOS A/B here. The recorded experiment used the original Linux runner with explicit core pinning. Reproduction commands above use a fresh output directory so the committed logs are never overwritten.
