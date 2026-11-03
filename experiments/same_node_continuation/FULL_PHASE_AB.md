# Full cut-phase continuation A/B

Protocol committed before launch at `fda28a6922b2c7661142264a93676c2b5948e31c` on `research/same-node-continuation`. Solver code remains based on cognac `7392c10`; cognac was checked again at `d88b3e0`, whose intervening changes are a search experiment script, parameters and partial log, with no DNN implementation change.

## Predeclared comparison

Two nodes, each run with the existing Z-only initialization and the opt-in Y/Z initialization. Within each pair the command differs only by `-warmY`. The cut pool starts empty. Four processes run concurrently, each pinned to a distinct physical core with one BLAS thread.

| Node | Zero-based assignments | Grid | Incumbent / cutoff |
|---|---|---|---|
| wil50 node 13 | `21:28,8:41,28:19` | 5 x 10, side 0 | 48,816 / 48,814.05 |
| sko49 depth-4 node 13 | `7:48,26:41,2:12,44:30` | 7 x 7, side 0 | 23,386 / 23,384.05 |

Protocol: 400 initial uncut iterations, 12 cut rounds of 300, 1,000 new regional triangle cuts per round, uncapped Polyak target best + 50, theta 1, and up to 1,500 quiet iterations from the best round's multipliers. Cut-round stalling is disabled; valid early fathoming remains enabled. Maximum total is 5,500 DNN iterations per arm.

Use the existing dnncuts eigensolver settings (tracked mode, degree 4, PRSM, default tight-lambda 0 and inner-LAP auto policy). These differ from the earlier short diagnostic probes (degree 3, tight-lambda 2, exact inner LAP, checkpoints every five); the full-phase comparison matches the established experiment protocol instead.

Primary outcomes: final best certified lower bound and fathoming. Secondary outcomes: total iterations and whole-process CPU, including cut separation and objective assembly. Timing is conditional on this four-process workload: shared-memory contention, core differences, and different finishing times limit interpretation. No serial production benchmark or whole-search speedup is inferred.

The uncut first bounds must match within each pair. Cut selection and multiplier trajectories after the first restart may differ; that is part of the policy being tested. Both arms retain the strongest certificate, and the quiet phase uses the corresponding best-round state.

## Running and reproducing

Build `dnncuts` with the same LAPACK settings as the formal checks, then:

```bash
python3 experiments/same_node_continuation/run_cut_phase_ab.py
python3 experiments/same_node_continuation/summarize_cut_phase_ab.py experiments/same_node_continuation/results/full_cut_phase_ab
```

The runner refuses an existing output directory; use `--out` for a fresh repetition. The wall deadline is 3.75 hours, including all four arms, to stay within the user's four-hour authorization. Every arm retains its log even if interrupted. The manifest records the exact source commit, binary SHA-256, commands, core assignments, elapsed wall time, whole-process CPU and maximum RSS.

## Results

Pending completion; no outcome is inferred from partial rounds.
