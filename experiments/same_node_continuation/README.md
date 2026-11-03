# Same-node continuation: formal contract and tests

Research branch based on cognac `7392c10` (2 November 2026, Los Angeles). The experiment is opt-in. The search defaults and child warm-start policy remain unchanged.

## Question and conclusion boundary

Does restoring the primal Y, as well as the dual Z, give a true continuation of the same-node ADMM/PRSM iteration? Does that improve finite-budget progress with the tracked eigensolver?

These are separate questions. The full-spectrum state-identity test has a mathematical justification and a deterministic numerical regression. The tracked solver also has eigenbasis/history state that is not restored by this first implementation. Large-instance probes test progress; they do not establish a whole-search speedup.

## Mathematical contract

For a nontrivial node with m >= 2, fix the node, assignment indexing, lifted objective L and its scale, elimination mask, beta, gamma and gamma1. Disable Halpern, Anderson, and momentum acceleration, CPU/heuristic/convergence stops, and use a full spectral projection. Let P denote the existing box/gangster projection. One iteration is the deterministic map T on the pair (Y,Z):

1. R = projection_PSD(V' (Y+Z/beta) V), and H = V R V'.
2. Z_half = Z + gamma1 beta (Y-H).
3. Y_next = P(H - (L/scale + Z_half)/beta).
4. Z_next = Z_half + gamma beta (Y_next-H).

Consequently T^(a+b)(Y0,Z0) = T^b(T^a(Y0,Z0)). Restoring both matrices implements the same mathematical recurrence. Restoring Z and replacing Y by its cold initialization generally does not.

This statement is about iterates. To compare returned running-best certificates, align certification checkpoints and take the maximum over all chunks; an individual continuation call does not retain the previous call's best certificate. For finite-precision equality, also control BLAS threads, compiler arithmetic, and full-eigensolver mode. The certificate's adaptive Cholesky hint is reset between calls, so bitwise certificate equality on all possible inputs is not claimed; it is checked on the specified fixtures.

For a changed cut objective on the same feasible set, restoring Y is a warm start, not continuation of a single fixed map. Preserve the unscaled multiplier by loading Z_new = Z_old * scale_old / scale_new. The existing certificate applies to the resulting iterate as usual. A valid modified objective plus its correctly adjusted shift still produces a lower bound on the original node. This does not guarantee faster convergence or a better finite-budget bound.

## Minimal opt-in implementation

- `DNNParams::warm_primal`, default 0, loads matching finite Y and Z when explicitly requested. On the iterative path it checks dimensions, scale, finiteness and scaled-Z overflow; malformed requested state returns DNN_ERROR with no finite bound.
- The caller must ensure the same node and indexing. Equal dimensions alone do not identify a node. This flag is not enabled for parent-to-child starts.
- `DNNCutParams::warm_primal`, default 0, applies the option between cut rounds and during the quiet round. Existing `best_mu`/`best_save` selection supplies the matching best-round state.
- `dnncuts -warmY` opts into the experiment. No bnb environment flag or production search default is changed.
- The spectral workspace, acceleration history, stopping counters, and previous running best are not persisted. In tracked mode this is Y/Z warm initialization, not a complete resumable solver session.

## Deterministic test suite

Build with the same BLAS/LAPACK settings used for the solver, then run:

```bash
make continuationcheck
```

`tests/dnn_continuation.C` constructs 18 seeded symmetric QAPs of order 3–8 with signed entries, nonzero linear costs and constant shifts. It enumerates every permutation to obtain exact integer optima.

For both ADMM and PRSM and both double/single full-spectrum projections it compares uninterrupted 45 iterations with chunks 10+15+20, certifying every five iterations. It compares every entry of final Y/Z and the running-best certificate, and runs the legacy Z-only negative control. The declared tolerance is 1e-12 for state, 1e-7 objective units for certificates; raw maximum errors are printed.

The suite separately changes the lifted objective with a valid inequality, forcing a different scale, and compares automatic rescaling with an independently pre-rescaled multiplier. Every checked certificate must stay below the enumerated optimum; changing-cut-pool runs use the optimum as cutoff and must never fathom.

Coverage:

- 72 full-spectrum identity cases and 72 Z-only negative controls;
- 648 checked certificate values, including 72 objective-rescaling pairs;
- 108 actual changing-cut-pool cases, both aggregation axes and all three eigensolver modes, with best-round quiet restoration;
- seven malformed/misapplied-start cases, including scaled-multiplier overflow and a freshly constructed child warm start without Y.

Passing enumeration checks is strong regression evidence, not a substitute for the existing numerical certificate's proof. The new code changes only initialization; no certificate formula or rounding allowance is changed.

## Paired large-instance probes

The diagnostic driver runs four arms in one process: 200 iterations to a common state; 100 more with legacy Z-only restart; 100 more with YZ restart from that same state; and a separate uninterrupted 300-iteration solve. All use one BLAS thread, PRSM, subspace degree 3, tight-lambda 2, exact inner LAPs, and no early stopping. Compare `max(first_lb, continuation_lb)` for each split arm.

```bash
make -C latest dnncontinuation
./latest/dnncontinuation --probe qaplib/wil50.dat '21:28,8:41,28:19' 200 100 2 1
./latest/dnncontinuation --probe qaplib/sko49.dat '7:48,26:41,2:12,44:30' 200 100 2 1
python3 experiments/same_node_continuation/materialize_fixed_cuts.py /tmp/fixed.L.bin
./latest/dnncontinuation --probe qaplib/wil50.dat '34:40,28:22,14:29' 200 100 2 1 /tmp/fixed.L.bin -261
```

The third case holds the archived 1,000-triangle, multiplier-1 objective fixed throughout all arms. The regenerating script uses the already versioned Astra cut specification; the dense 39 MB matrix is not committed. Its SHA-256 is `21dd5cb2947628085ec03ecdad14ed5130c5f2276202621b9e6c11b34d520151`.

These probes ran concurrently with validation work. Iteration budgets are matched; CPU timings are diagnostic and should not be used as controlled speedup estimates. They are short continuation probes, not the full 12-round adaptive cut protocol or a complete subtree A/B.

For an adaptive-pool A/B, run the same existing `dnncuts` command once with and once without `-warmY`; retain the strongest certificate in both arms and charge the entire phase cost. Such an A/B also changes subsequent cut selection because Y changes, so the frozen-objective test above must be interpreted first.

## Results

See [RESULTS.md](RESULTS.md) and retained scalar logs under `results/`.

The subsequent [full cut-phase A/B](FULL_PHASE_AB.md) is complete: +7.939138 on wil50 node 13 and +0.864434 on sko49 depth-4 node 13 after all cut and quiet rounds, with no additional fathoms. Use that follow-up when assessing the full policy; the short probes above measure a different budget. The runner and comparison validator are linked in the follow-up report.

The [shorter-budget follow-up](budget_study/RESULTS.md) is also complete: Y12s (`-warmY`, 12 x 200 plus quiet 1,000) passed all three gates on six wil50 nodes with 64.17% of baseline CPU. The report contains all 30 runs and conditional roadmap cost projections; production defaults remain unchanged.
