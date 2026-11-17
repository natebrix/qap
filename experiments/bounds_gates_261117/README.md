# Cut and ADMM experiment gates — 17 November 2026

All three features were implemented and tested. None demonstrated a material end-to-end improvement sufficient to recommend a production change or a better whole-instance solve estimate. All experiment jobs have finished. Production was not changed.

Baseline: cognac `ba47c211` (solver source and recipes identical to `7fdac382`). The study executes the coarse-label cut, CP congruence scaling and regional second-moment equality proposals from the November 17 review.

| Track | Most informative result | Decision |
|---|---|---|
| Coarse-label cuts | Selected tho40 pools lose at equal CPU; 600-step tai60b comparison gains 105946 at equal iterations but loses 51005 at equal CPU | No production promotion |
| CP congruence scaling | Same crossing check as stock for endpoint targets in all three screens | Leave off; no material benefit measured |
| Regional equalities | One node fathoms in 663 vs 687 CP iterations; complete subtree follows same 38 nodes with 1962 vs 1811 total DNN iterations | Leave off; phase improvement did not reduce combined work |

See [RESULTS_LEDGER.md](RESULTS_LEDGER.md) for exact figures, costs and interpretation. [verification/README.md](verification/README.md) documents the independent exact checks. [REPLAY.md](REPLAY.md) gives common warm-state generation and replay commands. Detailed implementation reports and raw measurements are on the individual experiment branches.

## Experimental branches

- `research/coarse-label-cuts-261117`: final head `26260591`; measured core solver implementation `6e53f4ff`.
- `research/cp-congruence-261117`: final head `79c5dc37`; primary implementation `f6b2bbe0`, safe homogenizing-weight fallback `d175ac1f`.
- `research/cp-region-equalities-261117`: final head `57e6a285`; primary implementation `0fc7508b`, imported-region capacity fix `f734af31`.

Each idea is off by default or restricted to the experimental driver. They were tested independently. Combining congruence scaling with the equality block requires a weighted equality Gram derivation; overlapping patches must not be combined mechanically.

## Protocol and controls

All experiments used one BLAS thread. Algebraic and feature-on certificate checks preceded performance interpretation. Native warm comparisons received identical serialized parent state and cut pool. A separate forced-line seed build generated those inputs; evaluated candidates did not contain the seed patch. The state history begins at a controlled cold anchor, not the live production run's exact root state.

Source was committed before measurements. The coarse experiments used shared same-node checkpoints, production numerical settings, explicit iteration budgets and matching regional refresh schedules. Those driver experiments omit native branch-and-bound stopping and inheritance. Native warm replays and the completed 38-node subtree used the committed production recipe. Separation and candidate setup were charged in equal-CPU cut comparisons.

Development and screening jobs overlapped. The 600-step tai arms ran concurrently. Timings are screening observations, not idle-machine speedup estimates; small differences were not promoted. The completed subtree's identical ordered nodes and increased total iteration count are independent of this timing caveat.

## Baseline and correctness

The baseline full `make test` passed, including 4,536 brute-force cases, resume, parallel, coordinator, settings, exact regression and estimator checks. The optional AWS mock scenario skipped because boto3/moto were unavailable. Candidate branches passed full default tests, additional brute-force seeds and feature-specific exact checks; see their reports and the independent verification directory.

Compiler overrides `DNN_LAPACK` and `DNN_LAPACK_DEFS` must not reach solver runtime: strict settings rejects unknown `DNN_*` names. For `make test`, use `PYTHON='env -u DNN_LAPACK -u DNN_LAPACK_DEFS python3'` while retaining those make overrides for linking. An initial baseline invocation failed solely because make exported the compiler variables; the clean-runtime test passed.

## Publication

All work is committed locally. Automatic approval review rejected the initial protocol-branch push because it interpreted the user's authorization as local experiments, not GitHub publication. No remote branch was created. Explicit permission is required before publishing the four research branches.
