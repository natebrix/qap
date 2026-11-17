# Cut and ADMM experiment gates — 17 November 2026

Status: experiments in progress; no performance result yet.

Baseline: cognac `ba47c211` (solver source identical to `7fdac382`). This study executes the coarse-label cut, CP congruence scaling, and regional second-moment equality proposals from the 17 November review. Production is not modified.

## Experimental branches

- `research/coarse-label-cuts-261117`: exact local multi-label inequalities in CP.
- `research/cp-congruence-261117`: CP-compatible tensor congruence scaling.
- `research/cp-region-equalities-261117`: affine regional consistency projection.

Each idea is opt-in. A separate forced-line seed build generates replay inputs only; evaluated candidates do not contain the forced-line patch. All compared variants must receive identical serialized node input and parent warm state/pool.

## Protocol

1. Verify explicit projection algebra, cut validity on permutations, and certificate accounting; run existing correctness checks with switches off and feature-specific checks with switches on.
2. Commit implementation before final measurement. Use single-thread SciPy OpenBLAS consistently.
3. Screen matched certificate progress on tho40t hard depth-3/depth-5 nodes and tai60b's critical depth-2 branch; add near-tight nodes as controls.
4. Preserve production cut regions, penalty, continuation policy, cutoff, and eigen settings. Diagnostic-driver runs with simplified stopping are labeled separately.
5. Run final baseline/candidate timings sequentially. Intermediate concurrent development timing is not a speedup claim.
6. Advance only promising variants to complete subtree comparisons. Report all tested variants, negative results, overhead and correctness limitations.

## Environment note

The installed SciPy OpenBLAS exports `scipy_`-prefixed BLAS symbols. Build overrides `DNN_LAPACK` and `DNN_LAPACK_DEFS` must not reach solver runtime: strict settings correctly rejects unregistered `DNN_*` names. `make test` can use `PYTHON='env -u DNN_LAPACK -u DNN_LAPACK_DEFS python3'` while retaining these make overrides for linking. An initial baseline test failed solely because those compiler variables were exported by make; after removing them from runtime, all 4,536 brute-force cases and all five resume cases passed. The full baseline `make test` passed, including coordinator, settings, regression and estimator checks. The optional AWS mock scenario was skipped because boto3/moto were unavailable.
