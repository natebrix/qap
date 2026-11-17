# Working status

This file is a checkpoint, not the final result report.

- Baseline solver source: `7fdac382`, identical to solver source in base `ba47c211`.
- Baseline full `make test` passed. Optional AWS mock skipped: dependency unavailable.
- Scaling implementation: `f6b2bbe0`, branch `research/cp-congruence-261117`; full make test passed, dense reference and enabled exact comparisons pass.
- Equality implementation: `0fc7508b`, branch `research/cp-region-equalities-261117`; full make test and brute-force seeds 2/3 passed, enabled exact comparisons pass.
- Coarse-label implementation: branch `research/coarse-label-cuts-261117`, implementation/testing in progress. Sparse generic cut half-spaces, exact dyadic RHS validation, same-node replay. CP+Ladd explicitly rejected; generic cuts use the validated hook instead.
- Independent review: three variants passed direct feature-on exact comparisons including executed lazy-exit recovery; verification scripts/results being packaged.
- Preliminary equality cold tho40 depth-5 diagnostic: +27.17 bound units at 400 CP iterations, no fathom. Development timings are not speedup evidence.
- Warm tai35b scaling screen running; common warm tho40/tai60b production-recipe corpus generation in progress.

The seed branch's forced-line patch is exclusively for creating shared inputs. It must not be in benchmark binaries. Generated states originate from controlled recipe runs, not exports from the live cloud production run.

Local workspace staging uses sibling directories `qap`, `qap-coarse`, `qap-scaling`, `qap-equalities`, `qap-experiment-seed` and `experiment_run_261117`. The replay manifest/commands make the states reproducible.

Remote publication is pending explicit approval: an automatic approval review rejected the initial protocol-branch push because experiment authorization was not interpreted as permission to publish a new GitHub branch. No remote branch was created. Local experiments and commits continue.
