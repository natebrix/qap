# Prelaunch validation

- `make-test.log`: full initial `make test` passed (4,536 brute-force runs, resume and parallel checks, seven-instance exact regression).
- `make-test-final.log`: final brute-force, resume and parallel checks passed; the regression freshness guard then caught a header documentation edit newer than the binary.
- `regression-final.log`: rebuilt binary, final exact regression passed. No numerical regression failure occurred.
- `pool-regression.log`: both-axis round-trip objective and DNN value equality, 72 child remaps / 8,640 permutation lifts, two one-solve frozen integrations, and 20 invalid-input rejections.
- `continuation.log`: existing continuation, certificate safety and rescaling regression passed.
- `smoke-{old,new,save,load}.log`: prior/default/save driver sequences match; loading at a child performs exactly one quiet solve. `smoke-parent.pool` is the saved fixture.
- `analysis.log`: queue, parser and gate checks passed.
