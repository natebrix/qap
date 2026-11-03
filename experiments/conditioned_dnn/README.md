# Conditioned DNN with reuse of child computations

Isolated experiment off cognac `6293cb7a26544ee62bbfde46186ac1b8629b943f`,
3 November 2026. No production bound, certificate, branching rule, or default
parameter changes. The driver calls the existing `dnn_bound`, `child_lb`,
`dnn_warm_child`, and opt-in same-node `warm_primal` APIs.

## Predeclared experiment

Does a **complete** facility-to-location disjunction strengthen a parent
bound at an acceptable cost, and can the resulting child computations be
used when the parent branches?

For a parent P and a selected free facility i, enumerate **every** remaining
location j, without symmetry reduction or sampling. Each child receives the
parent's mapped Z (the current production initialization). Let bP be the
certified parent bound and c_j its cheap certified child bound. After a
short child solve with bound p_j, store

    b_j = max(bP, c_j, p_j),  b_conditioned = min_j b_j.

Store each child's final Y/Z and its own running-best bound on disk. Bind
the snapshot to its exact reduced objective, original assignment labels,
and fixed solver configuration. On reuse, load that same child's Y/Z with
`warm_primal=1` and keep `max(b_j, new_certificate)`. Never use the primal
objective as a bound. Stop a child only after its certificate exceeds the
specified cutoff; a non-fathomed child receives its full iteration budget.

PRSM settings: gamma1=gamma=1, beta=.3, exact inner LAPs, check every 10,
subspace degree 3, tight-lambda 2, refresh 40; no convergence, stall, or
hopeless stopping. One BLAS thread per process. The tracked eigensolver's
workspace is **not** serialized, so reuse is a Y/Z warm restart, not an
identical continuation of its full numerical state.

Cases selected before measuring:

| Case | Parent fixings (zero based) | Facility | Parent iterations | Probe | Additional reuse | Fresh direct child |
|---|---|---:|---:|---:|---:|---:|
| nug20 validation/pilot | root | 0 | 100 | 50 | 100 | 150 |
| wil50 existing near-miss node 9 | 26:27,8:38,34:28 | 37 | 200 | 50 | 100 | 150 |

Facility 37 is the already recorded branch facility for node 9 in
`wil50_feasibility/tuning/giveup/all32.txt`. This experiment includes all
47 children, not just the eight depth-4 examples in that file. Cutoffs are
2568.05 for nug20 and 48814.05 for wil50, matching the repository's
even-integer objective convention and incumbent values 2570/48816.

## Four comparisons and honest accounting

1. **Probe:** parent solve plus all short child solves. Report the minimum
   child bound, weakest location, children fathomed, and total CPU.
2. **Reuse:** probe, then restore surviving child Y/Z for the extra budget.
   Charge the complete parent and probe cost. Already fathomed children
   require no extra solve.
3. **Discard:** same probe, then recompute surviving children from the
   parent's Z for probe+extra iterations. This measures the waste of throwing
   away work done by lookahead. Retain previous certified bounds.
4. **Direct:** branch immediately after the parent, with probe+extra
   iterations per child, without a preliminary probing pass. This is the
   essential control: reuse is not automatically faster than ordinary B&B.

Also run an **unsplit-parent CPU control** from cold, giving it the total
solver CPU of parent+complete probe. Its CPU stopping check is every ten
iterations, so report actual checkpoint overshoot rather than claiming
exact equality. It tests the first conditioning pass, not the longer
probe+continuation policy. The control has a 100,000-iteration cap, checked
in the report if it binds.

Raw solver CPU is the primary cost measurement. Also retain process wall
times, snapshot write times, and bytes. Probe workers run concurrently and
the second phase overlaps the parent control: summed process wall times
are not a sequential speed benchmark. Checkpoint read/write and hashing
are real overhead, outside solver CPU. No whole-search speedup is claimed.

Promotion requires a useful bound/fathoming improvement over the CPU
control, then comparable or better child decisions and total cost than
**direct** branching. Savings against **discard** alone only demonstrate
the benefit of avoiding known duplicate work.

## Validation and running

Build with the repository's supported BLAS/LAPACK settings:

```sh
make CC=gcc
make -C latest dnnconditioned dnncontinuation
./latest/dnnconditioned selftest /tmp/conditioned-selftest.bin
./latest/dnncontinuation
make test CC=gcc
```

`selftest` enumerates all permutations of ten signed, asymmetric QAPs of
orders 3–7, with nonzero linear costs and constant shifts. It checks every
child optimum, complete disjunction coverage, inherited cheap certificates,
probe/reuse/fresh certificates, full-spectrum split-versus-uninterrupted
Y/Z and bound equality, exact snapshot roundtrips, and wrong-node cache
rejection. The production certificate and its numerical margins are used
unchanged. The driver explicitly preserves source shifts because
`reduceQAP` is designed to reduce an original zero-shift instance.

Commit the experiment before measuring, then:

```sh
python3 experiments/conditioned_dnn/run.py --instance qaplib/nug20.dat \
  --facility 0 --cutoff 2568.05 --parent 100 \
  --out experiments/conditioned_dnn/results/nug20
python3 experiments/conditioned_dnn/run.py --instance qaplib/wil50.dat \
  --fixes 26:27,8:38,34:28 --facility 37 --cutoff 48814.05 \
  --out experiments/conditioned_dnn/results/wil50_node9
```

The runner refuses to overwrite an attempt. Every invocation produces a
raw log and JSON result; manifests record commit/binary/input hashes and
cache hashes. A summary is produced only after complete child coverage and
both controls finish. Large regenerable Y/Z caches are ignored by git;
their hashes, byte counts, and all scalar results are versioned. Snapshots
are a local native-binary experiment format, not portable B&B checkpoints.

Results will be recorded separately in `RESULTS.md`.
