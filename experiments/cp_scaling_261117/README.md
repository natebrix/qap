# CP tensor congruence scaling, November 17, 2026

Implementation `f6b2bbe0`, based on `ba47c211`. The measured stock binary
was built from `7fdac382`, which has the same solver and recipes. The opt-in experiment
uses `DNN_CP_DSCALE` for the row-norm exponent and
`DNN_CP_DSCALE_UNTIL` for the window within each CP call. An exponent
of zero selects the existing CP metric. A window of zero uses the metric
through the call. The tested recipes use `w=0.36,wc=2.88`.
An enabled CP experiment also uses the existing `dnn.dscale0` setting
for the homogenizing coordinate. A homogenizing scale whose square
underflows selects the plain metric. The two experiment settings are present
in the state hash, so native replays across builds use
`--allow-settings-change`.

## Projection and certificate

With `Yt = D Y D`, the objective and cut normals transform by inverse
congruence. The CP projector clips at `d_r*d_c`, fixes `Yt[0,0]=d0*d0`,
and retains the fixed zero entries. Region moments remain in the original
units. The diagonal overlap term uses the sum of inverse fourth powers
over the intersection. The pair overlap term uses the sum of products
of inverse squares over distinct unordered pairs in the intersection.
Each term also uses the inverse weight of the corresponding components
on the other side. This keeps the region compression and computes the
correct Frobenius cut norm in the scaled metric.

At a certificate check `Zc` is mapped by `D Zc_t D`. Cut multipliers
remain scalars in the original units. The certificate, including its lazy
snapshots, assembles the original normals and original right-hand sides.
Saved `Y,Z` and shipped pools are in their existing units. When a
window ends the live state and extrapolation checkpoints are converted
and the eigenspace tracker refreshes.

## Reproduction

Build the baseline and experiment with the same BLAS and compiler. The
Linux runs used one BLAS thread and the runtime SciPy OpenBLAS.

```sh
blas=/opt/codex/runtimes/codex-primary-runtime/dependencies/python/lib/python3.12/site-packages/scipy.libs
make -C latest -j2 bnb dnncuts dnncontinuation \
  DNN_LAPACK="-L$blas -l:libscipy_openblas-6cdc3b4a.so -Wl,-rpath,$blas" \
  DNN_LAPACK_DEFS='-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS'
OPENBLAS_NUM_THREADS=1 ./latest/dnncontinuation
OPENBLAS_NUM_THREADS=1 make test \
  PYTHON='env -u DNN_LAPACK -u DNN_LAPACK_DEFS python3' \
  DNN_LAPACK="-L$blas -l:libscipy_openblas-6cdc3b4a.so -Wl,-rpath,$blas" \
  DNN_LAPACK_DEFS='-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS'
OPENBLAS_NUM_THREADS=1 python3 tests/brute_force.py --seed 2 --jobs 2
OPENBLAS_NUM_THREADS=1 python3 tests/brute_force.py --seed 3 --jobs 2
```

The native replay harness captures the parent once, stores the solver's
fp32 warm state and pool, and extracts a one-node state for all arms.
With the harness, baseline, seed, and experiment checkouts beside one
another, the tai35b commands are as follows.

```sh
replay=qap-scaling/experiments/cp_scaling_261117/tai35b_replay.py
harness=experiment_run_261117/replay.py
python3 "$replay" "$harness" capture qap-experiment-seed/latest/bnbdrv --cap 200
python3 "$replay" "$harness" run qap/latest/bnbdrv --arm scaling-tai35b-stock --cap 200
python3 "$replay" "$harness" run qap-scaling/latest/bnbdrv --arm scaling-tai35b-025 --cap 200 --exponent 0.25
python3 "$replay" "$harness" run qap-scaling/latest/bnbdrv --arm scaling-tai35b-050 --cap 200 --exponent 0.5
```

The tai35b parent was seeded with the tai50b production recipe, the comb
regions used in the tai35b solve, and a 200-iteration phase cap. The
timed node `19:12,18:26` is the mirror of the first two assignments
in the known optimum. Each replay retained 833 parent cuts and began
its CP phase from the same 110-iteration ordinary solve with certificate
278,672,093.3741. CPU to the target is the native trace's phase CPU
at the actual crossing check. The timed arms ran consecutively.

## Correctness checks

- Dense reference projection on both lifted axes, with different positive
  factors, a nonunit homogenizing factor, eliminated assignments, and an
  additional cut installed while scaled. Max iterate error was
  (3.33\times 10^{-16}). Final KKT residual was at most
  (3.74\times 10^{-17}).
- The continuation harness passed 29 checks against exact optima, including
  12 unscaled resumes from a scaled CP state and pool. It separated 247
  new cuts during the checks, exercised the metric switch with guarded
  extrapolation, and validated the saved primal's original box.
- `make test` passed. The coordinator's AWS scenario skipped because
  boto3 and moto were unavailable. The seven regression reference cases
  matched exactly. The default and two additional brute force seeds
  passed 4,536 solver checks each.
- An additional 192 exact-enumeration searches on random sizes 5–8
  passed with the CP experiment active at exponent .50 for 20 iterations
  per call. They covered both seeds 2 and 3 and produced 662 CP phases.
  The fixture used `cut.cp = 10:20`, `cut.cp_iters = 30`, and `cut.zone = 1`.
- Complete nug12 solves at exponent 0 and .25 both proved the known optimum
  578 with 37 visited nodes, including 24 DNN nodes and nine CP phases.
  The branch with scaling off reproduced the stock DNN trace after
  removing PID and CPU fields. The stock and branch traces also matched
  with the existing ordinary-DNN exponent .50 and a 40-iteration window.

## tai35b result

The shared certified target was the stock check-200 value
280,065,070.2172. Both scaled arms crossed it at check 200.

| Exponent | Certified bound | CP CPU to target | First plus CP CPU |
|---:|---:|---:|---:|
| 0 | 280,065,070.2172 | 16.28 s | 25.62 s |
| .25 | 280,066,383.9761 | 16.04 s | 24.07 s |
| .50 | 280,067,658.0827 | 16.47 s | 25.01 s |

The fixed first solve varied from 8.03 to 9.34 CPU seconds despite
identical iteration counts and certificates. The small CPU differences
in the phase do not meet the 15% promotion gate.

## tho40 pooled depth-3 result

The node was `26:10,8:0,2:32` on the transposed instance. The seed used
the production `tho40c.par` and the depth-2 parent
`26:10,8:0` with its ordinary and cut phases. The replay consumed
the same parent fp32 `Y/Z` and cut pool in each arm. The 110-iteration
first solve certified 237,871.8951 with 1,482 inherited cuts retained.
The numbers are a 200-iteration screen under shared-machine load.

```sh
harness=experiment_run_261117/replay.py
python3 "$harness" run scaling-tho40-stock tho40-d3-hard qap/latest/bnbdrv --cap 200
python3 "$harness" run scaling-tho40-025 tho40-d3-hard qap-scaling/latest/bnbdrv --cap 200 --env DNN_CP_DSCALE=0.25
python3 "$harness" run scaling-tho40-050 tho40-d3-hard qap-scaling/latest/bnbdrv --cap 200 --env DNN_CP_DSCALE=0.5
```

The shared target was 238,277.7695. Both scaled arms reached it at
check 200.

| Exponent | Certified bound | CP CPU to target | First plus CP CPU |
|---:|---:|---:|---:|
| 0 | 238,277.7695 | 30.44 s | 47.88 s |
| .25 | 238,282.4542 | 32.97 s | 50.04 s |
| .50 | 238,280.4580 | 30.25 s | 48.78 s |

Both pooled screens left the promotion gate unmet. Scaling is opt-in.
The measurements cover early phase progress at the stated depths.
