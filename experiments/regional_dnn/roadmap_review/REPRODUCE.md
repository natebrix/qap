# Reproduce the small research checks

These are diagnostics, not a wil50/sko49 performance implementation. Use a separate checkout of cognac commit `ca50467835d640a252f1f4d687618f593a395d7c` for headers and cached data. The surrounding research branch deliberately retains an older solver baseline; do not compile these diagnostics against that older `dnn.h`.

From this directory, set `QAP_REPO` to the pinned cognac checkout. Python checks require NumPy; cache analysis also requires SciPy.

```bash
QAP_REPO=/absolute/path/to/pinned-cognac
python3 bound_audit_checks.py
python3 tree_cache_analysis.py "$QAP_REPO"
```

Expected exact checks: 49 positive/16 zero PSD pivots, regional triangle maximum -1/16, conditioned five-point violation 1/8, and 600 parent-child objective identities. Cache output is retained in `tree_cache_analysis_results.txt`; it analyzes existing measurements and surrogate costs, not new solver runs.

## Same-node continuation

The included `speed_dnn_experiment.C` is an isolated copy of the pinned DNN implementation, with an environment-switched Y reload and the original LowerSum definitions appended for standalone linking. Do not also link `assign.C`. This is not a production patch and does not alter the child warm-start policy.

Set `SCIPY_BLAS` to an installed LP64 SciPy OpenBLAS shared library exporting the `scipy_`-prefixed symbols. The tested library was `libscipy_openblas-6cdc3b4a.so` in the SciPy wheel's `scipy.libs` directory. Other BLAS distributions require matching link flags and symbol-prefix definitions.

```bash
SCIPY_BLAS=/absolute/path/to/scipy.libs/libscipy_openblas-6cdc3b4a.so
g++ -O3 -march=native -ffp-contract=off -DBNB_STANDALONE \
  -DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS \
  -I"$QAP_REPO/latest" -I"$QAP_REPO/src/MESCHACH" \
  speed_dnn_experiment.C speed_continue.cpp "$SCIPY_BLAS" \
  -Wl,-rpath,"$(dirname "$SCIPY_BLAS")" -o /tmp/qap_speed_continue
/tmp/qap_speed_continue "$QAP_REPO/qaplib/nug20.dat" 2
/tmp/qap_speed_continue "$QAP_REPO/qaplib/nug20.dat" 1
/tmp/qap_speed_continue "$QAP_REPO/qaplib/nug30.dat" 2
```

Mode 2 uses the tracked subspace, mode 1 the full single-precision spectrum. The harness compares 200+100 iterations with Z-only and YZ restarts against uninterrupted 300. In the full-spectrum nug20 control, YZ and uninterrupted bounds/primal values matched to printed precision; tracked-subspace state is not yet retained, so those trajectories differ. Timings are illustrative only; no wil50 speedup is inferred.

```bash
g++ -O3 -march=native -ffp-contract=off speed_yz_microbench.cpp \
  -o /tmp/qap_speed_yz_microbench
/tmp/qap_speed_yz_microbench
```

This isolated bookkeeping check verifies Y/Z and check-residual equality under the tested conditions. Its timing is noisy and is not an end-to-end speed claim.

## Validation during packaging

The exact bound audit and cache analysis were rerun; the cache output matched the retained text exactly. The standalone continuation harness was rebuilt and the full-spectrum nug20 control rerun. The research update changes no production solver source, so the full solver regression suite was not repeated for this documentation bundle.
