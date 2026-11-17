# Independent certificate verification

Three CP implementations ran with their experimental features enabled on four asymmetric integer QAPs of order 7. Two fixtures were evaluated at the root and two after fixing facility 0 to location 4. Both lift orientations were exercised for each fixture. Every finite certificate checked was at or below the optimum found by enumeration.

| Track | Activated solves | Finite values checked | CP lazy checks | Solves with lazy recovery | Largest precise bound minus optimum |
|---|---:|---:|---:|---:|---:|
| Coarse-label cuts | 8 | 33 | 366 | 1 | −2.584e−7 |
| Congruence scaling | 8 | 33 | 362 | 1 | −1.598e−7 |
| Regional equalities | 8 | 33 | 384 | 1 | −9.610e−8 |

The finite counts include trace values, the CP summary, the quiet certificate, and recovered values printed at four decimal places. Some lines summarize the same certificate. The run checked 24 active coarse rows on every fixture, nonzero cardinality and second-moment residuals on every equality fixture, and four scaled-to-plain transitions at iteration 23. Dumped primal states passed the original lifted box and homogenizing-entry checks. [Per-case data, commands, and build fingerprints](certificates.json)

## Reproduce

From the results checkout, with the three built experiment checkouts alongside it, run the commands below. Build each checkout with its usual `make -C latest dnncuts` configuration first. The verifier links temporary copies of the driver to the existing solver objects. It sets the `dnn.dscale0`, `dnn.lazycert`, and `cut.cp_lazycert` policy keys and uses a CP convergence tolerance of 0.005 to exercise lazy recovery. Each solve uses a cutoff 1,000 above the exact optimum so it reaches CP and cannot legitimately fathom. The verifier leaves the solver checkouts unmodified.

```bash
python3 experiments/bounds_gates_261117/verification/verify_exact.py \
  --coarse ../qap-coarse --scaling ../qap-scaling --equalities ../qap-equalities
python3 experiments/bounds_gates_261117/verification/verify_geometry.py
```

The second command uses NumPy and passed 2,690 comparisons. It builds explicit dense lifted normals to check scaled overlap kernels, regional identities on all 120 permutations at order 5, and the exact full norm and dyadic RHS of a cross-label cut on groups of unequal capacity.

## Certificate review

- **Coarse-label cuts.** The loader recomputes the RHS in integer arithmetic over every capacity-feasible local coarse assignment. The symmetric matrix uses two half-coefficients per off-diagonal placement pair. Sparse projection, regional overlap corrections, and certificate assembly use this same normal. The exact fixtures exercise both signs of cross-label assignment identities, including RHS ±1. They check validity and activation. The benchmark separator is responsible for finding inequalities that add geometric strength. Coarse rows stay at the selected node.
- **Congruence scaling.** The iterative coordinates use the scaled clipping box and transformed normals. The certificate uses the original cut normal with multiplier `scale * beta * nu` and transforms the PSD multiplier back via `D * Z_t * D`. Saved Y and Z and pending extrapolation states convert at the end of the scaling window. Exact tests used nonuniform scaling with exponent 0.5 and homogenizing factor 2.5. The review applies to the tested fixed metric and settings. Existing adaptive beta residual balance uses a mix of original and scaled units if enabled.
- **Regional equalities.** The cached pseudoinverse uses the lifted Frobenius Gram metric and handles dependent rows. Capacities come from each region's self-overlap, which counts distinct indices exactly as the projector does. The Dykstra corrections are signed. Regular, lazy, and quiet certificates include the signed dense term, its absolute contribution to the rounding allowance, and the signed cardinality RHS. The solver's separate projection self-test covers dependent equations, repeated imported indices, and negative RHS multipliers against dense direct calculations.
- **External terms and checkpoint provenance.** On the stock `nug12` driver, adding an external constant contribution of −0.75 lowered the first round by 0.75 and left the CP and quiet bounds unchanged. The coarse experimental driver rejects external terms in CP mode and integrates the validated coarse cuts directly. The warm replay driver's checkpoint includes a previously certified bound supplied by the local file. Use trusted checkpoints from the corresponding instance and node. The independent exact verifier uses newly computed CP and quiet certificates throughout.

These checks cover the three experimental arms independently. They establish certificate behavior on the stated fixtures and do not estimate speed, whole-tree performance, or the strength of generated cuts.

## Stock external-term regression

At stock commit `ba47c211`, the external objective `−Y[0,0] − (−0.25)` contributes exactly −0.75 because the projector fixes `Y[0,0]` at 1. The CP certificate was 459.9095451953 and the quiet certificate was 486.0732292342 in both stock runs. Only the first bound shifted. Run from a stock checkout with `dnncuts` built.

```bash
fixture_dir=$(mktemp -d)
python3 - "$fixture_dir/external.bin" <<'PY'
from array import array
import sys
term = array("d", [0.0]) * (145 ** 2)
term[0] = -1.0
with open(sys.argv[1], "wb") as stream:
    term.tofile(stream)
PY
experiments/dnn_cuts/dnncuts -I 10 -k 10 -cpsep 5:10 -cpiters 10 -cpquiet 10 \
  -grid 3:4 -side 1 qaplib/nug12.dat
experiments/dnn_cuts/dnncuts -I 10 -k 10 -cpsep 5:10 -cpiters 10 -cpquiet 10 \
  -grid 3:4 -side 1 -Ladd "$fixture_dir/external.bin:-0.25" qaplib/nug12.dat
```
