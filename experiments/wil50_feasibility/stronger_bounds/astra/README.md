# Reproducible stronger-bound probes

These are research drivers, not changes to the branch-and-bound solver. See the accompanying research memo for proofs, interpretation, and limitations. The regional order-two SDP is formulated in the memo but is **not implemented** in this bundle.

## Source snapshot and dependencies

- Repository: https://github.com/natebrix/qap
- Commit: `cd43c7645d8af01a2ec22aa0ae9bba9fafea7022`
- C++ compiler, the repository's Meschach library, and BLAS/LAPACK as required by its makefile.
- Python: NumPy and SciPy; SymPy for exact witness verification.
- `region_probe.C` uses the existing `DNNParams::L_add` interface and returns the existing lower certificate.

The driver requires residual order at least two. It preserves qap's native index order: X[grid location, flow facility] for wil50. Fixes are in the repository's original zero-based matrix-index convention. The Python geometric separators are specifically for wil50's 5×10 grid. Do not apply their coordinate formula to another instance without changing it.

## Build

Build the repository normally first, including Meschach. With this bundle extracted to an absolute directory, run from `qap/latest`:

```bash
make -f makefile -f /absolute/path/wil50_research/probe.mk \
  RESEARCH_DIR=/absolute/path/wil50_research region_probe
```

Pass the same BLAS/LAPACK make variables used for the normal solver build. In this environment, the SciPy wheel's 32-bit OpenBLAS library was linked using `DNN_LAPACK_DEFS='-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS'` and the library's absolute path and rpath. The existing makefile also supports Accelerate on macOS and `scipy-openblas32`.

## Exact witness

```bash
python3 verify_witness.py
```

This reconstructs the rational DNN witness, verifies all lifted equalities and atomic Boolean triangles, obtains its exact eigenvalues, checks the regional cut coefficient construction on 720 permutations, and verifies the degree-four Boolean SOS identity.

`toy6.dat` has both data matrices equal to the adjacency matrix of K3,3. Its optimum is 10; DNN is 9. To recreate the modified toy objective:

```python
from region_cuts import cut_matrix
G,h = cut_matrix(dict(locations=[0,1,2], facilities=[0,1,2], type=3),6)
(8*G).tofile('toy6.L.bin')
```

```bash
PROBE_DOUBLE=1 PROBE_ADMM=1 ./region_probe toy6.dat - toy6_base 2000
PROBE_DOUBLE=1 PROBE_ADMM=1 ./region_probe toy6.dat - toy6_cut 2000 toy6.L.bin -8 -
```

## wil50 baseline and fixed-multiplier trials

The difficult sampled node is `34:40,28:22,14:29`. Replace `/path/to/wil50.dat` below with the file from the pinned repository. Work in this directory so checkpoints and logs have the expected prefixes.

```bash
./region_probe /path/to/wil50.dat '34:40,28:22,14:29' node1 1500 > node1.log
python3 region_cuts.py node1 --fixes '34:40,28:22,14:29' --mode threshold --keep 100
python3 build_trials.py node1 --fixes '34:40,28:22,14:29'

PROBE_ADMM=1 PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_refined 700 - 0 node1.warm.bin > node1_refined.log

PROBE_ADMM=1 PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_tri100_mu2 700 node1_tri100_mu2.L.bin -46 node1.warm.bin \
  > node1_tri100_mu2.log

PROBE_ADMM=1 PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_tri100_mu10 700 node1_tri100_mu10.L.bin -230 node1.warm.bin \
  > node1_tri100_mu10.log
```

The shift values above correspond to the retained experiment's 23 type-3 cuts. If another machine selects a different set of nearly tied cuts, use the shifts printed by `build_trials.py` instead.

The larger-pool trials use PRSM and a matched uncut PRSM continuation:

```bash
python3 build_trials.py node1 --fixes '34:40,28:22,14:29' --keep 1000 --multipliers 1,3
PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_prsm_control 700 - 0 node1.warm.bin > node1_prsm_control.log
PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_tri1000_mu1 700 node1_tri1000_mu1.L.bin -261 node1.warm.bin \
  > node1_tri1000_mu1.log
PROBE_THREADS=2 ./region_probe /path/to/wil50.dat \
  '34:40,28:22,14:29' node1_tri1000_mu3 700 node1_tri1000_mu3.L.bin -783 node1.warm.bin \
  > node1_tri1000_mu3.log
```

Again, use generated shifts if cut selection differs. The controls and trials all restart from the same original multiplier; retain the original 1,500-iteration bound if a continuation returns a weaker certificate. The prototype warm-start interface saves Z but resets Y, so a warm continuation is not identical to an uninterrupted solve.

By default, the driver uses PRSM parameters gamma=gamma1=1, degree-four subspace tracking, exact inner GL assignment solves, one BLAS thread, and certification every 25 iterations. `PROBE_ADMM=1` changes the dual updates to gamma=1.618, gamma1=0. `PROBE_THREADS` controls BLAS threads; `PROBE_DOUBLE=1` chooses double-precision full eigendecomposition for small validation problems. The driver always saves the final approximate Y and the warm-start multiplier. Only `RESULT ... lb=...` is the lower bound.

The reported `cpu` field is inherited from the repository and uses its timing convention. Do not equate that field to total wall time or total CPU consumption of all BLAS threads in concurrent runs.

## More separation

```bash
python3 region_cuts.py node1 --fixes '34:40,28:22,14:29' --mode mixed --keep 100
python3 region_cuts.py node1 --fixes '34:40,28:22,14:29' --mode rectangles --keep 100
python3 triplet_separate.py node1 --fixes '34:40,28:22,14:29' --facilities 4,25,46
```

The triplet separator takes **local free facility indices**. For this node, `4,25,46` correspond to original B-matrix indices `4,26,49`. Its `.npz` file stores integer potentials and denominator 1024; the `.G.bin` file stores the symmetric lifted cut with right-hand side zero. Every injective placement is checked in integer arithmetic after rounding and repair.

`.Y.bin`, `.warm.bin`, `.L.bin`, and `.G.bin` are native little-endian float64, row-major arrays. Warm files begin with one float64 objective scale. The large checkpoint and objective files are omitted from the archive; regenerate them using these commands. The exact 37×37 witness and small toy instance are included.

Logs and JSON files distinguish diagnostics from certificates. An infeasible approximate Y can violate constraints because of numerical error. A violated cut, even when valid, is not a demonstrated bound improvement until the modified problem produces a better lower certificate. Always compare against an uncut control and the strongest prior certificate.
