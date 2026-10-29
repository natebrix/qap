# A certified DNN bound in C++ (phase 1: library, driver, validation, benchmark)

**Short answer.** `latest/dnn.C` computes the same certified DNN lower bound as the Python
prototype. It is **8x faster per ADMM iteration at m = 28** (43 ms against 339 ms for dnn.py,
measured back to back on the same loaded machine). It also needs **about 6x fewer iterations
to certify**. On the 162 nug20/nug30 nodes of the two studies it makes **the same fathoming
decision as dnn.py on every node**, with **13-36x less CPU per node**:

* nug30 depth 2: 9.6 s per node on average, against 139 s.
* nug30 depth 3: 7.4 s per node, against 121 s.
* The 8 lowest-QPB nug30 misses: 11 s each, against 397 s.

The bound stays at or below the optimum on 1512 brute-force node bounds. In prototype mode it
reproduces dnn.py's certified values to 1e-7 relative. Results are bitwise reproducible, and
valgrind reports no errors or leaks. Nothing is wired into `bnb.C`.

## 1. What changed relative to dnn.py

The relaxation, the ADMM (OWX: R-step projection, Y-step clip, Z-step with gamma 1.618) and the
meaning of the certified value are the same as in dnn.py. Five things differ (details in the
`dnn.C` header comment):

1. **Products with V in O(m^4).** The basis of e-perp is the first m-1 columns of the Householder
   reflector H with H e_m = e/sqrt(m). So V'MV and VRV' are a reflection along each of the four
   tensor modes plus a 2x2 rotation of the coordinate pair (0, (m,m)). No V is formed. The
   R-step and the bound do not depend on the choice of orthonormal basis.
2. **Certificate without an extra eigendecomposition.** The certificate uses the R-step
   multiplier Z' = Z + beta (Y - VRV'). Its projection V'Z'V = beta (W - R) is negative
   semidefinite by construction. lambda_max(V'Z'V) is bounded rigorously by a Cholesky
   factorization of delta I - V'Z'V (Rump's criterion, about 7 ms at m = 28), and the bound pays
   (m+1) delta, using trace(R) = m+1 on the feasible set. dnn.py instead projects V'ZV with one
   more full eigendecomposition per check.
3. **A stronger polyhedral part** (Gilmore-Lawler on the residual, `DNN_POLY_GL`). Every column
   of a feasible Y lies in range(V). So the m x m block Y_(kl),(ij) has all row and column sums
   equal to x_ij = Y_0,(ij), and x is doubly stochastic. Minimizing <L+Z, Y> over that set
   instead of dnn.py's box gives M_00 + LAP(c), with c_ij = 2 M_0,(ij) + M_(ij),(ij) + the LAP of
   the (m-1) x (m-1) minor. It is never below the box value and has the same limit. Every LAP
   value used is the objective of a dual-feasible point, repaired explicitly, so inexact inner
   LAPs stay valid. On nug20 nodes this certifies a given value 1.5-6x sooner, for example
   2587 instead of 1704 after 50 iterations on node {0->19, 1->0}.
4. **R-step eigensolver.** W has hundreds of eigenvalues of each sign (about 70% positive), so
   a partial spectrum with dsyevr does not help. dsytrd alone costs 40 ms at K = 730. The default
   (`DNN_EIG_SUBSPACE`) works as follows:
   * It tracks the negative eigenspace of W from one iteration to the next. A degree-2 Chebyshev
     filter, CholQR2 and Rayleigh-Ritz run in single precision with BLAS-3.
   * A full ssyevd runs every 20 iterations, or when the block looks too small.
   * R is formed as W minus the computed negative part, so W - R stays exactly negative
     semidefinite and the certificate in item 2 is unaffected.

   Only the ADMM step is inexact, and iteration counts did not change in the tests. The exact
   modes (dsyevd as dnn.py, or ssyevd) are still available.
5. **Step size.** beta = 0.3 x (dnn.py's m/3). This certifies 2-3x sooner on hard nodes (nug30
   {4->10, 6->8}: 480 iterations instead of 1110; nug20 {1->7, 19->11}: 270 instead of 640).
   Easy nodes need 10-20 more iterations.

Also new:
* An early stop `DNN_HOPELESS`, after 200 iterations, when the certified bound extrapolated
  generously (5x a geometric tail) stays below the cutoff. It is a heuristic and the returned lb
  remains valid.
* An optional warm start (`DNNWarm`, `dnn_warm_child`): a child starts from its parent's
  multiplier Z, mapped like the objective.
* A prototype mode (`dnntest -P`: box part, dnn.py's certificate, dsyevd, beta 1, checks every
  50 iterations) that reproduces dnn.py.

Floating point: every certified value subtracts a margin built from standard error bounds with
generous constants. The parts are the Cholesky criterion, the Householder transform error
(64 m u ||Z||_F) and the summation error of the polyhedral part. The largest margin seen at
m = 28 is 1.5e-7 in objective units. This is not interval arithmetic.

## 2. Files, build and dependencies

**New files:**
* `latest/dnn.h`, `latest/dnn.C`: the library. There is no global state and the workspace is
  allocated per call (about 80 MB peak RSS at m = 28).
* `latest/dnntest.C`: the driver.
* `tests/dnn_test.py`: the validation.
* `bench/dnn_bench.py`: the benchmark.

**Modified files:**
* `latest/makefile`: new `dnntest`, `dnn.o` and `dnntest.o` rules, and `clean` also removes
  `dnntest`. The default `bnbdrv` build is unchanged.
* `Makefile`: new `dnntest` and `dnncheck` targets. `make test` is unchanged.
* `.gitignore`: adds `latest/dnntest`.

Solver sources are untouched.

**Linkage.** `latest/makefile` asks Python for the scipy-openblas32 wheel
(`python3 -c "import scipy_openblas32"`):
* If it is installed, `dnntest` links its `libscipy_openblas.so` with an rpath, compiles
  `dnn.C` with `-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS` (the wheel prefixes its symbols, e.g.
  `scipy_dsyevd_`), and uses 32-bit integers.
* Otherwise it links the system `-l:liblapack.so.3 -l:libblas.so.3`.
* Both can be overridden: `make dnntest DNN_LAPACK="..." DNN_LAPACK_DEFS="..."`. Run
  `make clean` after switching.

I installed the wheel in this container with `pip install scipy-openblas32` (into
`/usr/local/lib/python3.13/dist-packages`; no root or apt needed). numpy's own OpenBLAS64 is
bundled inside numpy and is not linkable. The reference LAPACK fallback was tested: same
iterations, 7x slower per iteration (40 ms against 5.7 ms at m = 18). The driver sets OpenBLAS
to one thread.

**Running the driver:**

```sh
make dnntest                      # or: make -C latest dnntest
cd latest
./dnntest -u 6124 -x 4:10,6:8 ../qaplib/nug30.dat   # node {4->10, 6->8}, cutoff 6124 - 1.95
./dnntest -P -x 0:19,1:0 ../qaplib/nug20.dat        # prototype mode (dnn.py's numbers)
./dnntest -N nodes.txt -u 2570 ../qaplib/nug20.dat  # many nodes, one "i:j,... [cutoff]" per line
./dnntest -v -v ...                                 # bound at every check + CPU profile
```

Other flags:
* `-e 0|1|2`: eigensolver.
* `-b`, `-l 0|1|2`: box part, or inner LAP mode.
* `-B`, `-g`: beta factor and gamma.
* `-H`: when the hopeless test starts.
* `-K i:j`: bound a child cold and warm-started.

The last line is `RESULT key=value ...`.

## 3. Validation (`python3 tests/dnn_test.py`, 1 min 54 s on one core; `make dnncheck`)

* **A. Brute force.**
  * Instances: 48 random instances of the four `tests/brute_force.py` families (n = 3-8). Each
    gets its root and 4 random subproblems with 1 to n-1 facilities fixed (m = 1-8), and the
    optimum comes from enumeration.
  * Variants: default, ssyevd, dsyevd, inner LAPs reduce and exact, box, and prototype. That is
    1512 node bounds.
  * lb <= optimum always. Over m >= 2, max(lb - optimum) is -1.65e-13.
  * No fathom was claimed with the cutoff at the optimum. With cutoff = optimum - 1, every
    fathom has lb > cutoff.
  * The bound is within 1e-3 of the optimum on 1487 of the 1512 node bounds: DNN is nearly
    exact at this size.
* **B. dnn.py agreement.** dnntest -P was compared with dnn.py at every check (every 50
  iterations): nug12 root, two nug12 subproblems, and nug20 nodes {0->19, 1->0}, {1->7, 19->11}
  and {1->1, 12->19, 16->15}. That is 44 values; the largest relative difference is 8.6e-8.
  The nug12 root gives 567.9802 in both, after 900 iterations, converged.
* **C. nug20 nodes with exact subproblem optima** (from the study's bnbdrv runs). The default
  bound stays below the optimum with the cutoff at 2568.05 and at the optimum itself. It fathoms
  exactly the 5 nodes the study fathomed. A second run reproduces every result bit for bit.
* **Fathoming decisions on all study nodes** (`bench/dnn_bench.py study`): identical to dnn.py's
  on all 162 nodes:

  | instance | set | fathomed by both |
  |---|---|---|
  | nug20 | depth 2 | 96/98 |
  | nug30 | depth 2, random | 34/36 |
  | nug30 | depth 2, lowest QPB | 0/8 |
  | nug30 | depth 3 | 19/20 |

* **valgrind** (separate `-O1 -g` build, every eigensolver, certificate and LAP mode, and
  m = 1, 2, 10, 12): 0 errors and no leaks.

## 4. Speed (one core; the machine was loaded by the other benchmarks throughout)

**CPU per ADMM iteration** (`bench/dnn_bench.py iter`, measured back to back, 40 C iterations
and 13 dnn.py iterations each). The dnn.py column includes one certificate per 13 iterations
(about 7%).

| m | node | C default (subspace) | C ssyevd | C dsyevd | dnn.py | speedup |
|---|---|---|---|---|---|---|
| 18 | nug20 {0->19, 1->0} | 5.7 ms | 9.7 | 13.4 | 26.9 ms | 4.8x |
| 22 | nug24 {0->0, 1->1} | 13.9 ms | 24.3 | 44.3 | 73.6 ms | 5.3x |
| 25 | nug27 {0->0, 1->1} | 21.1 ms | 50.1 | 75.4 | 146.6 ms | 7.0x |
| 28 | nug30 {4->29, 6->7} | 43.4 ms | 103.8 | 136.3 | 338.6 ms | 7.8x |

dnn.py ran at 135 ms per iteration at m = 28 in the nug30 study, so the machine is now about
2.5x slower than it was then. On that basis the C bound would be roughly 17-20 ms per iteration
at m = 28; this is an extrapolation, not a measurement.

Breakdown at m = 28: eigensolver about 60%, the two transforms about 25%, Y/Z-step about 6%,
certificates about 4% (Cholesky, plus exact inner LAPs only near the cutoff).

**CPU per node**, cutoff optimum - 1.95. The dnn.py figures for nug20 and nug30 come from the
studies, which ran on the same machine under lighter load, so these ratios understate the
speedup. For nug24 and nug27 they are prototype-mode iteration counts times dnn.py's measured
CPU per iteration.

| set (m) | nodes | fathomed C / dnn.py | C mean (median iters) | dnn.py mean (median iters) | ratio |
|---|---|---|---|---|---|
| nug20 depth 2, all (18) | 98 | 96 / 96 | 0.62 s (70) | 8.35 s (500) | 13x |
| nug24 depth 2, random QPB-open (22) | 8 | 8 / 8 | 1.75 s (70) | 46 s est. (425) | 27x |
| nug27 depth 2, random QPB-open (25) | 8 | 7 / 7 | 3.4 s (75); 3.1 s excl. the miss | 132 s est. (650), excl. one 4000-iteration miss | 43x excl. the miss |
| nug30 depth 2, random (28) | 36 | 34 / 34 | 9.6 s (130 if fathomed) | 139 s (750) | 14.5x |
| nug30 depth 2, 8 lowest QPB (28) | 8 | 0 / 0 | 11.0 s (hopeless at 200) | 397 s (run to convergence) | 36x |
| nug30 depth 3, 3 samples (27) | 20 | 19 / 19 | 7.4 s | 121 s | 16x |

Per-node percentiles (0/10/25/50/75/90/100) for fathomed nug30 depth-2 nodes:
* CPU: 2.8/3.1/4.1/7.0/9.7/26.9/29.5 s, against 67/71/88/102/137/233/375 s for dnn.py.
* Iterations: 60/70/80/130/170/430/500, against 500/550/650/750/1000/1700/2750.

Where the gain comes from at m = 28:
* Per iteration: 7.8x. Single precision gives 1.3x over dsyevd, and subspace tracking 2.4x more.
* Iterations to certify: about 6x at the median. The GL certificate and beta = 0.3 contribute
  about equally.

Warm start (`-K`), on 2 fathomed nug30 depth-3 children of the miss {4->5, 6->13}: 40 instead of
100 iterations, and 30 instead of 70. On the unfathomable child: 450 instead of 540. Starting Y
from the parent (conditioned on x_ij = 1) made a child slower, so only Z is carried over. That is
3 nodes: a hint, not a measurement.

## 5. Open issues

* **Cross-machine reproducibility.** OpenBLAS (DYNAMIC_ARCH) picks kernels per CPU. The bits of
  a DNN bound, and so a fathom decision that falls within ~1e-6 of the cutoff, can differ between
  machines. That would break the "node counts are machine-independent" rule of CONTRIBUTING.md
  once the bound is in the search. Pinning `OPENBLAS_CORETYPE` would probably fix this; I have
  not tested it. Results are bitwise identical on one machine.
* **The hopeless stop is a heuristic.** With the defaults it never gave up on a fathomable node
  in 178 nodes. During tuning, though, it did fire wrongly at 200 iterations with a degree-1
  filter, where the ADMM was slowed down. For nug30 depth 2, where one wrong give-up costs hours
  of QPB, use a later or no hopeless stop and select nodes by QPB (section 6).
* **The subspace eigensolver makes the ADMM inexact.** Convergence with inexact projections is
  only shown empirically: same iteration counts on every node tested, and a full ssyevd every
  20 iterations. Certificates are unaffected. `eig_mode = 1/0` gives exact projections.
* **Rigor of the margins.** They come from standard error bounds with generous constants: Rump's
  criterion for the Cholesky, 64 m u ||Z||_F for the transforms, and 8 (m+2) u sum|M| for the
  polyhedral part. This is not verified interval arithmetic. The margins are 1e-9 to 1.5e-7 in
  objective units, far below the 1.95 fathoming tolerance and the 0.001-0.1 certified margins
  seen.
* **Tested data.** Only nug instances and the small random families were tested. Large-magnitude
  data (tai*b) and instances where both A and B are asymmetric are untested; L is symmetrized as
  (A kron B + A' kron B')/2, which is valid. With `STORE_J_MATRIX` (off) the 1e5 penalties in C
  would inflate the scaling.
* **Remaining cost.** About 60% of an iteration is in the eigensolver. Further options: a
  block-Krylov or LOBPCG refinement, a lower refresh rate, and fusing the expand transform into
  the Y-step. The certified lb still trails the primal on close nodes. Extrapolating the
  multiplier between checks, then projecting it, is untried.
* **kra30a/kra30b root bounds (Oct 29, held out of the report's Table 4).** With the default
  parameters (beta = 0.3 m/3 on the problem scaled to max |L| = 1, 4,000 iterations, no cutoff)
  the ADMM does not converge on kra30a/kra30b: the certified root bounds are 17,368 and 17,572,
  far below QPB (about 68,600 and 69,000). The bounds are valid, just useless; the relaxation
  itself is not expected to be weak there. Suspects: beta, tuned on nug, against kra's distance
  scale (3-D grid, unit costs 50/50/115) and sparse flows. To do: rerun with `dnntest -B` (beta
  factor) and more iterations, and check the primal/dual residuals (`-v`).
* **Per-call allocation.** Each call allocates about 80 MB at m = 28. That is negligible against
  seconds of work, but a reusable workspace would help if the bound is called millions of times
  (it should not be).

## 6. Recommended wiring into bnb.C (phase 2)

**Where.** In `bnb()`, right after the QPB bound and `node.bound = max(...)` (bnb.C around
line 607), before the fathom/branch test. `qr` (reduced A, B, C, shift) is already there:

```cpp
if (dnn_wanted(param, node, m, bound, s->getObjective())) {
  DNNResult r = dnn_bound(qr->A, qr->B, qr->C, qr->shift, s->getObjective() - DELTA, dnnpar);
  stats->dnn(node.depth, r);                          // calls, fathoms, CPU per depth
  node.bound = std::max(node.bound, r.lb);
  if (r.lb > s->getObjective() - DELTA) continue;     // fathomed (count as such)
}
```

**Keep the QPB value `bound` for `branch(..., bound + qpbp->lapU, ...)`.** Child elimination by
reduced costs U is valid only relative to the QPB solution that produced U. Passing the DNN lb
there would be invalid. The DNN lb can still raise `node.bound`, which the children inherit.

**Which nodes.**
* By depth and size, set per level in the params file (as a new file, not `default.par`, per
  CONTRIBUTING.md):
  * nug30-class: depths 2 and 3.
  * nug24-28: depth 2.
  * Not depth 1 (DNN sits 70-120 below the optimum there), and not m <= 18. At nug20, depth 2
    now costs 61 s for all 98 nodes plus about 60 s of QPB for the 2 surviving subtrees, against
    about 180 s QPB-only: a small win. Depth 3 is still a loss.
* By QPB, as in the nug30 study. The node's QPB predicted the DNN outcome perfectly: DNN is
  about QPB + 480-580 at depth 2. Skip nodes whose relative gap (inc - QPB)/inc is above a
  per-depth threshold (about 0.094 at nug30 depth 2). Run the rest with `hopeless_after` larger
  or 0, because there a wrong give-up costs a whole subtree.
  * With the default hopeless stop, a miss costs only about 11 s at m = 28, so a simpler rule
    (DNN on every node at the chosen depths) also works. The QPB filter then saves about 40 x
    11 s per nug30 depth-2 frontier.
* Estimated nug30 hybrid cost, using the study's counts and this code's measured per-node means:
  * DNN at depth 2: 185 x 8 s + 40 x 11 s, about 0.5 CPU-h.
  * DNN at depth 3 under the 40 survivors: 1096 x about 10 s, about 3 CPU-h.
  * Below the about 25 depth-3 survivors: DNN at depth 4, about 650 x 7 s, about 1.3 CPU-h, plus
    QPB below the few depth-4 survivors.
  * In total about 5-10 CPU-h, against 60-120 CPU-h with dnn.py and 33-120 CPU-days QPB-only.
    This is an estimate built on the study's counts and the solver's subtree estimates, not a
    measurement.

**Warm starts.**
* When DNN misses at depth d, keep the node's `DNNWarm` (`save` argument; 4.9 MB at m = 28).
  Hang it on the children like the Frank-Wolfe `WarmStart` (refcounted, released with the last
  child).
* Start each child from `dnn_warm_child(parent, i_local, j_local, child)`.
* Store Z at most on the active depth-2 and depth-3 path, not in state files: it is only a warm
  start and does not affect validity. If it is not stored, resumed runs give different iteration
  counts but the same bounds' validity.

**Other wiring notes.**
* Call `dnn_set_blas_threads(1)` once at startup. `bench/parallel.py` uses processes, so that
  is enough.
* Report DNN calls, fathoms and CPU per depth in the results block.
* Add a brute-force run of the solver with DNN enabled at every depth (m <= 8 is fast) to
  `tests/brute_force.py`.
