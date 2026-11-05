# Removing Meschach from the solver

What it would take to remove the Meschach matrix library (`src/MESCHACH/`)
from the branch-and-bound solver in `latest/`: what the solver uses, what
replaces each piece, what that does to node counts, how much time Meschach
takes, and a staged plan. Research only: nothing in `latest/` or `src/` was
changed. Prototypes and measurement logs are in
[`docs/meschach_removal/`](meschach_removal/).

**Status (Nov 5, late): stage 0 is done** (`src/MESCHACH/makefile`: `CC ?= cc`,
`-std=gnu17 -ffp-contract=off`; `bench.py` honors `MESCHACH_CC`; `make` rebuilds
Meschach when its makefile changes). Checked on the cloud box: the quick
suite matches `bench/reference/` exactly when built with gcc 13, gcc 11, clang
18, and clang or gcc with `-march=native` (FMA hardware). Every gate passes.
The control, clang with FMA contraction on, reproduces the M2 Pro's numbers
(nug12 264, nug14 2,317, nug15 1,935 nodes) and fails the check. Confirmation
on an M2 itself is pending.

Every number is marked **[M]** (measured, with the log or command) or
**[E]** (estimated). Measurements were taken on the 4-vCPU container
(Intel Xeon @ 2.10 GHz, g++ 13.3, clang 18.1, scipy-openblas32), at most two
compute processes at a time, at `cognac` 45d05a2.

## Summary

- The solver calls 17 Meschach functions [M]. Fourteen are allocation, copying
  and trivial arithmetic. Three (`Hfactor`, `makeHQ`, `trieig`) are the
  symmetric eigendecomposition behind the QP bound (`my_svd` in
  `latest/linalg.C`), which runs twice per bounded node and once per lookahead
  child. The DNN code uses only the `MAT` struct's fields.
- On a QPB-only nug20 run, Meschach takes about 20% of the CPU time [M]: 13.3%
  in the eigendecomposition and 7% in allocation and copying. In DNN runs the
  share drops to 5% with `-x 2:2` and 0.07% with `-x 2:3` [M]. It is lower
  still in the nug30 record setting, where the DNN bound takes 99% of the
  time [E].
- LAPACK is no faster at these sizes (n = 10 to 16): dsyev took 12.8 µs per
  call against Meschach's 13.3, and dsyevr 17.8 [M]. It also changes the
  search. On the standard suite every instance moved, with node ratios
  between 0.870 and 1.114 [M].
- A plain-array port of Meschach's own eigen code, done in this study,
  matches Meschach bit for bit [M]. It agreed on 656 of 656 sampled
  matrices, on the standard suite and on the nug20 tree, and the brute-force
  test passed. It is 1.5x faster than Meschach on the eigendecomposition and
  6.1% faster end to end on nug20, in two paired runs [M].
- Compiling Meschach with FMA contraction (clang's default on arm64) gives
  exactly the M2 Pro discrepancy that RUNNING.md reports (nug15 root bound
  984.79, nug14 +6%, nug15 -1.5%, nug16b -1.5%) [M, simulated on x86].
  Adding `-ffp-contract=off` to Meschach's CFLAGS should make the M2 match the
  Xeon reference [E, not run on an M2].
- Meschach does not compile under C23 with either gcc 13 or clang 18 [M].
  GCC 15 defaults to C23, so a newer toolchain will break the build unless
  `-std=gnu17` is passed.
- Plan: stage 0 (build flags), stage 1 (eigen port), stage 2 (data
  structures) and stage 3 (build cleanup) can each be bit-identical on x86.
  Together they take 4 to 5.5 engineering days [E]. Recommendation: do
  stage 0 now (half a day). Do stage 1 when convenient (1 to 1.5 days, ~6%
  faster QPB-only runs, no change in results). Do stages 2 and 3 later, as
  hygiene. Do not move the QP bound's eigensolver to LAPACK.

## 1. Inventory

### What is linked

The Meschach symbols each solver object references, taken from `nm` on the
`.o` files
([`symbols_by_object.sh`](meschach_removal/symbols_by_object.sh),
[log](meschach_removal/logs/symbols_by_object.txt)) [M]:

| object | Meschach functions referenced |
|---|---|
| linalg.o | `Hfactor makeHQ trieig m_transp _m_copy m_get m_free m_resize m_zero v_get v_free v_resize` |
| qpbnd.o | `_in_prod _m_copy _v_copy m_get m_free m_resize px_get px_free px_rows v_get v_free` |
| qpbnd2.o | `m_get m_free px_get px_free px_rows v_free` |
| branch.o | `m_get m_free m_resize v_get v_free` |
| console.o | `_in_prod _m_copy m_get m_free m_resize m_zero v_get v_free` |
| assign.o | `_m_copy m_get m_free m_resize` |
| bnb.o | `m_get m_free m_resize` |
| fw.o, qap.o | `_m_copy m_get m_free` |
| state.o | `m_get m_free` |
| svd.o | `_m_copy _v_norm_inf ev_err get_col get_row m_free m_ident mem_stat_reg_list v_resize v_zero` (dead code, below) |
| dnn.o, dnncut.o, lap.o, dlap.o, util.o, driver.o, others | none |

Without `svd.o` that is 17 distinct functions [M]. `m_copy`, `v_copy` and
`in_prod` are macros over `_m_copy`, `_v_copy` and `_in_prod`. `M_FREE`,
`V_FREE` and `PX_FREE` are macros over the `*_free` functions.
`dnntest.o`, `tests/dnn_continuation.C` and `experiments/dnn_cuts/dnncuts.C`
reference none [M]. They link `meschach.a` only because they link the
solver's objects.

`latest/svd.C` (586 lines) is dead code [M]. It holds a modified copy of
Meschach's SVD (`svd`, `bisvd`, `hhvec`), and nothing calls `svd()`.
`my_svd` replaced it with a symmetric eigensolver when sign guessing was
found to give invalid bounds (see the comment in `linalg.C`). Its `hhvec`
is a C++ symbol (`_Z5hhvecP3VECjPdS1_`), so Meschach's `Hfactor` still uses
Meschach's own `hhvec`.

### Identifiers by category

Counts from [`inventory.py`](meschach_removal/inventory.py)
([log](meschach_removal/logs/inventory.txt)) [M], over 27 source files
outside `src/MESCHACH` (1,784 occurrences, 227 of them in the dead `svd.C`):

| category | identifiers (occurrences) | where |
|---|---|---|
| data structures | `MAT` 574, `VEC` 99, `PERM` 4. Field access: `->me` 502, `->ve` 151, `->pe` 4. `->max_*` and `->base` are never used. | qpbnd.C, branch.C, fw.C, linalg.C, headers. dnn.C/dnncut.C/dnn.h use `MAT` 15 times, reading `->m` and `->me` only. |
| allocation | `m_get` 73, `v_get` 28, `px_get` 3, `M_FREE` 111, `V_FREE` 38, `PX_FREE` 3, `m_resize` 34, `v_resize` 8 | everywhere in the bound and branching code |
| copy and init | `m_copy` 25, `v_copy` 2, `m_zero` 5, `v_zero`/`m_ident`/`get_col`/`get_row` (svd.C only) | assign.C, fw.C, qpbnd.C, linalg.C |
| arithmetic | `in_prod` 10, `px_rows` 5, `m_transp` 1, plus svd.C only `v_norm_inf`/`v_mltadd` | qpbnd.C (the S,T updates), qpbnd2.C, linalg.C |
| linear algebra | `Hfactor`, `makeHQ`, `trieig`, one call site each, in `my_svd` | linalg.C |
| I/O and errors | `m_output`, `v_output` only in comments or dead code. `error(E_*)` 18 and `MACHEPS` 4, all in svd.C | none live |

The solver already does its own matrix arithmetic on `MAT` storage:
`my_m_mlt`, `mtvm_mlt`, `vtav`, `vavt`, `permrows` and `permcols` in
`linalg.C`, and the Frank-Wolfe gradient in `fw.C`. None of it calls
Meschach.

Outside `latest/`:

- `tests/dnn_continuation.C` and
  `experiments/wil50_feasibility/speed/beta_sched/beta_change_check.C` build
  a `MAT` by hand (`MAT a{}` plus a row-pointer vector) and use its fields.
- `experiments/wil50_feasibility/new_cuts/facility/src` has `MAT *` in
  signatures, and its `build.sh` links `/home/user/qap/src/MESCHACH/meschach.a`
  by absolute path. `experiments/wil50_feasibility/roadmap/convergence2/build.sh`
  builds Meschach.
- `bench/bench.py` (`build()`, used by `--build` and `--rev`) runs
  `make basic` in `src/MESCHACH`, with the makefile's default `CC = clang`.
  Without clang it fails.
- `infra/phase0/bench.sh` installs clang or falls back to `CC=gcc`.
- `Makefile`, `latest/makefile` (`INCLUDES`, `LDLIBS`) and RUNNING.md ("Build").
- `python/` and `coord/` do not use Meschach.

### Hot paths

Call counts and sizes on nug20, QPB only (`default.par`, `-u 2570`,
526,205 nodes), from the timing build described in section 4 [M]:

| path | when | Meschach work per call |
|---|---|---|
| `qpbnd_parametric` → `init_ST` (qpbnd.C) | every bounded node | 2× `my_svd` of size m−1 (m = free facilities), 2× `px_get`/`px_rows`, 4× `v_get`, 4× `m_get` |
| `init_w0` + `qpbnd2` (qpbnd2.C, via branch.C rules 3, 4 and 7) | every lookahead child | 1× `my_svd` of size m−2, `px_rows`, `m_get` for G and friends |
| `update_STG` (qpbnd.C) | each S,T update inside the parametric bound | `v_copy`, `in_prod`, `m_copy`, `m_get` |
| `reduceQAP` (assign.C), `bnb.C` | every node | `M_FREE` + `m_copy`/`submatrix` (`m_get`) of A, B, C, and `m_resize` of X and U |
| `fw` (fw.C) | every bound | `m_get`/`m_copy` of the best X and the gradient step |
| Gilmore-Lawler (`glbnd`/`leader_matrix`, bound 4; `qp_glb_bnd`, bound 2) | per node when selected; `default.par` uses bound 1 | `m_get`/`M_FREE` of the leader-matrix work arrays; no eigendecomposition |
| qap.C (read), state.C (resume), console.C (`-c`), `pbnd` (bound 5) | setup or off by default | `m_get`, `m_copy`; `pbnd` also 2× `my_svd` |

There were 1,311,622 `my_svd` calls in all (2.5 per node), 87% of them
at n = 11 to 14. There were 18.3M `m_get` and 24.3M `m_free` calls, about
35 and 46 per node [M,
[log](meschach_removal/logs/nug20_timer.err)].

## 2. Replacements

| Meschach | replacement | same bits? |
|---|---|---|
| `MAT` {m, n, me, base, max_*} | A solver-owned struct with the same `m`, `n` and `me` fields (row pointers into one contiguous block), so the 502 `->me` accesses and every `MAT *` signature, including `dnn.h`'s `const MAT *`, stay as they are. About 150 lines with the functions below [E]. | yes (storage only) |
| `VEC` {dim, ve} | Same treatment. `pbnd` writes `gsort->dim--` directly, so `dim` must stay a plain field. | yes |
| `PERM`, `px_get`, `px_rows`, `px_free` | `int` arrays and the existing `permrows()` (linalg.C). Three call sites. | yes (copies) |
| `m_get`, `v_get` | Allocate zeroed memory. Meschach uses calloc, and callers rely on it: `pbnd` accumulates into a fresh `v_get` vector. | yes |
| `m_resize` | Must keep Meschach's semantics: the top-left block is kept and new rows and columns are zeroed. bnb.C and branch.C resize X and U per node. | yes |
| `m_copy`, `v_copy`, `m_zero`, `m_transp` | Loops or memcpy. | yes |
| `in_prod` | Sequential loop. g++ does not vectorize a floating-point reduction without `-ffast-math`, so the summation order is kept. | yes |
| `Hfactor` + `makeHQ` + `trieig` | Either the port in [`eig_port.C`](meschach_removal/eig_port.C) (167 lines, same operations in the same order) or LAPACK `dsyev`/`dsyevd`/`dsyevr`. | port: yes [M]. LAPACK: no [M]. |
| `m_output`, `error`, `MACHEPS`, `mem_stat_*` | Live code doesn't use them. Delete `svd.C`. | n/a |

**No routine lacks a LAPACK equivalent.** `Hfactor` + `makeHQ` corresponds
to `dsytrd` + `dorgtr`, and `trieig` to `dsteqr`, so `dsyev` does what
`my_svd` does. But no LAPACK routine reproduces Meschach's bits. `Hfactor`
applies the full two-sided Householder update without exploiting symmetry,
`makeHQ` builds Q one column at a time, and `trieig` uses its own Wilkinson
shift and deflation test, while `dsteqr` chooses between QL and QR. At
n ≤ 25, `dsyevd` falls back to `dsteqr`. Its results matched `dsyev`'s
exactly in the benchmark below [M].

**Does `bnbdrv` already link LAPACK?** Yes, by default. `latest/makefile`
has `DNN ?= 1`, and `ldd latest/bnbdrv` shows `libscipy_openblas.so` [M].
Only `make DNN=0` builds without LAPACK. A LAPACK eigensolver would make
LAPACK mandatory and remove the `DNN=0` build. The port adds no dependency.

A side finding: because every `bnbdrv` links OpenBLAS, a QPB-only run starts
OpenBLAS's thread pool, and those threads spin at startup. nug12 reports
0.537 s of CPU against 0.140 s with `OPENBLAS_NUM_THREADS=1`, on the same
246-node tree [M, [log](meschach_removal/logs/openblas_startup.txt)]. The
quick suite showed the same ~0.35 s on every instance. `bench/reference`
does not include this overhead (nug12 0.141 s), but a current `bench.py run`
does, so CPU ratios on sub-second instances are skewed against a build
linked differently. Setting `OPENBLAS_NUM_THREADS=1` in `bench.py`'s child
environment would remove it [E]. This is independent of Meschach.

## 3. Numerical reproducibility

### Experiments

Each variant is a scratch build of `cognac` with only the stated change
([`build_variant.sh`](meschach_removal/build_variant.sh),
[`build_lapack.sh`](meschach_removal/build_lapack.sh),
[`build_port.sh`](meschach_removal/build_port.sh)), checked with
`bench.py check --tree` against `bench/reference/` (exact match required):

| variant | quick suite | standard suite (19) | brute force | source |
|---|---|---|---|---|
| Meschach built by gcc `-O3 -march=native -ffp-contract=off` | exact | exact | – | [M] [std_gccO3](meschach_removal/logs/std_gccO3.txt) |
| `my_svd` on the port (`eig_port.C`) | exact | exact, and the nug20 tree is identical (526,205 nodes, 19,167,276 FW iterations) | 4,536/4,536 | [M] [std_port](meschach_removal/logs/std_port.txt), [bf](meschach_removal/logs/brute_force_port.txt) |
| Meschach built by clang `-O -march=native` (FMA contraction on) | nug12 +7.3%, nug14 +6.3%, nug15 −1.5%, nug16b −1.5%, nug15 root 984.793 | – | – | [M] [quick_fma](meschach_removal/logs/quick_meschach_fma.txt) |
| `my_svd` on LAPACK `dsyev` (OpenBLAS) | all 7 differ | all 19 differ. Node ratio geometric mean 0.997 (min 0.870 tai12a, max 1.114 nug12). Root bounds move by up to 4.4% (scr12) and 0.8% (nug14). | 4,536/4,536 | [M] [std_dsyev](meschach_removal/logs/std_dsyev_ratios.txt) |
| `my_svd` on LAPACK `dsyevr` | all 7 differ, and differently from dsyev (nug14 +4.5%, had16 −1.2%) | – | 4,536/4,536 | [M] [quick_dsyevr](meschach_removal/logs/quick_dsyevr.txt) |

On nug20, QPB only, `dsyev` gave 525,115 nodes and root 2243.88, against
526,205 and 2243.94 [M, [log](meschach_removal/logs/nug20_dsyev.out)].

**The M2 Pro discrepancy.** RUNNING.md says that on an Apple M2 Pro "nug15's
root QPB bound is 984.79 against 986.55 and the quick suite's node counts
differ by up to 8% (nug14 +6%, nug15 -1.5%, nug16b -1.5%)". The FMA build
above reproduces all four figures on x86. The cause: `src/MESCHACH/makefile`
compiles with `clang -O` and no `-ffp-contract=off`. Clang's default
(`on`) fuses multiply-adds wherever the target has FMA, which arm64 always
does. Meanwhile `latest/` is built with `-ffp-contract=off`, so only the
Meschach part of the arithmetic differs. Adding `-ffp-contract=off` to
Meschach's CFLAGS, or using the port, which is built with the solver's
flags, should make the M2 match the Xeon reference on QPB-only runs [E].
It has to be verified on an M2. Apple's libm could still differ in
functions other than `sqrt`.

### Why the search is so sensitive

Meschach and LAPACK agree on the eigenvalues to about 4e-15 relative to
max |λ| [M, [eigbench](meschach_removal/logs/eigbench_nug20.txt)], yet root
bounds move by up to 4%. There are two reasons [E, from the code]:

- The eigenvectors are not unique where eigenvalues repeat. Grid distance
  matrices (nug, had) and sparse flows (scr, chr) have repeated eigenvalues,
  and QPB builds S = V·W·diag(s)·W'·V' from the eigenvectors. That depends
  on the basis chosen inside an eigenspace whenever `find_basis` assigns
  different s values within it. So the QP bound is a function of the
  eigensolver, not of A and B alone. Every valid choice gives a valid
  bound, which is why brute force passes with LAPACK.
- The bound is the Frank-Wolfe value after a fixed iteration budget, and
  each FW step solves a LAP whose optimum can jump on a tiny perturbation.

So a different eigensolver changes the bound at most nodes by much more
than rounding, and node counts move on every instance, not only near the
fathoming threshold.

### What would change if the eigenvalues came from LAPACK

- `bench/reference/` (all 19 instances) and the shared baseline
  `bench/results/ref-cognac-b` would need `bench.py check --update`.
- Every QPB-only node count quoted in RUNNING.md: nug20 526,205, nug24
  ~25M, nug25 52,989,998 and 40,183,792 (nug_deep).
- The DNN-run counts in RUNNING.md and `experiments/nug30_solve` (nug20
  `-x 2:2` 84,592 and `-x 2:3` 141, nug24 299, nug25 717/761, nug27 675,
  nug28 1,736, nug30 3,820/3,844) depend on QPB through the gap rule
  `(inc - QPB)/(inc - root QPB)`, the lookahead and the nodes QPB fathoms
  before the DNN runs. They would likely move, but less [E]. RUNNING.md
  records that the M2 Pro reproduced nug27 and nug28 `-x 2:5` exactly,
  even though its QPB values differ from the Xeon's (the FMA experiment
  explains why).
- Unaffected: `bench/LOG.csv` rows (history, and `--rev` baselines rebuild
  the old code with its own Meschach). Also unaffected are the pass criteria
  of `make bruteforce`, `resume`, `parallel`, `coord` and `shipcheck`, which
  compare against enumeration or against the same binary, and `make
  dnncheck`, because `dnntest` never calls the QP bound (no `qpbnd`,
  `init_ST` or `my_svd` reference).
- Cross-machine reproducibility would get worse. OpenBLAS picks CPU-specific
  kernels (with FMA on Haswell and later) inside `dsytrd`, so QPB node counts
  would become machine-dependent, as the DNN's already are [E]. Accelerate
  on macOS would differ again.

### Staging and bit-identity

| stage | same bits on x86? | evidence |
|---|---|---|
| 0. Meschach CFLAGS (`-ffp-contract=off -std=gnu17`), compiler-neutral `CC` | yes | gcc `-O3 -march=native -ffp-contract=off` matches the standard suite [M] |
| 1. Eigen port | yes | 656/656 matrices bitwise, standard suite, nug20 tree, brute force [M] |
| 2. Own `MAT`/`VEC` and allocation | yes (copies and zero fills only) | by construction [E]; acceptance below |
| 3. Remove `src/MESCHACH` from builds | yes | no arithmetic involved |
| LAPACK eigensolver (not proposed) | **no**, under any LAPACK | table above [M] |

On arm64, stages 0 and 1 are intentionally not bit-identical to today's
M2 results. They move the M2 onto the Xeon's numbers, which also changes
earlier laptop measurements such as the M2 rows in RUNNING.md and the
laptop nug30 runs. Bit-identity is impossible only if the eigen algorithm
or its operation order changes.

**Validating a non-identical stage** (if the LAPACK route is ever taken, or
the port is later restructured, for example to exploit symmetry):

1. `tests/brute_force.py` with `--seed 1`, `2` and `3`, and with
   `--solver-args "-x 0:20"`, all passing.
2. `make dnncheck`, `resume`, `parallel` and `shipcheck` unchanged.
3. Standard suite (`bench.py run`, both binaries back to back on the same
   machine): geometric-mean node ratio within 0.98 to 1.02, and no instance
   outside 0.85 to 1.15. The `dsyev` variant measured 0.997 with
   0.870 to 1.114 [M], so per-instance swings of ±13% are normal and not a
   sign of a bug.
4. `large` suite, plus nug24 with `--jobs 2` (or `bench/paired_sample.py`
   extended to take two binaries, about half a day [E]) for a CPU ratio
   with a confidence interval. Single nug20 runs varied by 5% between
   repeats here [M: 125.8, 120.0 and 121.1 s].
5. Bound quality on fixed nodes: dump a depth-3 frontier with `-d 3 -o`,
   compute each node's QPB with both binaries, and require a mean difference
   indistinguishable from zero (no systematic weakening).
6. nug25 `nug_deep.par -x 2:5` within a few percent of 717 nodes, then
   `bench.py check --update` and a LOG.csv row explaining the move.

## 4. Performance

### Method

[`build_timer.sh`](meschach_removal/build_timer.sh) builds a scratch
`bnbdrv` linked with `-Wl,--wrap` for each of the 17 functions.
[`mes_wrap.C`](meschach_removal/mes_wrap.C) reads the TSC around the
outermost Meschach call and times `my_svd` as a whole, with a size histogram.
With `MES_DUMP` it also writes every 2000th `my_svd` input to a file for the
benchmarks below. A gprof build (`-pg` on the solver and Meschach) was run
as a cross-check. perf is not installed in the container. All runs used
`OPENBLAS_NUM_THREADS=1`, and node counts equal the uninstrumented build's.

### Share of a run [M]

| run | total CPU | all Meschach | eigen (Hfactor+makeHQ+trieig) | alloc, copy, permute | log |
|---|---|---|---|---|---|
| nug20 QPB only, default.par, `-u 2570`, 526,205 nodes | 127.5 s (uninstrumented 120.0 to 125.8 s) | 26.2 s, 20.6% | 17.0 s, 13.3% (my_svd with its own allocations: 18.8 s, 14.8%) | 9.2 s, 7.3% (`m_get` 3.6 s, `m_free` 1.8 s, `m_copy` 0.8 s) | [timer](meschach_removal/logs/nug20_timer.err) |
| nug20 `-x 2:2` (84,592 nodes, DNN on 97) | 53.4 s (DNN 41.7 s) | 2.73 s, 5.1% | 1.83 s, 3.4% | 0.9 s | [timer](meschach_removal/logs/nug20_x22_timer.err) |
| nug20 `-x 2:3` (141 nodes, DNN on 133) | 48.6 s (DNN 47.9 s) | 0.035 s, 0.07% | 0.015 s | 0.02 s | [timer](meschach_removal/logs/nug20_x23_timer.err) |
| nug30 record setting | 3,270 s, DNN 3,248 s (RUNNING.md) | under 0.7% [E: all non-DNN time is 22 s] | | | |

The wrappers add about 15 to 20 ns to each of the ~90M timed calls on
nug20, so the allocation row is overstated by up to ~1.8 s [E]. The
corrected Meschach share is 19 to 21%. The eigen row has only 3.9M timed
calls and is accurate.

gprof on the same nug20 run put `dlap` (the LAP solver) at 53.8% and `fw`
at 18.1% of the sampled time, and Meschach's own code at 9.0% [M,
[log](meschach_removal/logs/nug20_gprof_share.txt)]. gprof undercounts
Meschach. Its eigen code makes over a billion tiny internal calls
(`__mltadd__` 581M, `__ip__` 404M, `rot_cols` 196M, `get_col` 209M). Under
`-pg` each pays for an `mcount` call, and that time, like libc's `malloc`
and `free`, lands outside the sampled text: 36 of the run's 147 s were
unsampled. The wrapper numbers are the ones to use.

### Eigensolvers on the actual matrices [M]

[`eigbench.C`](meschach_removal/eigbench.C) runs 656 `my_svd` inputs
sampled from the nug20 run, 50 passes each, single thread, and reports
microseconds per call
([log](meschach_removal/logs/eigbench_nug20.txt)).
[`eigport_test.C`](meschach_removal/eigport_test.C) times the port
([log](meschach_removal/logs/eigport_nug20.txt)). The port column comes from
eigport_test's own run of the same matrices, whose Meschach timings were
within 20% of eigbench's. Its port/Meschach time ratio over all matrices was
0.648, and 0.658 on a rerun:

| n | calls in run | Meschach | dsyev | dsyevd | dsyevr | port |
|---|---|---|---|---|---|---|
| 10 | 68,666 | 8.6 | 8.3 | 9.9 | 15.8 | 5.8 |
| 11 | 199,618 | 10.1 | 10.9 | 10.8 | 14.7 | 6.3 |
| 12 | 334,948 | 11.7 | 11.8 | 11.8 | 17.0 | 7.5 |
| 13 | 296,794 | 13.5 | 13.1 | 13.2 | 18.4 | 8.9 |
| 14 | 315,408 | 15.7 | 14.3 | 14.9 | 20.1 | 9.8 |
| 15 | 68,756 | 19.9 | 17.0 | 16.5 | 20.1 | 11.5 |
| 16 | 6,360 | 23.1 | 18.5 | 17.5 | 16.6 | 12.7 |
| sample mean | | 13.3 | 12.8 | 13.0 | 17.8 | 0.65x Meschach |

LAPACK includes the copy into its column-major buffer, and its workspace is
queried once per size. At these sizes it gains 4% over Meschach (dsyev) or
loses 34% (dsyevr). The port gains 35% while keeping every bit. Its speedup
comes from reading the Householder vectors once instead of once per column
of Q, rotating rows of Qᵀ in contiguous memory instead of strided columns,
dropping the final transpose, and not allocating per call.

### End to end [M]

nug20 QPB only, baseline and port run at the same time on separate cores,
twice ([logs](meschach_removal/logs/)):

| pair | baseline | port | ratio | tree |
|---|---|---|---|---|
| 1 | 120.03 s | 112.75 s | 0.939 | identical (526,205 nodes, 19,167,276 FW iterations, root 2243.94) |
| 2 | 121.06 s | 113.69 s | 0.939 | identical |

That is consistent with 13.3% eigen time × 35% saved, plus the
allocations `my_svd` no longer makes [E]. The remaining allocation churn
(~35 `m_get` and ~46 `m_free` per node, ~5 s or 4% of nug20) could be cut
by reusing buffers once the solver owns its allocator. That is a separate,
bit-identical optimization, with gain estimated at 3 to 4% [E].

## 5. Plan and cost

Effort is in engineering days for someone who knows the code [E].

| stage | work | days | acceptance test |
|---|---|---|---|
| 0. Build hygiene | `src/MESCHACH/makefile`: `CC ?= cc`, `CFLAGS = -O -std=gnu17 -ffp-contract=off`. `bench.py build()` honors `CC`/`MESCHACH_CC`. Drop the clang install from `infra/phase0/bench.sh` and the clang requirement from RUNNING.md. | 0.5 | `make test` exact with gcc and with clang on x86 (gcc `-O3` variant measured exact). On an M2, `bench.py check` (quick) expected to match the reference, then delete the M2 caveat in RUNNING.md. |
| 1. Eigen port | Move `docs/meschach_removal/eig_port.C` into `latest/` and have `my_svd` call it. Drop `Hfactor`, `makeHQ`, `trieig` and `m_transp`. Rerun eigport_test on dumps from a few families (nug, scr, chr, tai, brute-force sizes). | 1–1.5 | `bench.py check --suite standard` exact. Brute force with seeds 1–3. `resume`, `parallel`, `coord`. nug20 tree identical. Paired standard-suite and nug20 timing logged to LOG.csv (expect nodes 1.000 and nug20 time ≈ 0.94). |
| 2. Data structures | `latest/mat.h/.C`: `MAT`/`VEC` with the same field names and `m_get`, `m_free`, `m_resize` (same top-left semantics), `m_copy`, `m_zero`, `v_*`, `in_prod`, `M_FREE`, `V_FREE`. Replace `PERM` with `int *`. Delete `svd.C`. Remove `#include "meschach.h"`. A rename to a modern API (a `Matrix` class, std::vector) touches ~1,500 lines and would take 4–6 days, so do it later if at all. | 2–3 | All gates exact (`make test`, `dnncheck`, `shipcheck`, `continuationcheck`). Quick suite under valgrind or ASan with no new leaks (the `init_w0` leak fix and `m_resize` semantics are the risk). |
| 3. Remove from builds | `Makefile` and `latest/makefile` (`INCLUDES`, `LDLIBS`). `bench.py` builds Meschach only if `src/MESCHACH` exists, so `--rev` of older commits still works. Update `infra/phase0/bench.sh`, RUNNING.md and CONTRIBUTING.md, then delete `src/MESCHACH/` (107 tracked files). Experiment build scripts stay pinned to the commits they document. | 0.5 | A fresh clone builds with `make`, `g++` and (for `-x`) LAPACK only. `make test`. `bench.py run --rev <pre-removal commit>` still builds. |
| (4. Allocation reuse) | Per-depth work buffers instead of `m_get`/`m_free` per node. | 1–2 | Exact gates; ~3–4% QPB-only time [E]. |
| (not proposed) LAPACK eigen | `dsyev` in `my_svd`, plus a reference update | 1 + 2–3 validation | Section 3 plan. Gains 0–4% eigen time, changes every node count, machine-dependent. |

Stages 0–3 total 4 to 5.5 days [E].

**Risks.**

- `m_resize` and zero-fill semantics in stage 2. Callers rely on calloc'd
  storage, and a wrong `m_resize` corrupts X and U silently. Brute force
  catches invalid bounds, while the exact regression check catches
  everything else.
- The port must keep the operation order. Any "cleanup" of `eig_port.C`
  (exploiting symmetry, fusing loops, reordering sums) breaks bit-identity.
  Keep a bitwise unit test (eigport_test on a committed dump) in `make test`.
- On arm64 the results change by design (stage 0/1). Earlier laptop numbers
  stop being reproducible at HEAD.
- Meschach's error paths (`error()` → longjmp/exit) disappear, so a negative
  size would no longer abort with a message. Add asserts.
- Doing nothing has its own risks. Meschach is 1994 K&R C. It fails to
  compile under C23 with gcc 13 and clang 18 [M,
  [log](meschach_removal/logs/meschach_c23_build.txt)], and GCC 15 defaults
  to C23. It needed the `isatty` fix to build at all. Its static work vectors
  (`Hfactor`, `makeHQ`, `hhtrcols`) make the QP bound non-reentrant, which
  matters only if the solver ever goes multithreaded in one process.

**Recommendation.** Do stage 0 now. It is half a day, removes the clang
requirement properly (the AWS kit's fallback becomes the default), protects
against the C23 default, and very likely makes M2 QPB results match the
Xeon reference. Do stage 1 when convenient. It is the only stage with a
speed payoff (6% on QPB-only runs [M]) and changes no results. The
prototype and its bitwise test exist. Stages 2–3 are hygiene: they delete a
vendored library and the last C toolchain quirk, but buy no speed, so take
them when someone touches the bound code anyway. The current work (nug30,
wil50) spends under 1% of its time in Meschach. Do not replace the QP
bound's eigensolver with LAPACK: it is no faster at n ≤ 30, changes every
node count, and ties QPB results to the BLAS build.

## Files

- [`symbols_by_object.sh`](meschach_removal/symbols_by_object.sh) and
  [`inventory.py`](meschach_removal/inventory.py): the inventory.
- [`build_timer.sh`](meschach_removal/build_timer.sh) and
  [`mes_wrap.C`](meschach_removal/mes_wrap.C): the timing build, including
  the matrix dump.
- [`gprof_share.py`](meschach_removal/gprof_share.py): the gprof summary.
- [`eigbench.C`](meschach_removal/eigbench.C): Meschach against LAPACK on
  dumped matrices.
- [`eig_port.C`](meschach_removal/eig_port.C) and
  [`eig_port.h`](meschach_removal/eig_port.h): the port.
  [`eigport_test.C`](meschach_removal/eigport_test.C) is its bitwise and
  timing test.
- [`build_variant.sh`](meschach_removal/build_variant.sh),
  [`build_lapack.sh`](meschach_removal/build_lapack.sh) (with
  [`my_svd_lapack.C`](meschach_removal/my_svd_lapack.C)) and
  [`build_port.sh`](meschach_removal/build_port.sh): the scratch builds.
  Each copies `latest/` and `src/MESCHACH/` to a destination directory and
  leaves the checkout alone.
- [`ratio.py`](meschach_removal/ratio.py): node ratios against
  `bench/reference`.
- [`logs/`](meschach_removal/logs/): every run quoted above.
