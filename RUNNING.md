# Building, running and benchmarking

How to build the branch-and-bound solver in `latest/`, run it on QAPLIB
instances, and compare performance between code or parameter changes.

## Build

Requires `make`, `g++` and a C compiler (`cc`: gcc or clang; `CC=clang make`
picks one). For the DNN bound (`-x`), LAPACK too (see `latest/makefile`).

```sh
make          # from the repo root: builds src/MESCHACH/meschach.a, then latest/bnbdrv
make test     # brute-force tests + regression check (see "Correctness checks")
```

Or build the two parts by hand:

```sh
cd src/MESCHACH && make basic     # builds meschach.a (old-style C warnings are expected)
cd ../../latest && make           # builds latest/bnbdrv
```

After changing solver code, run `make` again. It rebuilds whatever depends on
the files you changed, headers included (the objects' `.d` files list them).
`bench.py` refuses to run a binary that is older than the sources in `latest/`,
since it records results under the current commit.

`latest/makefile` compiles with `-O3 -march=native -ffp-contract=off`, which
is about 2% faster than plain `-O3`. Don't drop `-ffp-contract=off`. Without
it, the compiler fuses multiply-adds where the CPU supports them, which
changes rounding, and then node counts depend on the machine.

## Run one instance

Always run from `latest/`. The solver looks for its parameter and symmetry
files relative to the current directory.

```sh
cd latest
./bnbdrv ../qaplib/nug12.dat                        # instance's own defaults
./bnbdrv -p PARAM/default.par ../qaplib/nug12.dat   # a specific params file
./bnbdrv -p PARAM/default.par -u 578 ../qaplib/nug12.dat   # start from a known upper bound
```

With no flags, the solver uses:

| What | Default | Notes |
|---|---|---|
| Params file | `./PARAM/<name>.par` | If this file doesn't exist, the solver stops with an error. The reader also rejects a file whose level count or layout is off (it used to run on whatever values ended up in each field). Only a few instances (nug*, kra30b, esc32d, rou15, scr15, tai17a, tai20a, ste36a, …) have their own params file, so pass `-p` for anything else. |
| Symmetry file | `./SYM/<name>.sym` | Optional. If missing, the solver prints "assuming no symmetry" and continues. |
| Initial upper bound | 100000 | Effectively none. The solver finds incumbents itself and prints `New objective value: …` each time it improves. |

Other flags (`./bnbdrv` with no arguments prints help):

| Flag | Meaning |
|---|---|
| `-p F` | params file |
| `-s F` | symmetry file |
| `-u X` | initial upper bound (incumbent) |
| `-d X` | stop at tree depth X (overrides the params file) |
| `-f X` | scale reported times by X |
| `-h` | run the heuristic only (no branch and bound) |
| `-e` | estimate tree size (reads `estimate.par`) |
| `-b` | benchmark mode |
| `-c` / `-a F` | interactive console / run console commands from file F |
| `-l F` | log file |
| `-o F` | write the search state to F: at checkpoints, on `-t`, and at the end (with `-d D`, the open nodes at depth D) |
| `-k X` | with `-o`, checkpoint every X CPU seconds |
| `-t X` | with `-o`, stop after X CPU seconds, saving the state; exit code 3 |
| `-r F` | resume from state file F instead of the root |
| `-x D1:D2[:G[:H[:W]]]` | DNN bound (latest/dnn.C) on nodes QPB does not fathom at depths D1-D2 with relative gap (inc - QPB)/(inc - root QPB) at most G (default 1); H: give up on a node after H iterations when hopeless (default 0, never); W=1: start a node's ADMM from its parent's final multiplier when the parent's DNN bound missed (default 0; on nug25 -x 2:5:1:300 it cut iterations per fathomed node ~40% but DNN CPU only 4%). Needs LAPACK (see below) |
| `-w W` | Frank-Wolfe patience at lookahead/strong levels (default 0 = off; see `bench/LOG.csv`) |

### Long runs: checkpoint and resume

```sh
./bnbdrv -p PARAM/default.par -u 3488 -o nug24.state -k 600 ../qaplib/nug24.dat   # checkpoint every 10 min
./bnbdrv -p PARAM/default.par -o nug24.state -k 600 -r nug24.state ../qaplib/nug24.dat  # after a crash
```

The state file (format in `latest/state.h`) holds the open nodes, the
incumbent, the root bound and all statistics, plus the Frank-Wolfe point of
each open node's parent (a child starts its bound there), with doubles at full
precision,
so a resumed run makes exactly the same decisions: the final node counts, FW
iterations and per-level tables equal an uninterrupted run
(`tests/resume_test.py`). A SIGUSR1 to a run started with `-o` also saves the
state and stops it at the next node.

With the DNN bound (`-x`), a branched node's final ADMM state (Z and Y), from
which its children start, is held in memory only. It goes into the state file
(version 6) only with `QAP_DNN_SHIP=16|32|64` (bf16, fp32 or the exact
doubles). With 64, a stopped and resumed DNN search is the uninterrupted
search node for node (`make shipcheck`). Without it, nodes restored from a
file start their ADMM cold, which costs iterations but nothing else. This
matters for `bench/parallel.py`, whose tasks are state files: on nug30,
`QAP_DNN_SHIP=32` took 15% off the CPU on an identical tree.

### Several cores: bench/parallel.py

```sh
python3 bench/parallel.py qaplib/nug20.dat --jobs 4 --ub 2570
python3 bench/parallel.py qaplib/nug24.dat --jobs 4 --ub 3488 --checkpoint nug24.ckpt --checkpoint-every 600
python3 bench/parallel.py qaplib/nug24.dat --jobs 4 --ub 3488 --resume nug24.ckpt   # after a restart
```

For scale: nug25 with the known optimum as the incumbent (`--ub 3744
--jobs 4`, default.par at db8df94, 4 cores at 2.8 GHz) took 52,989,998 nodes,
18,353 CPU seconds and 78 minutes wall. `bench/sample_estimate.py` had
predicted 18,652 +- 11,627 CPU seconds from 600 sampled depth-3 subtrees. With
`--params PARAM/nug_deep.par` (rule 6 down to depth 6 and fewer "can't fathom"
FW iterations; tuned for grid instances, slower on esc16a and tai17a, see
bench/LOG.csv) the same run took 40,183,792 nodes, 12,362 CPU seconds and
53 minutes.

The root is solved serially and its children go into a pool; up to `--jobs`
solver processes each take a chunk of nodes. When a worker is idle and the
pool is empty, the busiest worker is sent SIGUSR1 and its remaining nodes are
split among the idle ones. Statistics are summed and printed in the solver's
own format (`Total Time` is CPU summed over workers). With a fixed upper bound
such as the optimum, node and FW totals are identical to a serial run
(`tests/parallel_test.py`). `--checkpoint F` periodically stops all workers and
merges every open node into F, an ordinary state file: `--resume F` continues
it, or `bnbdrv -r F` finishes it serially. `bench.py run --jobs N` runs each
instance this way; compare such runs with `compare --time wall`.

On this 4-vCPU container nug16a ran 3.2x faster in wall time with 4 workers and
nug17 2.7x; total CPU rises somewhat under full load, so expect less than
N-fold on shared machines.

Known optimal solutions are in `QAPSOLN/<name>.sln`. The first line is `n` and
the optimal objective, and the second line is the permutation. A few files are
broken: kra30a, kra30b, kra32, esc32d, tai25a and tho30 contain `<value> ??`,
and nug27's claims an optimum of 5. The tools treat those optima as unknown.

When the solver starts, it prints two facts about the instance that affect
the search:

- `objective values are even integers: fathoming tolerance 1.95`. A node is
  fathomed when its bound exceeds `incumbent - tolerance`. The tolerance is
  1.95 only when every objective value is provably even (A and B symmetric
  integer matrices, even a_ii·b_kk, even C), 0.95 for other integer data, and
  0 otherwise. This used to be 1.95 for every instance, which could cut off
  optima with odd values.
- `symmetrized B`. The bound assumes symmetric A and B. When only one of them
  is symmetric, the other is replaced by (M + M')/2, which leaves every
  objective value unchanged. If neither is symmetric (for example the bur
  family), the solver prints a warning and its results aren't reliable.

### Output

The run ends with a `---- RESULTS ----` block:

```
    Total Time: 0.243727 seconds        <- CPU time (clock()), not wall time
         Nodes: 2455 / 1635 / 65        <- nodes evaluated / fathomed / enumerated
    root bound: 477.493
 fw iterations: 132715                  <- total Frank-Wolfe iterations
```

After that come a per-depth table (nodes, fraction fathomed, fraction of
columns eliminated, time), a breakdown of which params level handled each
depth, and relative-gap statistics.

### Params file format

Each value goes on its own line, and anything after the value is a comment.
Example from `PARAM/default.par`:

```
0     has the location of a symmetry file been provided?   (if 1, next line is the path)
-1    restriction on tree depth (-1 means none)
3     number of branching levels (L)
--- repeated L times; the last level omits the first two lines ---
0.20  use this level when the relative gap (opt-bnd)/opt is worse than...
2     ...and depth is at most this
4     branching rule: 1 U-sum, 2 U-fathom, 3 strong, 4 lookahead, 5 learned (BUFF), 6 eccentric, 7 eccentric lookahead, 8 DNN marginals, 9 lookahead between the two marginal lines (8 and 9 need -x; lookahead where no DNN bound ran)
150   Frank-Wolfe iterations for the bound
100   Frank-Wolfe iterations when we know we can't fathom
25    Frank-Wolfe iterations for row selection (lookahead children; each starts from the parent's X)
1     bound: 0 QP, 1 QP parametric, 2 QP+GLB, 3 QP improved, 4 GLB, 5 PB, 6 EVB3, 7 BUFF
30    step for parametric bound
30    compute B for the best _ rows of U
```

Rule 5 (`BRANCH_BUFF`) reads coefficients from `latest/lgb_init.txt`, which
isn't in the repo. Don't use rule 5 without that file.

Rule 6 branches on the free location farthest from the other free locations
(largest B entry), at no per-node cost. `PARAM/eccentric.par` is
`default.par` with rule 6 instead of lookahead at depths 3-4. It is much
faster on dense instances (nug20 2.1x, nug21 1.7x, tai17a 2.0x, standard
suite 1.3x) but slower where the flow matrix is sparse and lookahead pays
for itself (chr15a 1.4x slower, esc16a more than 2.9x slower), so it is not
the default. See `bench/LOG.csv`.

Rule 7 is lookahead on a smaller candidate set: the `use_best` most eccentric
free locations (columns) and the `use_best`-1 rows with the largest U sums,
instead of rule 4's `use_best` best rows and columns by U. `default.par`
uses it at depths 3-4 with `use_best` 3. Measured against rule 4 there
(5 rows + 5 columns): standard suite time 0.938, nug24 0.937, large set
about 0.84 (esc16a 0.69, chr15a 1.13 slower). It keeps lookahead's payoff
on sparse-flow instances, which rule 6 loses.

Rules 8 and 9 need the DNN bound (`-x`) and fall back to lookahead where it
did not run.
- **Rule 8** branches on the line (row or column) of the DNN solution's
  marginals x̄ with the largest sum of squares: the line where the relaxation
  is most decided. On nug30 it took 16% off the CPU (4,904 to 4,133 s, depth-3
  misses 85 to 59). `PARAM/nug_deep_xbar_all.par` uses it at every node
  through depth 6; with the DNN bound from depth 1 it is the nug30 record
  setting.
- **On nug28 rule 8 lost 37%.** Rule 9 is a hedge: lookahead between the two
  marginal lines (rule 8's row and its column), chosen by the children's QPB
  bounds. It gives nug30 -13% and nug28 +8%
  (`experiments/wil50_feasibility/tuning`).

### The DNN bound (-x)

`latest/dnn.C` computes a certified lower bound from the doubly nonnegative
(DNN) SDP relaxation (experiments/dnn_bound). It is far stronger than QPB and
far more expensive, so it pays only near the top of the tree: nug20 with
`-x 2:2` takes 84,592 nodes and 117 s (QPB alone 526,205 nodes, 197 s), with
`-x 2:3` 141 nodes and 105 s. With 4 workers (bench/parallel.py --extra "-x ..."),
the known optimum as incumbent, 2.8 GHz:

| run | nodes | CPU s | wall |
|---|---|---|---|
| nug24 default.par, QPB only | ~25M | ~5,000 | ~25 min |
| nug24 default.par -x 2:3 | 299 | 796 | 4.4 min |
| nug25 default.par, QPB only | 52,989,998 | 18,353 | 78 min |
| nug25 nug_deep.par -x 2:3 | 10,234,599 | 5,065 | 22 min |
| nug25 nug_deep.par -x 2:4 | 4,622,159 | 4,055 | 18 min |
| nug25 nug_deep.par -x 2:5 | 717 (DNN: 710, 678 fathomed, 3,256 s) | 3,279 | 15 min |
| nug27 nug_deep.par -x 2:5 | 675 (DNN: 663, 641 fathomed, 4,405 s) | 4,428 | 20 min |
| nug28 nug_deep.par -x 2:5 | 1,736 (DNN: 1,727, 1,664 fathomed, 12,926 s) | 12,976 | 62 min |
| nug30 nug_deep.par -x 2:5:1:300, then 2:6:1:300 | 45,486,332 (DNN: 5,615, 5,403 fathomed, 26,770 s) | 40,555 | ~3.0 h (experiments/nug30_solve) |
| nug30 nug_deep.par -x 2:6:1:300 (clean run) | 6,532,419 (DNN: 5,764, 5,551 fathomed, 17,192 s) | 18,614 | ~1.3 h (experiments/nug30_solve) |
| nug30 nug_deep.par -x 2:6:1:300 (A/B, ff1f044, 2.8 GHz container) | 6,532,419 (DNN: 5,764, 5,551 fathomed, 26,726 s) | 28,428 | ~2.1 h |
| nug30 nug_deep.par -x 2:7:1:300:1 (same machine) | 5,817 (DNN: 5,807, 5,595 fathomed, 23,537 s) | 23,618 | 1.7 h (experiments/nug30_solve) |
| nug30 nug_deep.par -x 2:7:1:300:1 + Peaceman-Rachford, degree 4 (Nov 1, 2.1 GHz) | 4,528 (DNN: 4,518, 4,355 fathomed, 11,677 s) | 11,748 | 50.5 min (experiments/wil50_feasibility/stronger_bounds/search) |
| **nug30 record setting (Nov 5, 2.1 GHz, below)** | **3,820 (DNN: 3,819, 3,678 fathomed, 3,248 s)** | **3,270** | **14.2 min** (experiments/nug30_solve/nov5_cloud) |
| nug25 nug_deep.par -x 2:5:1:300 (at 0a179fc, 2.8 GHz) | 761 (DNN: 754, 720 fathomed, 1,086 s) | 1,100 | 4.7 min |
| nug25 nug_deep.par -x 2:5:1:300:1 (warm start) | 717 (DNN: 710, 678 fathomed, 1,039 s) | 1,055 | 4.5 min |

The same runs on an Apple M2 Pro laptop (6 P + 4 E cores, Accelerate,
`--jobs 8`, 2026-10-29). CPU is summed over 8 workers, so compare wall:

| run (M2 Pro, 8 workers) | nodes | CPU s | wall |
|---|---|---|---|
| nug25 nug_deep.par -x 2:5:1:300:1 (4 workers, OpenBLAS) | 750 (DNN: 743, 709 fathomed, 554 s) | 684 | 2.6 min |
| nug27 nug_deep.par -x 2:5 | 676 (DNN: 663, 641 fathomed, 1,357 s) | 1,363 | 3.3 min |
| nug28 nug_deep.par -x 2:5 | 1,736 (DNN: 1,727, 1,664 fathomed, 3,994 s) | 4,003 | 9.2 min |
| nug30 nug_deep.par -x 2:7:1:300:1 | 585,550 (DNN: 5,834, 5,619 fathomed, 8,240 s) | 8,355 | 19.7 min (experiments/nug30_solve) |
| nug30 nug_deep.par -x 2:7:1:300:1 (Apple M1, 4 P + 4 E cores, 8 workers) | 536,075 (DNN: 5,858, 5,643 fathomed, 18,231 s) | 18,383 | 42.9 min (experiments/nug30_solve) |
| **nug30 record setting (Nov 4, below)** | **3,844 (DNN: 3,843, 3,701 fathomed, 1,826 s)** | **1,838** | **5.4 min** (experiments/nug30_solve) |

With plain `-x 2:5` the trees match the Xeon runs exactly (nug27, nug28).
With the give-up threshold and warm starts (`-x ...:300:1`) a few decisions
at the cutoff go the other way: nug25 took 750 nodes instead of 717, and on
nug30 two depth-7 nodes missed where the Xeon fathomed all 96, so 580,000
QPB-only nodes followed below depth 7 (84 CPU s, 1% of the run). Per DNN
node the laptop is about 3x faster than the 2.8 GHz Xeon (depth 3: 1.9 s
against 5.8; depth 4: 1.2 against 3.5).

**The nug30 record setting** (Nov 4: 1,838 CPU s and 5 min 22 s on the
M2 Pro with 8 workers; Nov 5: 3,270 CPU s and 14 min 9 s on 4 cores of a
2.1 GHz Xeon; how each switch got there is in
`experiments/nug30_solve/README.md`):

```sh
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
       QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 \
       QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32
python3 bench/parallel.py qaplib/nug30.dat --jobs 8 --ub 6124 \
    --params PARAM/nug_deep_xbar_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint nug30.ckpt --checkpoint-every 300
```

It is tuned on nug30. Rule 8 lost 37% on nug28, so for other sizes
compare against `nug_deep.par` or rule 9 first.

It needs LAPACK. On macOS the makefile links Apple's Accelerate framework
(on the M2 Pro, 30 ms per ADMM iteration on a nug30 m=28 node against 40 ms
with scipy-openblas32 and 73 ms on the Xeon; Accelerate's per-process speed
drops about 1.5x with 10 processes running, since the AMX units are shared,
but total throughput still rises up to 10). The solver cannot set
Accelerate's thread count, so `export VECLIB_MAXIMUM_THREADS=1` before
running several solver processes by hand; `bench/parallel.py` sets it for
its workers, and `bnbdrv -x` warns when it is unset. Elsewhere, `pip
install scipy-openblas32` gives an OpenBLAS the makefile finds (it falls
back to the system liblapack/libblas, ~7x slower); `make DNN=0` builds
without it (then `-x` is an error). Rebuild after installing or removing
the OpenBLAS wheel. Its floating-point results depend on the BLAS kernels
chosen for the CPU, so node counts with `-x` can differ between machines
(rarely: only a bound within ~1e-6 of the cutoff decides differently).
`make dnncheck` validates the bound itself; `python3 tests/brute_force.py
--solver-args "-x 0:20"` checks it inside the search.

Opt-in DNN settings, as environment variables read by `bnbdrv -x` (all off
by default; each has its experiment report): `QAP_DNN_GAMMA1=1.0
QAP_DNN_GAMMA=1.0` the Peaceman-Rachford step and `QAP_DNN_CHEB=4` the
filter degree (experiments/eigensolver: 1.5x less DNN CPU on nug30);
`QAP_DNN_TIGHT_LAM=2` with `QAP_DNN_CHEB=3` (degree-4 strength at 85% of the
cost); `QAP_DNN_GAP_RATIO=0.65` the gap-ratio give-up rule
(experiments/wil50_feasibility/tuning); `QAP_DNN_EXT=20` keeps the DNN bound
up to 20 levels below a last-level miss instead of handing its subtree to
QPB alone, and `QAP_DNN_GU_ZONE=0.003` makes a give-up within 0.3% of the
cutoff wait 100 iterations for confirmation (both from J.-P. Goux's nug30
solver; experiments/wil50_feasibility/tuning/goux). `QAP_DNN_CUT_ZONE=0.002`
adds a cut phase after a DNN miss that ends within that fraction of the
cutoff: regional triangle cuts (latest/dnncut.h) from the node's own iterate,
`QAP_DNN_CUT_GRID=R:C` the grid whose row and column prefixes are the
regions (wil50 `5:10`; nug30 `5:6` with `QAP_DNN_CUT_SIDE=1`, the second
matrix), `QAP_DNN_CUT_ROUNDS` (12), `_ITERS` (300), `_K` (1000), `_DELTA`
(50), `_K5` (0: also this many of the regional 5-point count inequalities a
round, latest/dnncut.h); every value stays certified
(experiments/wil50_feasibility/stronger_bounds). On nug30 the cut phase
costs more than it saves (+34% CPU on the laptop, +42% on the cloud), so the
record setting leaves it off.

Later switches, all used in the nug30 record setting except the last two
items (defaults in parentheses; experiment reports in
`experiments/wil50_feasibility/`):
- `QAP_DNN_BETA=F`: the ADMM penalty is F·m/3 (0.3).
- `QAP_DNN_REFRESH=N`: a full eigendecomposition every N iterations (20).
  Goux's 0.21 and 120 together took 20% off nug30 (`tuning/goux`); on wil50
  they do not pay.
- `QAP_DNN_WARM_Y=2`: a child starts from its parent's final Y, conditioned
  on its new assignment, as well as its Z (0: Z only; `warm_y/`). With it, a
  child's `-x ...:W` warm start carries the primal too.
- `DNN_CHOLQR_PASSES=1`: one CholQR pass in the subspace step (2;
  `roadmap/node_speed.md`).
- `QAP_DNN_SHIP=32`: carry the parent's ADMM state in state files, as
  above (0).
- `DNN_AMX=1`: Intel AMX bf16 kernels for the subspace products, on CPUs
  with AMX-BF16 only (`amx/`).
- The wil50 cut-phase switches `QAP_DNN_POOL`, `QAP_DNN_CUT_WARMY` and
  `QAP_DNN_CUT_DAMP1` (`pool_search/`, `pipeline_ab/`). Unused on nug30.
- `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7` (unset): a schedule for the ADMM
  penalty. Cold calls start at 1.6667 x the base penalty and drop to the base
  at the first full eigendecomposition at or after iteration 60. Warm-started
  calls run at 0.7 x the base throughout (`speed/beta_sched/`). On wil50: 13%
  less CPU on 32 cold nodes, and on the node-10 subtree with the production
  candidate 1.14x fewer DNN iterations on 232 matched nodes, 1.25x end to end
  (`pipeline_ab/`, arm `fast`). Not yet tried in a nug30 search, whose base
  penalty is 0.21.
- `DNN_R32_NONCERT=1` (unset): R's rank term in fp32 except on certificate
  iterations (`speed/r32_noncert/`). 4.2% less per wil50 iteration without
  AMX, 2.5% on nug30; about 1% with `DNN_AMX=1` (the `fast` arm above).

Which branching rules the gap thresholds actually reach matters at the DNN
depths. The gap thresholds in `nug_deep.par` send nearly every depth-3 and
depth-4 node (QPB relative gap about 0.02) to the last block, rule 2, so
changing the depth 3-4 block does nothing there. `PARAM/nug_la4.par` ignores the gap and uses rule 4 lookahead (30
candidates) on every node through depth 4, rule 6 to depth 6: on nug28 with
`-x 2:5` it changed the DNN misses by one node at depth 4 (5 to 4) and one
at depth 5 (0 to 1), DNN CPU -2%, within run-to-run variation. A miss at
depth d always has n-d children, so only the number of misses can change.
More QPB lookahead is not a lever. Branching on the DNN marginals is: rule 8
in every block (`nug_deep_xbar_all.par`) cut nug30's depth-3 misses from 85
to 59 and its CPU by 16% (above), though it lost on nug28.

## Correctness checks

`make test` runs both of these. Run it before measuring any change.

**Brute force** (`tests/brute_force.py`, about 10 s). This generates small
random instances (n = 2 to 8) and finds each optimum by enumerating every
permutation. It covers four families: nug-like grids, heavy ties, asymmetric
B, and a linear cost term. It runs the solver under 7 params settings
(`PARAM/default.par`, `rule1`–`rule4`, `glb`, and `tests/params/qp0.par` for
the plain QP bound), each with no upper bound, the optimum as upper bound, and
optimum + 2. Every run must report a root bound no greater than the optimum,
never report an objective below it, report only objectives some permutation
attains, and finish on the optimum. Instances come from a fixed seed. Use
`--seed N` for different ones and `--solver PATH` to test another build.

**Regression** (`bench/bench.py check`, about 15 s on the quick suite). This
reruns a suite and requires the objective, node, fathomed and enumerated
counts, Frank-Wolfe iterations and root bound to match `bench/reference/`
exactly.

- `--tolerance 0.02` allows relative differences, for example when a change
  can resolve ties differently.
- `--suite standard` checks every instance in the reference.
- `--tree DIR` checks a build in another source tree.
- When a change is *meant* to alter the search, regenerate the reference
  with `bench.py check --update` (runs the standard suite), commit
  `bench/reference/` with the change, and say why the counts moved.
- The reference was made on an x86 Xeon. Until Nov 5 an Apple M2 Pro did
  not match it: nug15's root QPB bound was 984.79 against 986.55 and the
  quick suite's node counts differed by up to 8% (nug14 +6%, nug15 -1.5%,
  nug16b -1.5%). The cause was Meschach compiled with fused multiply-adds,
  which ARM always has: an x86 build with them reproduces the M2's numbers
  (`docs/meschach_removal.md`). Meschach is now built with
  `-ffp-contract=off`, which should make the M2 match the reference; until
  `bench.py check` on an M2 confirms it, use `--tolerance` there, or compare
  two builds on the same machine with `bench.py compare`.

The two checks catch different things. We planted two bugs: lower bounds
raised by 3, and an elimination test that was off by one. The brute-force
test failed 720 and 169 runs respectively. On QAPLIB with `--ub opt`, both
bugs *reduced* node counts, and every instance still reported `ok`. An
invalid bound looks like an improvement on real instances, so it's only the
brute-force test that can tell them apart.

## Benchmarking: comparing changes on an apples-to-apples basis

`bench/bench.py` (Python 3, standard library only) runs a suite of instances
under fixed settings, records results, and compares two recorded runs.

```sh
# 1. baseline: build and run the commit you branched from, in a separate worktree
python3 bench/bench.py run --suite quick --rev cognac --label base

# 2. candidate: build and run your working tree
python3 bench/bench.py run --suite quick --build --label mine

# 3. compare (A = baseline, B = candidate); --log records the result in bench/LOG.csv
python3 bench/bench.py compare base mine --log --note "rule 4 uses 30 FW iterations"
```

`--rev` checks out any commit or branch under `bench/.worktrees/` (gitignored),
builds it once, and reuses that binary on later runs. That way the baseline is
always built and run on the same machine as the candidate. Commits from before
the Meschach `isatty` fix (for example `master` at `16a72e8`) don't build with
modern g++. Commits before `0aa6ed8` build, but they have the bound bugs fixed
there, and their node counts include subtrees pruned wrongly. Use `0aa6ed8`
or later as a baseline. See [CONTRIBUTING.md](CONTRIBUTING.md) for when to
log results and how.

Sample comparison:

```
  instance     nodes A     nodes B     B/A      time A    time B     B/A  status
--------------------------------------------------------------------------------
     nug12         249         249   1.000    0.210337   0.20562   0.978  ok
...
geometric mean B/A over 9 instances solved in both:  nodes 1.000   solver time 0.978   (< 1 means B is better)
```

### What is held fixed, and why

1. **Same params file for every instance.** Runs default to `PARAM/default.par`.
   Many per-instance files were hand-tuned, so mixing them would compare
   tuning, not algorithms. To test a params change, point `--params` at a new
   file for the candidate run. `--params auto` reproduces the solver's own
   per-instance default and is meant for regression checks, not comparisons.
2. **Upper bound is set to the known optimum** (`--ub opt`, the default). This
   takes the heuristic's luck in finding good incumbents early out of the
   measurement, so node counts measure the bound plus branching only. This is
   the usual convention when comparing QAP branch-and-bound codes. Use
   `--ub none` to measure the full solve, including incumbent finding.
3. **Nodes are the primary metric.** Node counts are deterministic and don't
   depend on the machine, so you can compare them across machines and days.
   The harness warns if repeats disagree.
4. **Time is secondary and only valid on the same machine.** By default
   `compare` uses solver CPU time (`--time wall` for wall time). Run both sides
   on the same, otherwise idle machine and with the same build flags. Use
   `--repeats 3` to take the median wall time. `compare` warns when hostname,
   CPU, params, params-file contents or upper-bound mode differ between the two
   runs.
5. **Correctness is checked.** Each instance's status is `ok`, `timeout`,
   `error(...)`, or `wrong(...)`. `wrong(below-opt)` means a solution better
   than the known optimum was reported, and `wrong(not-opt)` (with `--ub none`)
   means the search finished without finding the optimum. The command exits
   nonzero on `wrong` or `error`. Ratios only include instances that are `ok`
   in both runs.
6. **Summary uses geometric means** of per-instance B/A ratios, so a 2x
   improvement on a small instance counts as much as one on a large instance.

### Suites

Suites are plain text files in `bench/suites/`, one instance name per line,
with `#` for comments. Times are CPU seconds with `PARAM/default.par` and
`--ub opt` on a 2.3 GHz Xeon.

| Suite | Use | Approx. total time |
|---|---|---|
| `quick` | smoke test after every change | ~15 s |
| `standard` | the run that goes in `bench/LOG.csv` (see CONTRIBUTING.md) | ~3 min |
| `large` | before merging branching, lookahead or bound-strategy changes | ~34 min |

Lookahead (rules 3 and 4 at the top levels) takes about 75% of the time on
`standard` but much less at larger sizes: 50% on nug20, 36% on nug21, and
14–19% on esc16a and chr15a. A speedup that only touches lookahead therefore
looks bigger on `standard` than it is. `large` sets its own 30-minute
per-instance timeout (`#! timeout 1800` in the suite file).

You can also pass instances directly: `bench.py run --instances nug12 had14`.

### What gets recorded

Each run writes `bench/results/<label>/`:

- `results.csv`: one row per instance (status, optimum, objective, nodes /
  fathomed / enumerated, root bound, FW iterations, solver time, median and
  min wall time, params file).
- `meta.json`: git commit and dirty flag, binary hash, params file hashes,
  hostname, CPU, compiler, settings, timestamps, and `--notes`.
- `logs/<instance>.log`: full solver output, starting with the exact command.

`bench/results/` is gitignored because timings only mean something on the
machine that produced them. To share a baseline, `git add -f` its folder.

`compare --log` appends one summary row to `bench/LOG.csv`, which *is*
committed. The row records the date, author, suite and settings, both commits,
how many instances were solved, geometric-mean node and time ratios, the
worst node ratio, any wrong or errored instances, the machine, and your note.
`--log` refuses to write a row when either run had uncommitted changes, or
when the runs differ in params, upper bound, instance list or machine. Pass
`--force-log` to record such a row anyway, and the warnings are saved with it.
