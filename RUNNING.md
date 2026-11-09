# Building, running and benchmarking

How to build the branch-and-bound solver in `latest/`, run it on QAPLIB
instances, and compare performance between code or parameter changes.

## Build

Requires `make` and `g++`. The DNN bound (`-x`), which is built by default,
also needs LAPACK (see `latest/makefile`) and a C compiler for
`latest/amx_gemm.c` (`cc`, which comes with g++ on the usual distributions).
`make DNN=0` builds without the DNN bound and needs neither.

```sh
make          # from the repo root: builds latest/bnbdrv (same as cd latest && make)
make test     # brute-force tests + regression check (see "Correctness checks")
```

To compile the solver with clang instead: `make CPP=clang++` (`latest/makefile`
calls the C++ compiler `CPP`). clang 18 crashes on `latest/amx_gemm.c` (an
LLVM bug with the AMX/bf16 flags), so `CC=clang` works only with `DNN=0`.
The solver no longer needs the Meschach library that used to be in
`src/MESCHACH/` (`docs/meschach_removal.md`): its matrix storage is
`latest/mat.h` and its eigensolver `latest/symeig.C`.

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
| Initial upper bound | none (1e100) |  The solver finds incumbents itself and prints `New objective value: …` each time it improves. |

Other flags (`./bnbdrv` with no arguments prints help):

| Flag | Meaning |
|---|---|
| `-p F` | params file |
| `-s F` | symmetry file |
| `-u X` | initial upper bound (incumbent); rounded up to an attainable objective value (even, or integer) when it is not one |
| `-U F` | incumbent file `obj [perm]`, re-read every 50 nodes and on SIGUSR2; a lower value is adopted |
| `-F i:j,...` | start at the node with these assignments fixed; its depth is the number of fixes, so `-d D` counts from the root, not from this node |
| `-d X` | stop at tree depth X (overrides the params file) |
| `-f X` | scale reported times by X |
| `-h` | run the heuristic only (no branch and bound) |
| `-e` | estimate tree size with Knuth-style dives (`latest/estimate.C`). Reads `./estimate.par`: six numbers, `start_dive_level trials print_every print_everything_every score_exponent ignore_depth` (the checked-in file is `2 1000 100 100 1.5 1000`). Without the file it silently uses the defaults in `estimate.h` (10,000 trials, exponent 1.0). Pass the known optimum with `-u`: dives are weighted by it. `rand()` is seeded only when `QAP_EST_SEED` is set, so set a different seed per process when running several in parallel; `QAP_EST_HIST=F` writes the per-dive history (depth, node estimate, time estimate) to F instead of `dive_history_nug25_2.txt` |
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
(`tests/parallel_test.py`). `--checkpoint F` periodically rewrites F, an
ordinary state file holding every open node (the pool plus the nodes the
running tasks were given) and the totals of the finished tasks, without
stopping any worker; `--resume F` continues it (redoing at most the tasks that
were in flight), or `bnbdrv -r F` finishes it serially. SIGTERM writes F the
same way before stopping the workers, and `--wall-limit S` stops the workers
at their next node after S seconds, writes F and exits 3. (Before Nov 8 2026
a checkpoint stopped every worker and launched nothing until the last had
exited, which on long-task instances such as tho40 idled 70% of the machine.) `bench.py run --jobs N` runs each
instance this way; compare such runs with `compare --time wall`.

On this 4-vCPU container nug16a ran 3.2x faster in wall time with 4 workers and
nug17 2.7x; total CPU rises somewhat under full load, so expect less than
N-fold on shared machines.

Known optimal solutions are in `QAPSOLN/<name>.sln`. The first line is `n` and
the optimal objective, and the second line is the permutation (absent for
esc32d, kra30a, kra30a-old, kra30b, kra30b-old, tai25a and tho30, whose files
held `<value> ??`; they, nug27, nug28 and nug30, which lacked n on the first
line, were repaired on Nov 8 2026, and kra32's was rebuilt from a tabu
solution checked against the data file). `bench.py` reads only a first line
of exactly `n optimum` and treats anything else as unknown.

Three data files are scaled variants, not the QAPLIB originals: `qaplib/kra30a.dat`,
`kra30b.dat` and `kra32.dat` have the flow and distance matrices swapped and every
distance divided by 5, so their optima are 17,780, 18,284 and 17,740 (QAPLIB's
88,900, 91,420 and 88,700 divided by 5) and every bound on them is one fifth of
the literature's. The unscaled originals are `kra30a-old.dat`, `kra30b-old.dat`
and `kra32-old.dat`. Pass the matching `--ub` (17780 for `kra30a.dat`, 88900 for
`kra30a-old.dat`).

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

When the DNN bound ran (`-x`), three more tables follow (Nov 9, state version
7, which saves their counts). From nug24 with the nug30 recipe:

```
DNN bound by level:
level    nodes	%fathom	%gaveup	      iters	 it/node	    time
 1           6	0.5000	0.5000	       1193	   198.8	   11.12
 2          69	0.9275	0.0725	       3017	    43.7	   26.47
 3         110	0.9727	0.0273	       3175	    28.9	   22.76
 4          63	1.0000	0.0000	       1404	    22.3	    8.00

Fathoms by cause:
level  fathoms	     qpb	     dnn	    cuts	  popped
 0           0	       0	       0	       0	       0
 1           3	       0	       3	       0	       0
 ...
```

`DNN bound by level`: the nodes the DNN bound ran on, the fractions it
fathomed and gave up on as hopeless (the `-x` give-up rule), its ADMM
iterations and CPU seconds. `Cut phase by level` (only with
`QAP_DNN_CUT_ZONE`): nodes, fraction fathomed, CPU. `Fathoms by cause` splits
each depth's fathoms: by QPB, by the DNN bound, by the cut phase, and nodes
popped after the incumbent improved past their bound (`qpb` is what is left
after the other three; a run resumed from a state file before version 7 has no
causes for its first part). Runs without `-x` print exactly what they did.

`bnbdrv -p P -r F instance` on a state file with no open nodes prints its
statistics as a finished run's report and stops. `bench/parallel.py` and the
coordinator (`coord/`) use it to print the solver's own tables after their
summed RESULTS head.

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
| nug30 record setting (Nov 4) | 3,844 (DNN: 3,843, 3,701 fathomed, 1,826 s) | 1,838 | 5.4 min (experiments/nug30_solve) |
| **nug30 record setting (Nov 7, below)** | **3,934 (DNN: 3,933, 3,788 fathomed, 1,358 s)** | **1,370** | **3.3 min** (experiments/nug30_solve) |
| kra30a, the nug30 record setting, `--ub 17780` (scaled file; Nov 7) | 198 (DNN: 197, 190 fathomed, 100 s) | 101 | 18 s (experiments/kra30_solve) |
| kra30b, the nug30 record setting, `--ub 18284` (scaled file; Nov 7) | 3,446 (DNN: 3,445, 3,315 fathomed, 1,252 s) | 1,261 | 3.1 min (experiments/kra30_solve) |
| tho30, the nug30 record setting, `--ub 149936` (Nov 7) | 18,410 (DNN: 18,403, 17,680 fathomed, 6,354 s) | 6,401 | 15.6 min (experiments/tho30_solve) |
| **nug30 record setting + rule 9 (`nug_deep_xbar9_all.par`; Nov 8)** | **3,705 (DNN: 3,704, 3,563 fathomed, 1,256 s)** | **1,270** | **3.0 min** (experiments/nug30_solve) |
| tho30, record setting + rule 9 (Nov 8) | 14,822 (DNN: 14,819, 14,213 fathomed, 4,683 s) | 4,733 | 10.9 min (experiments/tho30_solve) |
| tho30, record setting + side-1 cut phase through depth 2 (Nov 8) | 17,786 (DNN: 17,779, 17,056 fathomed, 6,113 s) | 7,172 | 17.1 min (experiments/tho30_solve) |
| tho30, the same with the CP form (`QAP_DNN_CUT_CP=100:1000`; Nov 8) | 17,786 (DNN: 17,779, 17,056 fathomed, 5,957 s) | 6,994 | 17.3 min (experiments/tho30_solve) |
| kra30b, record setting + rule 9, `--ub 18284` (Nov 8) | 1,686 (DNN: 1,685, 1,620 fathomed, 644 s) | 651 | 1.8 min (experiments/kra30_solve) |
| kra32, record setting + rule 9, `--ub 17740` (no symmetry file; Nov 8) | 10,807 (DNN: 10,806, 10,440 fathomed, 5,554 s) | 5,603 | 13.8 min (experiments/kra30_solve) |
| kra30a-old, record setting, `--ub 88900`, no symmetry file (Nov 8) | 873 (DNN: 872, 842 fathomed, 369 s) | 373 | 56 s (experiments/kra30_solve) |
| kra30b-old, record setting, `--ub 91420`, no symmetry file (Nov 8) | 28,112 (DNN: 28,111, 27,072 fathomed, 17,929 s) | 18,063 | 47.4 min (experiments/kra30_solve) |
| kra32, record setting, `--ub 17740` (scaled file, no symmetry file; Nov 8) | 10,852 (DNN: 10,851, 10,484 fathomed, 7,368 s) | 7,406 | 17.7 min (experiments/kra30_solve) |
| sko42t (grid second, `SYM/sko42t.sym`), record setting, `--ub 15812` (Nov 8, stopped after 5.5 h) | 9,651 bounded, 9,318 fathomed, 3,168 open and growing | 32,988 so far | checkpoint kept, not finished (experiments/sko42_solve) |
| tho40t (grid second, `SYM/tho40t.sym`), record setting + rule 9, `--ub 240516`, traced (Nov 8, stopped after 5 h) | 29,254 bounded, 28,271 fathomed, 5,440 open | 40,691 so far | not finished; with a depth-2 sweep, tree estimated at 6e5-3e6 nodes (experiments/tho40_estimate) |

With plain `-x 2:5` the trees match the Xeon runs exactly (nug27, nug28).
With the give-up threshold and warm starts (`-x ...:300:1`) a few decisions
at the cutoff go the other way: nug25 took 750 nodes instead of 717, and on
nug30 two depth-7 nodes missed where the Xeon fathomed all 96, so 580,000
QPB-only nodes followed below depth 7 (84 CPU s, 1% of the run). Per DNN
node the laptop is about 3x faster than the 2.8 GHz Xeon (depth 3: 1.9 s
against 5.8; depth 4: 1.2 against 3.5).

**The nug30 record setting** (Nov 7: 1,370 CPU s and 3 min 16 s on the
M2 Pro with 8 workers; the Nov 5 subset of it, without the last line of
switches, took 3,270 CPU s and 14 min 9 s on 4 cores of a 2.1 GHz Xeon; how
each switch got there is in `experiments/nug30_solve/README.md`) **is the
solver's default since Nov 8 2026**: with no `QAP_DNN_*`/`DNN_*` variable set,

```sh
python3 bench/parallel.py qaplib/nug30.dat --jobs 8 --ub 6124 \
    --params PARAM/nug_deep_xbar9_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint nug30.ckpt --checkpoint-every 300
```

gives the record tree (3,716 nodes with rule 9, verified node for node
against the explicit environment). The environment still overrides every
switch, and setting the recipe explicitly is harmless. The defaults that
changed, old -> new, for anyone whose A/B control was "unset":

| variable | old default | new default | variable | old default | new default |
|---|---|---|---|---|---|
| `QAP_DNN_GAMMA1` | 0 | 1.0 | `DNN_CHOLQR_PASSES` | 2 | 1 |
| `QAP_DNN_GAMMA` | 1.6 | 1.0 | `DNN_R32_NONCERT` | 0 | 1 |
| `QAP_DNN_CHEB` | 2 | 3 | `DNN_ADMM_CROSS` | 0 | 2 |
| `QAP_DNN_TIGHT_LAM` | 0 | 2 | `DNN_CORE_LAZYCERT` | 0 | 1 |
| `QAP_DNN_BETA` | 0.3 | 0.21 | `DNN_CORE_NOSYMV` | 0 | 1 |
| `QAP_DNN_REFRESH` | 20 | 120 | `DNN_NOSYMW` | 0 | 1 |
| `QAP_DNN_WARM_Y` | 0 | 2 | `DNN_CORE_REFRESH_PARTIAL` | 0 | 1 |
| `QAP_DNN_EXT` | 0 | 20 | `DNN_BETA_SCHED` | unset (off) | `w=0.36` (`""` turns it off) |
| `QAP_DNN_GU_ZONE` | 0 | 0.003 | `QAP_DNN_SHIP` | 0 | 32 |

The `QAP_DNN_*` defaults apply to the search (`bnbdrv -x`); the `DNN_*` ones
live in `dnn.C` and so also apply to `dnntest`, whose cold solves are
unaffected by `DNN_BETA_SCHED=w=...` (a warm-child factor). The cut phase
(`QAP_DNN_CUT_*`, `QAP_DNN_POOL_CP`, `DNN_CP_LAZYCERT`) stays opt-in. The
previous form of this block, for the record:

```sh
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
       QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 \
       QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32 DNN_R32_NONCERT=1 \
       DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 \
       DNN_CORE_REFRESH_PARTIAL=1 DNN_BETA_SCHED=w=0.36
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
  penalty is 0.21. The warm factor swept Nov 6 on 24 wil50 warm children
  (`speed/nov6_admm/`): the response is monotone from 0.85 down to 0.15, and
  `w=0.25` fathoms the same children in 43% fewer iterations than 0.7 with
  every miss higher; the key `wc=F` gives the rounds of a cut phase their own
  factor (unset: `w` for both; before the Nov 6 fix, pooled children took `wc`
  too). In the search (`pipeline_ab/`, laptop, Nov 6): `w=0.25,wc=0.7` gives
  1.24x fewer DNN iterations and 1.12x less DNN CPU on the node-10 subtree's
  232 matched nodes, so the wil50 production candidate is now
  `DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7`. `w=0.25` for the cut
  rounds as well (`wlow`) loses: the rounds take 1.5-1.9x the iterations.
- `DNN_ADMM_CROSS=2` (unset): between two regular certificates (every 10
  iterations), certify once more at the iteration where the running best is
  predicted to cross the cutoff (`speed/nov6_admm/`): fathoms come 0-9
  iterations earlier (4.2 per fathom on 24 warm children, 1.05x their CPU),
  nothing else changes. In the search: 1.035x in DNN iterations and 1.04x in
  DNN CPU (`pipeline_ab/`, Nov 6). `DNN_ADMM_CROSS_FIRST` and
  `DNN_ADMM_WARM_DAMP` are the same study's other probes (measured not worth
  it); all off by default.
- `DNN_ADMM_WARM_BETA=F1:T:F2` (unset): a node's own warm solve runs at F1 x
  the base penalty until the first full eigendecomposition at or after
  iteration T, then at F2 x it, in place of `DNN_BETA_SCHED`'s `w`; cut rounds
  keep `w` or `wc`. Under test Nov 6 as `0.25:60:0.7` (`pipeline_ab/` arm
  `wtwo`): children fathom as at 0.25, and a miss hands its cut phase a state
  closer to 0.7's.
- `DNN_R32_NONCERT=1` (unset): R's rank term in fp32 except on certificate
  iterations (`speed/r32_noncert/`). 4.2% less per wil50 iteration without
  AMX, 2.5% on nug30; about 1% with `DNN_AMX=1` (the `fast` arm above).
- `DNN_CORE_LAZYCERT=1` (unset): Rump's Cholesky test of the certificate only
  at checks that can decide something (the first, any whose polyhedral value
  could clear the cutoff, the last iteration, and the stopping check and the
  stored best at the exit); the same fathom decisions and returned lb by
  construction (`speed/nov6_core/`).
- `DNN_CORE_REFRESH_PARTIAL=1` (unset): the subspace tracker's full
  eigendecomposition back-transforms only the columns it reads (bit-identical
  eigenvectors, about a third of the refresh saved; `speed/nov6_core/`).
- `DNN_CORE_NOSYMV=1` (unset; with `DNN_NOSYMW=1`): no symmetrisation pass
  after the expand transform; the iterate is symmetric to rounding only and
  the certificate evaluates the symmetric part exactly (`speed/nov6_core/`).
  Changes the trajectory at the rounding level.
- `DNN_CORE_AMX_IB/KB/NT/ORDER` (unset): the AMX kernel's blocking, W load
  hint and product order, for microbenchmarks; every setting gives the same
  product bit for bit (`speed/nov6_core/`).

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
- The reference was made on an x86 Xeon, and an Apple M2 Pro matches it
  exactly (checked Nov 5). Before Nov 5 the M2 did not: Meschach was
  compiled with fused multiply-adds, which ARM always has, and nug15's root
  QPB bound came out 984.79 against 986.55. The fix was to build Meschach
  with `-ffp-contract=off` (`docs/meschach_removal.md`, stage 0); since
  Nov 6 Meschach is gone and its replacements are compiled with the
  solver's flags, `-ffp-contract=off` included. Laptop node
  counts from before Nov 5 (the nug30 record's 3,844, for one) came from
  the old build and will not reproduce exactly.

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
always built and run on the same machine as the candidate. For a commit from
before Meschach's removal (Nov 6), it builds `src/MESCHACH` first, with
`MESCHACH_CC` as the C compiler if that is set. Commits from before
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
   runs. The solver runs with one BLAS thread (`OPENBLAS_NUM_THREADS`,
   `OMP_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` set to 1 unless already set):
   until Nov 5 an idle OpenBLAS pool added about 0.36 s of CPU to every run,
   so CPU times of sub-second instances in older results are inflated by that.
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
