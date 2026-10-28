# Building, running and benchmarking

How to build the branch-and-bound solver in `latest/`, run it on QAPLIB
instances, and compare performance between code or parameter changes.

## Build

Requires `make`, `g++`, and `clang` (the Meschach makefile uses clang; pass
`CC=gcc` to use gcc instead).

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
| `-x D1:D2[:G[:H]]` | DNN bound (latest/dnn.C) on nodes QPB does not fathom at depths D1-D2 with relative gap (inc - QPB)/(inc - root QPB) at most G (default 1); H: give up on a node after H iterations when hopeless (default 0, never). Needs LAPACK (see below) |
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
4     branching rule: 1 U-sum, 2 U-fathom, 3 strong, 4 lookahead, 5 learned (BUFF), 6 eccentric, 7 eccentric lookahead
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

### The DNN bound (-x)

`latest/dnn.C` computes a certified lower bound from the doubly nonnegative
(DNN) SDP relaxation (experiments/dnn_bound). It is far stronger than QPB and
far more expensive, so it pays only near the top of the tree: nug20 with
`-x 2:2` takes 84,592 nodes and 117 s (QPB alone 526,205 nodes, 197 s), with
`-x 2:3` 141 nodes and 105 s. It needs LAPACK: `pip install scipy-openblas32`
gives an OpenBLAS the makefile finds (it falls back to the system
liblapack/libblas, ~7x slower); `make DNN=0` builds without it (then `-x`
is an error). Rebuild after installing or removing the OpenBLAS wheel. Its
floating-point results depend on the BLAS kernels OpenBLAS picks for the
CPU, so node counts with `-x` can differ between machines (rarely: only a
bound within ~1e-6 of the cutoff decides differently). `make dnncheck`
validates the bound itself; `python3 tests/brute_force.py --solver-args
"-x 0:20"` checks it inside the search.

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
