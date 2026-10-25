# Building, running and benchmarking

How to build the branch-and-bound solver in `latest/`, run it on QAPLIB
instances, and compare performance between code or parameter changes.

## Build

Requires `make`, `g++`, and `clang` (the Meschach makefile uses clang; pass
`CC=gcc` to use gcc instead).

```sh
cd src/MESCHACH && make basic     # builds meschach.a (old-style C warnings are expected)
cd ../../latest && make           # builds latest/bnbdrv
```

After changing solver code, run `make` in `latest/` again. The makefile has no
header dependencies, so after editing a `.h` file run `make clean && make`.

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
| Params file | `./PARAM/<name>.par` | **If this file doesn't exist, the solver does not complain.** It runs with uninitialized parameters. Only a few instances (nug*, kra30b, esc32d, rou15, scr15, tai17a, tai20a, ste36a, …) have their own params file, so pass `-p` for anything else. |
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

Known optimal solutions are in `QAPSOLN/<name>.sln`. The first line is `n` and
the optimal objective, and the second line is the permutation.

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
4     branching rule: 1 U-sum, 2 U-fathom, 3 strong, 4 lookahead, 5 learned (BUFF)
150   Frank-Wolfe iterations for the bound
100   Frank-Wolfe iterations when we know we can't fathom
50    Frank-Wolfe iterations for row selection
1     bound: 0 QP, 1 QP parametric, 2 QP+GLB, 3 QP improved, 4 GLB, 5 PB, 6 EVB3, 7 BUFF
30    step for parametric bound
30    compute B for the best _ rows of U
```

Rule 5 (`BRANCH_BUFF`) reads coefficients from `latest/lgb_init.txt`, which
isn't in the repo. Don't use rule 5 without that file.

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
modern g++, so use `cognac` or later as a baseline. See
[CONTRIBUTING.md](CONTRIBUTING.md) for when to log results and how.

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
| `standard` | the run that goes in `bench/LOG.csv` (see CONTRIBUTING.md) | ~5 min |

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
