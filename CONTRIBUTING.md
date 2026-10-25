# Contributing

Notes for collaborators working on the branch-and-bound code in `latest/`.
[RUNNING.md](RUNNING.md) covers building, running and the benchmark harness
in detail. This page covers how we work together and how we record
performance results.

## Branches

- `cognac` is the integration branch. Branch from it, and open pull requests
  against it.
- `master` holds the code as presented at OP23 (June 2023) and doesn't build
  with modern compilers. Leave it alone.
- Use one branch per idea, for example `rule-6-lookahead` or `fw-iterations-tuning`,
  so each logged result points at one change.

## Measuring the impact of a change

Every change that could affect the search (bounds, branching rules, params
files, Frank-Wolfe settings, fathoming, symmetry handling) needs a logged
benchmark result before it merges. Pure refactors should show a node ratio of
exactly 1.000, which is a cheap way to prove they really are refactors.

1. **Commit your change first.** `--log` won't record a run made with
   uncommitted changes, because the logged commit has to be the code that
   produced the numbers.

2. **Run the baseline and your branch on the same machine, back to back.**
   The baseline is the commit you branched from, usually `cognac`:

   ```sh
   python3 bench/bench.py run --suite standard --rev cognac --label base-cognac
   python3 bench/bench.py run --suite standard --build      --label my-branch
   ```

   Keep the machine otherwise idle, and don't run the two at the same time.
   You can reuse a baseline run for several candidates as long as it's from
   the same machine and the baseline commit hasn't moved.

3. **Compare and log:**

   ```sh
   python3 bench/bench.py compare base-cognac my-branch --log \
       --note "short description of what changed and why"
   ```

4. **Commit `bench/LOG.csv` on your branch**, together with the change, and
   paste the `compare` table into the pull request description.

Use the `quick` suite (about 15 s) as often as you like while developing. Only
the `standard` suite counts for the log entry that goes with a merge.

### Reading the results

- **Nodes are the headline number.** They're deterministic and don't depend
  on the machine, so a node ratio logged by one person can be compared with
  one logged by another. Lower is better.
- **Time is the tiebreaker**, and only means something between two runs on
  the same machine. A change that cuts nodes but costs more per node, such as
  a more expensive bound, shows up here.
- **Look at the per-instance table, not just the geometric mean.** The log's
  `worst_nodes_ratio` column flags a change that helps on average but blows
  up one instance.
- **Any `wrong(...)` status is a bug.** It means the search reported an
  objective better than the proven optimum, or missed the optimum. Don't merge
  until it's fixed.

### Rules that keep results comparable

- Don't edit `PARAM/default.par` or the suite files in the same branch as an
  algorithm change. If either changes, every earlier log row stops being
  comparable. Make that a separate PR and log a fresh baseline after it merges.
- To test a params change, create a new file (for example
  `PARAM/default-fw200.par`) and pass it with `--params` on the candidate run
  only. The log row will warn that params differ. That's expected, so pass
  `--force-log` and say so in the note.
- Keep the defaults (`--ub opt`, one params file for every instance). If you
  need a different setting, such as `--ub none` to measure the heuristic, run
  both sides with it and mention it in the note.
- Don't delete or rewrite rows in `bench/LOG.csv`. If a row turns out to be
  wrong, append a corrected row and explain it in the note.

## Code changes

- Match the surrounding style: C++ in `.C` / `.h` files, built with the plain
  makefile.
- After editing a header, run `make clean && make` in `latest/`, since the
  makefile doesn't track header dependencies.
- Don't commit build output (`*.o`, `*.a`, `bnbdrv`), `bench/results/` or
  `bench/.worktrees/`. `.gitignore` already excludes them.
