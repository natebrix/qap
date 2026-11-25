# Root runs the DNN+ paper needs

`dnn_plus.tex`'s results tables need root bounds computed under one rule on one named machine.
What exists was computed in different sessions with different limits (sources in the comments
under Table 2 of the .tex).

## Running them (the laptop, Nov 24)

On a checkout of `claude/dnn-plus-paper`, after `make`:

```
caffeinate -i sh -c "printf '%s\n' sko42 sko49 wil50 tho40 | xargs -P 4 -n 1 docs/dnn_plus/root_bound.sh"
caffeinate -i sh -c "printf '%s\n' tai25b tai30b tai35b tai40b tai50b tai60b nug30 | xargs -P 4 -n 1 docs/dnn_plus/root_bound.sh"
```

`root_bound.sh` runs, per problem, the plain bound (4,000 iterations, no give-up) and the bound
with cuts (the recipe's first solve, then 6,000 CP-ADMM iterations; `ITERS=` changes it), each on
one core through `bnbdrv -d 0`, and appends a line per arm to `docs/dnn_plus/runs/results.txt`.
The trace of the arm with cuts has the certified value at every check, so the stopping rule
below is applied afterwards, on the history. Commit `docs/dnn_plus/runs/` when done.

Tested in a cloud container (Nov 24): tai25b both arms (plain 337,736,517 at 4,000 iterations;
with cuts 340,674,316 after only 300 phase iterations), nug30 with cuts (5,952.45 after 200 phase
iterations, above the plain 5,950 of the Nov 5 table). The grid recipes' arms (sko42, sko49,
wil50, tho40) and tai35b-tai60b were not run there: only their params files were generated.
Not in the script: kra30a/b, tho30 (no recipe with a cut phase), tai80b.

## The rule (proposal, to be agreed)

- One machine and build, recorded with each run (commit, BLAS, threads; one core per solve).
- Cold start at the root, no cutoff (a cutoff far above the best known value), no give-up.
- Two values per instance, both certified:
  - **DNN**: the plain solve, stopped when the certified value gains less than 1e-6 of
    |best known| over 500 iterations, or at 10,000 iterations.
  - **DNN+**: the first solve as each recipe has it, then the cut phase in CP form
    (`cut.cp = 100:1000`), `cut.stall = 0`, stopped by the same rule on the certified value,
    or at 10,000 phase iterations.
- Record per run: the certified value at every check (for the bound-against-iterations
  figure), iterations, CPU seconds, the pool's final size, the regions file's hash.

## Runs

| instance | regions | exists today | needed |
|---|---|---|---|
| sko42 | `cut.grid = 6:7`, side 1 | nothing at the root with cuts | DNN, DNN+ |
| sko49 | `cut.grid = 7:7`, side 1 | nothing at the root with cuts | DNN, DNN+ |
| wil50 | `cut.grid = 5:10`, side 1 | node-level values only | DNN, DNN+ |
| tho40 | `cut.grid = 5:8`, side 1 | DNN+ 230,452.99 at 300 + 2,500 iterations, rising (`bound_lab/cpfast`) | DNN, DNN+ to the rule |
| nug30, tho30 | their grids | plain DNN in the 42-instance table | DNN+ (a column for Table 1) |
| tai25b, 30b, 35b, 40b | `bound_lab/regions/regions/*.comb*.txt` | 2,300-iteration values; 6,000-iteration roots for 35b, 40b in `bound_lab/tai_solves` | DNN, DNN+ to the rule |
| tai50b | `infra/phase1/recipes/tai50b.regions` | DNN+ 453,552,676 at 470 + 6,000 | DNN, DNN+ to the rule |
| tai60b | `infra/phase1/recipes/tai60b.regions` | DNN+ 605,315,696.64 (run of record's root) | DNN; DNN+ to the rule |
| kra30a, kra30b | none | held out of the Nov 5 table; `experiments/kra30_solve` quotes 86,837 and 87,858 | confirm, with the penalty that converges |

Cost [e]: 0.3-1.3 s an iteration at n = 40-60 on one core, so 1-4 hours a run at the limits;
about 25 runs.

## Open

- tai80b: the bound lab saw a root above the published 786,900,529 and rising. A row for it
  would be the largest instance in the table; memory and time not sized here.
- Whether a root with strong settings and one with the search's settings are both reported.

## Results (Nov 24, the laptop) and follow-ups

`runs/results.txt` has the eleven problems, both arms, and the tables of `dnn_plus.tex` are filled
from it. The stopping rule above never triggered: at 6,000 phase iterations every DNN+ value was
still gaining 1e-5 to 1e-4 of the best known value per 500 iterations, so the tables report the
value at the limit. Not blocking the first draft (Nathan, Nov 24); to do later:

- **Longer runs**: done Nov 25 (`runs_long/`): wil50 48,302.31, tho40 230,659.65, tai50b
  455,919,510.84 at 20,000 phase iterations; in the paper with `fig_history.py`'s figure.
- **Timings one at a time**: done Nov 25 (`runs_serial/`); the bounds are identical to Nov 24's.
- **A history for the plain arm**, for the bound-against-iterations figure (the plain solve does
  not record one; the cut phase does).
- **kra30a, kra30b, tho30**: done Nov 24-25 in a cloud container (Linux x86_64, commit 5b21644b;
  `runs_container/results.txt`). Plain DNN after 4,000 iterations: kra30a 86,836.56, kra30b
  87,856.27, tho30 143,573.65. DNN+ for tho30 (3 x 10 grid, 6,200 iterations): 146,123.70.
  The values are in the paper; the container's timings are not for the paper.
- **Warm starts and inherited pools**: to be redone with the current code (Nathan, Nov 25), by
  `warm_pool.sh` (the header has the arms): the children of one fixed node per problem, bounded
  cold, from the parent's Z, from Z and Y (no cuts anywhere), and then with and without the
  parent's cut pool after a 2,000-iteration phase at the parent. Cases: nug30 `0:14`, tho40
  `26:10,8:0,33:32`, wil50 `27:39,48:30,49:25`, tai35b's root (pool arms only).
  `warm_pool_summary.py` prints the tables. Done Nov 25 on the laptop (`runs_warm/`); in the paper as Table `tab:warm`.
- **QAPLIB's current listings**, by hand: the site is not reachable from cloud sessions.

## The 42-instance table on the laptop (proposed Nov 25)

`root_bound.sh` now takes the small and medium problems of the paper's first table (DNN arm only,
nug30's settings, 4,000 iterations). One command reruns the whole table on the laptop:

    OUT=$PWD/docs/dnn_plus/runs_small docs/dnn_plus/root_bound.sh esc16a esc16b esc16c esc16d esc16e \
      esc16g esc16h esc16i esc16j had12 had14 had16 had18 had20 kra30a kra30b nug12 nug14 nug15 nug16a \
      nug16b nug17 nug18 nug20 nug21 nug22 nug24 nug25 nug30 rou12 rou15 rou20 scr12 scr15 scr20 tai12a \
      tai15a tai17a tai20a tai25a tai30a tho30

Checked in a cloud container on esc16a, nug12, had12, scr12. With one setting for every problem some
values come out lower than the Nov 5 table's (esc16a 62.56, so 63 against 64): the Nov 5 runs used
other settings. Decide: one setting for all, or settings per family.
