# Root runs the DNN+ paper needs (spec, not launched)

`dnn_plus.tex`'s results tables need root bounds computed under one rule on one named machine.
What exists was computed in different sessions with different limits (sources in the comments
under Table 2 of the .tex). Nothing here has been run; the commands are to be written by
whoever launches, from the settings named.

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
