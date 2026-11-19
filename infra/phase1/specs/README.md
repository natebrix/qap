# Run specs: a launch from a committed file

`infra/phase1/start.sh --spec infra/phase1/specs/NAME.spec` launches what the file says, with no
questions but the last ("Launch?"; `--yes` answers it). The file must be committed as it is: the
run's record (`runs/RUN/spec.txt`, a copy of it) names it at that commit, and the record is
committed and pushed once the run is launched (`--no-commit`: not). Why: on Nov 18 an A/B arm went
out with a unit-sizing line cut short at a prompt, so its arms differed in tracing as well as in
`cut.zone` (`docs/review_nov18_runs.md`, items 2 and 5).

A spec is a shell file of `KEY=value` lines (sourced, as `recipes/*.sh`); an unknown key is an error.
A key it leaves out takes `start.sh`'s default, as Enter would.

| key | what | default |
|---|---|---|
| `PROBLEM` | the QAPLIB name (required) | |
| `RUN` | the run's name; with `ARMS`, each arm's is `RUN-ARM` | `PROBLEM-YYMMDD` |
| `PARAMS` | the params file, relative to the repository | `recipes/defaults.txt`'s |
| `ARMS` | an A/B: `"NAME=PARAMS NAME=PARAMS ..."`, one run per arm, identical but for the params file (instead of `PARAMS`) | |
| `INCUMBENT` | a number, or `none` | the best known + the grid's step |
| `DATA`, `FLAGS` | the data file; the solver flags (`none`: none) | the grid copy and its symmetry |
| `UNITS` | the unit sizing (`qapctl init`'s options) | 20-minute units, `--trace` |
| `VISIBILITY` | `public` or `private` | `public` |
| `ESTIMATE` | `SOLVE_HOURS:DIVE_SECONDS[:B[:SAMPLE[:unif\|qpb]]]`, or `none` for a solve | `none` |
| `ESTIMATE_ENV` | more estimator settings for every dives unit, `QAP_EST_X=V` separated by spaces: `QAP_EST_MIN_DIVES=5` (dives every branched open node gets whatever the budget), `QAP_EST_SAMPLE_DEPTHS=3:1,4:0.2,5:0.03` (the sampling rate by depth) | none |
| `WORKERS` | the worker group's size | 1 |
| `WORKER_TYPE`, `WORKER_TYPES` | the worker type; the group's types | `c7a.16xlarge`; it and `m7a`, `r7a` |
| `MARKET` | the workers: `spot` or `on-demand` | `spot` |
| `WORKER_SLOTS`, `COORD_TYPE`, `COORD_SLOTS` | slots per worker; the coordinator (on demand) and its slots | the cores (by memory); `c7a.2xlarge`, 7 |
| `HOURS`, `DISK` | the time limit; GB per instance | 72; 100 |

Before a spec launches anything: `start.sh --spec F --dry-run --yes` prints every AWS call it would
make. `collect.sh RUN` (or `RUN-ARM`) then shows the run, and its readout once it is done.

| spec | what |
|---|---|
| `nug30-smoke.spec` | the first run on worker groups: nug30 (about 1,300 CPU s) on one worker; checks the group's start, the worker joining, the group's 0 at the end |
| `tai50b-est.spec` | an estimate of tai50b under its current recipe (rule 9 at the bottom), whose tree is known: tai50b-r9's 6,629 nodes and 117 core-hours |
| `sko42-zone.spec` | the `cut.zone` A/B on sko42, 0.010 against 0.025, from scratch (the Nov 18 arms differed in tracing and were both reclaimed) |
| `tho40-est.spec` | a second estimate of tho40 under tho40c with the first's diagnosis applied: 5 dives at least per branched open node, sampling by depth |
