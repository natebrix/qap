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
| `FROM` | a paused run (its units back, queued) this run starts from, with `DIVES` or `SAMPLE`: its states and parents are read, never changed, so it can still be resumed | none |
| `DIVES` | with `FROM`: this run is a dives run of FROM's paused frontier, `DIVE_SECONDS[:B[:SAMPLE[:unif\|qpb]]]` as in `ESTIMATE` (no solve phase: a second estimate of the same open nodes, e.g. with the estimator fixed; FROM's params and incumbent, so no `PARAMS`, `ARMS` or `INCUMBENT`; `ESTIMATE_ENV` for the estimator's settings). Its readout is its own | none |
| `SAMPLE` | with `FROM`: this run starts from a sample of FROM's open nodes, `DEPTH:COUNT,...`, drawn in proportion to the square root of their subtrees' estimates in a dives run of FROM, and solves them to the end under `PARAMS` or each arm of `ARMS` (`coord/sample.py`: every arm the same nodes; `experiments/estimates_nov19/paired.py` compares them). The data file, flags and incumbent must be FROM's | none |
| `SAMPLE_BY`, `SAMPLE_SEED` | with `SAMPLE`: the dives run whose estimates size the nodes; the draw's seed | FROM's own dives run; 1 |
| `ARMS_ENV` | with `ARMS`: an arm's own solver environment, `NAME:K=V` items separated by spaces (an experiment setting no params file holds, e.g. `cross5:DNN_ADMM_WARM_BETA=0.1:30:0.36:1`) | none |
| `ESTIMATE_DECIDE` | the estimate's decision point (`coord/README.md`): when its pilot says the dives would cost over this share of the rest of the tree, they stop and the solve goes on to the end; `off`: the dives always run to the end | `0.15` |
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
| `sko42.spec` | sko42 to the end under `sko42.par` of Nov 19: `cut.zone 0.025` and the bound lab's sko49 policy block (also in `sko49.par`); a guide to sko49, and the block's whole-tree effect against `sko42-zone-261119-z25` (`experiments/estimates/sko42/zone_nov19/`) |
| `tho40-est.spec` | a second estimate of tho40 under tho40c with the first's diagnosis applied: 5 dives at least per branched open node, sampling by depth |
| `tho40-est3.spec` | tho40-est2-261119's paused frontier dived again with the estimator fixed (strata states spilled, not dropped; Nov 19): est2's settings otherwise, so the difference is the fix (est2: 239 c7a CPU-days read, 196 corrected) |
| `tho40-tune.spec` | tho40's first tuning round: four arms (tho40c; `branch.xbar_lines 4`; `dnn.cross_first 5` with a warm penalty schedule; `cut.zone 0.03`) on the same sample of tho40-est2's open nodes, 50 at depth 3 and 100 at depth 4, side by side; about 150 CPU h an arm |
