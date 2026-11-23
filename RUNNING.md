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
| `-e` | estimate tree size with Knuth-style dives (`latest/estimate.C`). Reads `./estimate.par`: six numbers, `start_dive_level trials print_every print_everything_every score_exponent ignore_depth` (the checked-in file is `2 1000 100 100 1.5 1000`). Without the file it silently uses the defaults in `estimate.h` (10,000 trials, exponent 1.0). Pass the known optimum with `-u`: dives are weighted by it. `rand()` is seeded only when `QAP_EST_SEED` is set, so set a different seed per process when running several in parallel; `QAP_EST_HIST=F` writes the per-dive history (depth, node estimate, time estimate) to F instead of `dive_history_nug25_2.txt`. `QAP_EST_SURV=b` (Nov 18) replaces Knuth's dives with survivor dives: every node of the start stack (`-e -r F`: the state's open nodes) is a stratum, bounded with its children; a dive goes on into a child that branches (weights exp(b x), x the shortfall over the depth's median), bounding all of each node's children, so it ends where the tree does. `QAP_EST_BUDGET` (wall seconds), `QAP_EST_DIVES` (JSON lines: coord/dives.py adds them up), `QAP_EST_WARM_KEEP`; experiments/estimator_replay for why |
| `-b` | benchmark mode |
| `-c` / `-a F` | interactive console / run console commands from file F |
| `-l F` | log file |
| `-o F` | write the search state to F: at checkpoints, on `-t`, and at the end (with `-d D`, the open nodes at depth D) |
| `-k X` | with `-o`, checkpoint every X CPU seconds |
| `-t X` | with `-o`, stop after X CPU seconds, saving the state; exit code 3 |
| `-r F` | resume from state file F instead of the root |
| `-x D1:D2[:G[:H[:W]]]` | DNN bound (latest/dnn.C) on nodes QPB does not fathom at depths D1-D2 with relative gap (inc - QPB)/(inc - root QPB) at most G (default 1); H: give up on a node after H iterations when hopeless (default 0, never); W=1: start a node's ADMM from its parent's final multiplier when the parent's DNN bound missed (default 0; on nug25 -x 2:5:1:300 it cut iterations per fathomed node ~40% but DNN CPU only 4%). Needs LAPACK (see below) |
| `-w W` | Frank-Wolfe patience at lookahead/strong levels (default 0 = off; see `bench/LOG.csv`) |
| `--print-settings` | print the resolved settings (below) and exit; needs no instance. With `-p F` the params file's hash is included, with `-x`/`-w` their values |
| `--allow-settings-change` | with `-r`: resume a state file written under other settings (a warning instead of a refusal) |

A flag without its argument, a number that does not parse (`-u 3O00`), an
unknown option, two instance files, `-F` with `-r`, and `-t` or `-k` without
`-o` are errors (exit 1): the run stops instead of doing something other than
what was asked.

### Settings

Every setting the solver reads is an entry of one table (`latest/settings.C`:
key and/or environment name, type, default, range, class, one line of
description; `docs/settings_plan.md` is the plan, `CONTRIBUTING.md` "Where
settings live" the rule). The table is read once at startup, so every reader
sees the same value. Where a setting lives follows from its class:

- **policy** (changes the search): a key in the params file's `[settings]`
  section, after the positional blocks ("Params file format" below). The
  families are `dnn.*`, `cut.*`, `pool.*`, `ship.*`, `branch.*`, `sb.*`, `est.*`,
  `fw.*`. The flags `-x` and `-w` are aliases of `dnn.depths`, `dnn.maxgap`,
  `dnn.giveup`, `dnn.warm` and `fw.patience`, and override the file.
- **experiment** (an A/B switch in flight): an environment variable (`DNN_*`).
  It graduates to a key or a default, or goes.
- **run** (this run's identity): the `-p` path, `QAP_EST_SEED`,
  `QAP_ALLOW_UNKNOWN`, `QAP_SB_PROGRESS`.
- **diagnostic** (output only): `QAP_DNN_TRACE`, `QAP_NODE_TRACE` (a line
  per node: format in `latest/bnb.C`), `QAP_DNN_POOL_CHECK`, `QAP_XBAR_DEBUG`,
  `QAP_EST_HIST`, `QAP_SB_TRACE`.

Each run prints the resolved table after it has read the params file, as the
dump `bnbdrv --print-settings [-p file.par] [-x ...]` gives on its own:

```
settings 64 hash ebd9848cb3f7e0d2
setting dnn.depths 1:4 par policy
setting dnn.maxgap 1 default policy
...
setting DNN_R32 0 default experiment
setting par.hash 3f2f813daeeb4a6b par policy
```

One line per entry: name, value (`""` when empty), source (`default`, `env`,
`flag`, `par`) and class. The hash covers every policy and experiment entry,
defaults included, so a changed default changes it.

- **An unknown name is an error.** A `[settings]` key not in the table, a
  variable in the environment whose name starts with `QAP_DNN_`, `QAP_EST_`,
  `QAP_XBAR_`, `QAP_NODE_`, `QAP_SB_` or `DNN_` and is not in the table, and a value that does not
  parse or is out of range all stop the run with exit code 2 and a message
  naming it. The environment name a policy setting had before Nov 9
  (`QAP_DNN_BETA`, `DNN_BETA_SCHED`, ...) is such an error, and the message
  names the key; `bnbdrv --settings-map` prints the whole map (one `OLD KEY`
  per line) and `bench/envpar.py` turns an old environment recipe into a
  params file with a `[settings]` section. `QAP_ALLOW_UNKNOWN=1` makes unknown
  and old names a warning, for bisecting old commits. Other `QAP_*` names
  (`QAP_ROLE`, `QAP_BUCKET`, ...) belong to the launchers and are not checked.
- **State files carry the record** (`QAP-STATE 8`: the hash and the entries
  set away from their defaults). `-r` refuses a file whose hash differs from
  the run's own settings, printing both records, unless
  `--allow-settings-change` is given. Version 7 and older files have no
  record and load as before. `bench/parallel.py` keeps the line verbatim in
  the task files it writes, so its workers are checked too.
- **Launchers check a golden dump.** `fullrun/settings.sh`'s `fullrun_check`
  diffs the binary's dump under `fullrun/fullrun.par` against the committed
  `fullrun/settings.golden`; any difference, a changed default included,
  fails it. The sampling launcher writes the dump to `settings.txt` beside
  each run's `env.json`.
- **`dnntest` reads the same keys** (`dnntest -p file.par`), its own flags
  (`-B`, `-g`, `-G`, `-H`, `-R`, `-M`, `-Z`, `-r`, `-d`, `-L`) overriding the
  file, so a cold node solve and the search can run from one file.

The table, as `bnbdrv --settings-table` prints it (`make test` fails when this
copy differs; `tests/settings_doc.py --update` rewrites it):

<!-- settings-table-begin -->
| setting | where | type | default | range | class | what |
|---|---|---|---|---|---|---|
| `dnn.depths` | par key or `-x` | text | unset | text | policy | -x dmin:dmax: the DNN bound's depth range (unset: no DNN bound) |
| `dnn.maxgap` | par key or `-x` | number | `1` | 0..inf | policy | -x third field: relative gap (inc - QPB)/(inc - root QPB) at most this |
| `dnn.giveup` | par key or `-x` | integer | `0` | 0..100000000 | policy | -x fourth field: give up on a hopeless node after this many iterations (0: never) |
| `dnn.warm` | par key or `-x` | integer | `0` | 0..1 | policy | -x fifth field: 1: a child's ADMM starts from its parent's final state |
| `dnn.giveup_warm` | par key | integer | `0` | 0..100000000 | policy | a warm-started node gives up after this many iterations instead of dnn.giveup (0: the same floor; the cut phase is unaffected) |
| `dnn.ext` | par key | integer | `20` | 0..1000 | policy | levels below the last DNN depth at which a child of a DNN miss keeps the bound |
| `dnn.gamma1` | par key | number | `1.0` | 0..2 | policy | Peaceman-Rachford first dual step (0: none) |
| `dnn.gamma` | par key | number | `1.0` | 0..2 | policy | ADMM dual step length |
| `dnn.cheb` | par key | integer | `3` | 1..12 | policy | Chebyshev filter degree of the subspace eigensolver |
| `dnn.tight_lam` | par key | integer | `2` | 0..10 | policy | the filter's upper end lam_hi from this many Lanczos steps (0: off) |
| `dnn.gap_ratio` | par key | number | `0` | 0..inf | policy | gap-ratio give-up rule in place of the extrapolation (0: off) |
| `dnn.gu_zone` | par key | number | `0.003` | 0..1 | policy | a give-up within this fraction of the cutoff is confirmed 100 iterations later (0: off) |
| `dnn.momentum` | par key | number | `0` | 0..10 | policy | fast ADMM with Goux's safeguards, restarted when the residual grows by this factor (1.05; 0: off) |
| `dnn.beta` | par key | number | `0.21` | 1e-06..100 | policy | ADMM penalty, as a multiple of dnn.py's m/3 |
| `dnn.refresh` | par key | integer | `120` | 1..1000000 | policy | full eigendecomposition every this many iterations |
| `dnn.warm_y` | par key | integer | `2` | 0..2 | policy | a child starts from its parent's Y too: 0 Z only, 1 restricted, 2 conditioned on its assignment |
| `dnn.beta_sched` | par key | text | `w=0.36` | text | policy | ADMM penalty schedule, e.g. two,f=1.6667,t=60,w=0.25,wc=0.7 ("": off) |
| `dnn.dscale` | par key | number | `0` | 0..4 | policy | diagonal congruence scaling of the lifted variable: exponent on the row norms of A and B (0: off) |
| `dnn.dscale0` | par key | number | `1` | 0..1000 | policy | with or without dnn.dscale: the scale of the homogenizing coordinate (1: unscaled) |
| `dnn.dscale_until` | par key | integer | `0` | 0..1000000 | policy | with dnn.dscale: the scaled iteration for this many iterations of a call, the plain one after (0: throughout) |
| `dnn.cross` | par key | integer | `2` | 0..100 | policy | extra certificates at the predicted cutoff crossing (0: off) |
| `dnn.cross_first` | par key | integer | `0` | 0..1000 | policy | with dnn.cross: one more extra certificate at this iteration of a warm call (0: off) |
| `dnn.lazycert` | par key | integer | `1` | 0..1 | policy | 1: Rump's test only at checks that can decide something |
| `dnn.amx` | par key | integer | `0` | 0..2 | policy | AMX bf16 kernels: 1 W*Y, 2 also CholQR and Rayleigh-Ritz (no effect without AMX) |
| `ship.bits` | par key | integer | `32` | 0..64 | policy | the parent's DNN state in state files: 0 none, 16 bf16, 32 fp32, 64 fp64 |
| `cut.zone` | par key | number | `0` | 0..1 | policy | run the cut phase when the uncut bound ends within this fraction of |incumbent| of the cutoff (0: off) |
| `cut.maxdepth` | par key | integer | `1073741824` | 0..inf | policy | the cut phase only at nodes of at most this depth |
| `cut.side` | par key | integer | `0` | 0..1 | policy | regions over the first (0) or the second (1) matrix |
| `cut.grid` | par key | text | unset | text | policy | R:C grid whose prefixes are the regions (unset: prefixes of the free elements) |
| `cut.regions` | par key | text | unset | text | policy | file of regions, one per line: original indices of cut.side (any sets; with cut.grid, added to its prefixes) |
| `cut.regions.hash` | par key | text | unset | text | policy | FNV-1a of the cut.regions file's bytes (set by bnbdrv; a resume with another file fails the -r check) |
| `cut.rounds` | par key | integer | `12` | 1..1000 | policy | cut rounds per phase |
| `cut.iters` | par key | integer | `300` | 1..1000000 | policy | ADMM iterations per cut round |
| `cut.k` | par key | integer | `1000` | 0..1000000 | policy | new triangle cuts per round |
| `cut.k5` | par key | integer | `0` | 0..1000000 | policy | 5-point count inequalities per round (0: off; rounds form only) |
| `cut.k7` | par key | integer | `0` | 0..1000000 | policy | local-7 inequalities per round (0: off; rounds form only) |
| `cut.sub7` | par key | integer | `1000` | 1..10000000 | policy | random 7-subsets sampled per region and round (rounds form only) |
| `cut.kbox` | par key | text | unset | text | policy | box inequalities per round, K[:QMAX] (unset: none; rounds form only) |
| `cut.delta` | par key | number | `50` | 0..inf | policy | Polyak target: best + delta (rounds form only) |
| `cut.stall` | par key | integer | `3` | 0..1000 | policy | after this many rounds, stop when the gap exceeds 3x the gain of the last that many (0: off) |
| `cut.pgate` | par key | text | unset | text | policy | primal gate R:M: stop after round R when the best primal is more than M below the cutoff (rounds form only) |
| `cut.nocap` | par key | integer | `1` | 0..1 | policy | 1: Polyak target uncapped at the cutoff; 0: capped (rounds form only) |
| `cut.freeze` | par key | integer | `600` | 0..1000000 | policy | iterations of the quiet round with frozen multipliers when a phase ends unfathomed with its primal above the cutoff (0: off) |
| `cut.warmy` | par key | integer | `0` | 0..1 | policy | 1: cut rounds restart from the previous round's Y as well as its Z |
| `cut.damp1` | par key | text | unset | text | policy | F[:G]: factor on the first (and later) Polyak steps of a phase that starts from a pool (rounds form only) |
| `cut.equip` | par key | text | unset | text | policy | G:R: a phase more than G below the cutoff runs at most R rounds, to equip the children |
| `cut.cp` | par key | text | unset | text | policy | k:K: the phase as one CP-ADMM solve separating K triangles every k iterations |
| `cut.giveup` | par key | text | unset | text | policy | T:G,...: from iteration T a CP-ADMM phase stops when its gap to the cutoff is still above G (unset: the stall rule only) |
| `cut.giveup_depth` | par key | integer | `0` | 0..inf | policy | cut.giveup applies at nodes of at least this depth (shallower phases stop by the stall rule only) |
| `cut.cp_iters` | par key | integer | `0` | 0..100000000 | policy | CP-ADMM iteration limit (0: rounds x iters + freeze) |
| `cut.handover` | par key | text | unset | text | policy | N0:A:B[:PHI[:W[:P[:KAPPA]]]]: a CP-ADMM phase stops when PHI x gap / rate, the iterations it still needs at the rate of its last W, exceeds C = m x (N0 + A g + B g^P), its children's cost at gap g in 0.1% of the cutoff, and (KAPPA > 0) rate x dC/dgap < KAPPA (unset: the stall rule) |
| `cut.cp_sweeps` | par key | text | `1:1` | text | policy | S[:P]: Dykstra sweeps and passes per CP-ADMM iteration |
| `cut.cp_extrap` | par key | text | unset | text | policy | E:GAMMA[:GUARD[:MODE]]: every E iterations a CP-ADMM phase extrapolates its state by GAMMA times its change over the last E/2, guarded (unset: off) |
| `cut.cp_wc_gap` | par key | number | `0` | 0..inf | policy | a CP-ADMM phase that starts within this of the cutoff keeps dnn.beta_sched's w instead of wc (0: off) |
| `cut.cp_wc_deficit` | par key | number | `0` | 0..1 | policy | a CP-ADMM phase takes wc only when the first solve's primal objective is at least this fraction of |cutoff| below its bound (0: off) |
| `cut.cp_lazycert` | par key | integer | `0` | 0..1 | policy | 1: the lazy certificate in CP-ADMM solves too |
| `cut.xbar` | par key | integer | `0` | 0..1 | policy | 1: rules 8/9 take their marginals from the cut phase's final iterate when a phase ran (0: the first solve's) |
| `pool.on` | par key | integer | `0` | 0..1 | policy | 1: a node hands its cut pool to its children |
| `dnn.retry` | par key | integer | `0` | 0..1 | policy | 1: a warm-started first solve that ends unfathomed more than 1% of the incumbent below the bound the node inherited is solved again, without its inherited pool and then cold (bound_lab/tai_solves; in the DNN trace as retry S rlb X) |
| `pool.cp` | par key | integer | `0` | 0..1 | policy | 1: a pooled child's first solve in CP form |
| `branch.weighted` | par key | integer | `1` | 0..1 | policy | rules 1, 3 and 4: weight the U sums by the row's gap fraction |
| `branch.weight` | par key | number | `0.5` | 0..1 | policy | the weight of that fraction |
| `branch.xbar_lines` | par key | integer | `1` | 1..64 | policy | rule 9: candidate rows and columns per side from the DNN marginals |
| `branch.la_exp` | par key | number | `0` | 0..500 | policy | rules 10 and 11: L, a line's score is the sum, over the children QPB does not drop, of exp(L x the child's QPB shortfall / the node's QPB gap); the rules need it set (32 in the runs measured); 0: unset |
| `sb.depths` | par key | text | unset | text | policy | DNN strong branching: dmin:dmax, the depths at which a branching node bounds the children of its candidate lines and branches on the line whose children predict the smallest tree (unset: off; needs dnn.depths over them and dnn.warm = 1) |
| `sb.k` | par key | integer | `4` | 1..1000 | policy | sb: the candidate rows, and the candidate columns where no symmetry is left, by marginal score |
| `sb.all` | par key | integer | `-1` | -1..inf | policy | sb: through this depth every free line is a candidate (-1: none) |
| `sb.iters` | par key | integer | `150` | 1..1000000 | policy | sb: warm ADMM iterations per child |
| `sb.power` | par key | text | `0:8,1:1` | text | policy | sb: P, or D:P[,D:P...]: from depth D on a line's score is the sum over its children of shortfall^P (1: the sum of shortfalls; the default: the 8th power at the root, the sum below) |
| `sb.gap` | par key | text | unset | text | policy | sb: D:G[,D:G...]: from depth D on only at nodes whose relative gap (inc - the node's DNN bound)/(inc - root QPB) is at least G (unset: every node in sb.depths) |
| `fw.patience` | par key or `-w` | integer | `0` | 0..100000000 | policy | -w: Frank-Wolfe patience at lookahead/strong levels (0: off) |
| `est.xbar` | par key | number | `0` | 0..inf | policy | estimator: child weights x_ij^a from the DNN marginals (0: QPB scores) |
| `est.xbar_mix` | par key | number | `0.1` | 0..1 | policy | estimator: uniform share mixed into the x_ij^a weights |
| `par.hash` | the `-p` file | text | unset | text | policy | FNV-1a of the params file's bytes |
| `DNN_R32` | environment | integer | `0` | 0..1 | experiment | 1: R's rank term in fp32 on certificate iterations too |
| `DNN_FUSE` | environment | integer | `0` | 0..1 | experiment | 1: fuse the Y/Z step with the next forward's row phase |
| `DNN_ADMM_WARM_BETA` | environment | text | unset | text | experiment | F1:T:F2[:force]: a warm call's own penalty schedule |
| `DNN_ADMM_WARM_DAMP` | environment | text | unset | text | experiment | D:K: both dual steps scaled by D for the first K iterations of a warm call |
| `QAP_EST_SEED` | environment | text | unset | text | run | estimator: seed for rand() (unset: the fixed sequence) |
| `QAP_EST_SURV` | environment | text | unset | text | run | estimator: b, survivor dives with weights exp(b x); kQ, Knuth's dives inside each start node, g^Q of the child's QPB bound (estimate.C; unset: Knuth's dives from the root) |
| `QAP_EST_SAMPLE` | environment | number | `1` | 0..1 | run | estimator, QAP_EST_SURV: the share of the start nodes sampled as strata (weighted 1/pi) |
| `QAP_EST_ALLOC` | environment | text | `auto` | text | run | estimator, QAP_EST_SURV: neyman (dives where the variance per CPU s is; reads low with Knuth's heavy tails), equal (round the strata in turn), auto (equal for kQ, neyman for survivor dives) |
| `QAP_EST_SAMPLE_W` | environment | text | `unif` | text | run | estimator, QAP_EST_SAMPLE: unif, or qpb (each start node weighted g^2 of its own QPB bound) |
| `QAP_EST_BUDGET` | environment | number | `0` | 0..inf | run | estimator, survivor dives: stop starting dives after this many wall seconds, counted from the end of the strata's setup (0: estimate.par's trials) |
| `QAP_EST_MIN_DIVES` | environment | integer | `1` | 1..1000000 | run | estimator, survivor dives: every stratum with children that branch gets this many dives whatever the budget (then QAP_EST_BUDGET) |
| `QAP_EST_SAMPLE_DEPTHS` | environment | text | unset | text | run | estimator, survivor dives: D:F,D:F,... a start node at depth D is a stratum with probability F, whatever QAP_EST_SAMPLE (unset: QAP_EST_SAMPLE for every depth) |
| `QAP_EST_WARM_KEEP` | environment | integer | `8` | 0..1000000 | run | estimator: the strata whose children keep their parent's DNN state in memory (the rest spill it to a file, QAP_EST_SPILL_MB) |
| `QAP_EST_SPILL_MB` | environment | number | `1000` | 0..1e+09 | run | estimator: MB of spilled strata states a process may write beside QAP_EST_DIVES (past it a stratum's children start cold, without the parent's pool) |
| `QAP_EST_DIVES` | environment | text | unset | text | run | estimator, survivor dives: append a JSON line per stratum and per dive to this file |
| `QAP_EST_DNN_TRACE` | environment | integer | `0` | 0..1 | run | dives runs: 1, the coordinator's worker keeps the DNN trace of every node a dives unit bounds (coord/worker.py sets QAP_DNN_TRACE for the unit and drops this one); the solver ignores it. Registered so that a coordinating instance whose environment has it (a spec's ESTIMATE_ENV) can still run the solver: sko42-est2-261121's gates failed on it, Nov 21 |
| `QAP_SB_PROGRESS` | environment | text | unset | text | run | DNN strong branching: append every child solve to this file and reuse the ones already in it (a node whose process stopped mid-way repeats only what is missing) |
| `QAP_DNN_WARM_DIR` | environment | text | unset | text | run | parent-state store: write each parent's DNN state to this directory and reference it from state files (unset: inline) |
| `QAP_ALLOW_UNKNOWN` | environment | integer | `0` | 0..1 | run | 1: an unknown or pre-Phase-2 QAP_DNN_/QAP_EST_/QAP_XBAR_/QAP_NODE_/DNN_ variable is a warning, not an error |
| `par.file` | `-p` | text | unset | text | run | -p: the params file |
| `QAP_DNN_TRACE` | environment | text | unset | text | diagnostic | append a line per DNN node to this file |
| `QAP_DNN_TRACE_HIST` | environment | integer | `0` | 0..1 | diagnostic | 1: the trace line carries the CP cut phase's running best at each check (| hist it:best ...) |
| `QAP_NODE_TRACE` | environment | text | unset | text | diagnostic | append a line per node taken off the stack (fixes, bounds, how it ended, the branching line; bnb.C) to this file |
| `QAP_SB_TRACE` | environment | text | unset | text | diagnostic | append a line per strong-branching node (the line chosen, every candidate's score, solves, iterations, CPU; bnb.C) to this file |
| `QAP_DNN_POOL_CHECK` | environment | integer | `0` | 0..1000000 | diagnostic | check each assembled pool term on this many random permutations |
| `QAP_EST_HIST` | environment | text | `dive_history_nug25_2.txt` | text | diagnostic | estimator: the per-dive history file |
| `QAP_EST_CHILD_DNN` | environment | integer | `0` | 0..100000 | run | estimator, Knuth's dives: weight each child a dive may take by its own short DNN solve of this many ADMM iterations, warm from its parent's state (experiments/estimates_nov19; 0: by QPB) |
| `QAP_EST_CHILD_QPB` | environment | integer | `0` | 0..1 | diagnostic | estimator: 1 weights each child a dive may take (and the first draw) by its own QPB bound, not the parent's prediction |
| `QAP_XBAR_DEBUG` | environment | integer | `0` | 0..2 | diagnostic | 1: print the rule 8/9 line choice per node to stderr; 2: also every free line ranked by its marginal score |
<!-- settings-table-end -->

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
state and stops it at the next node. Since version 8 the file also records
the settings the run resolved ("Settings" above), and `-r` refuses a file
written under other settings unless `--allow-settings-change` is given.

With the DNN bound (`-x`), a branched node's final ADMM state (Z and Y), from
which its children start, is held in memory only. It goes into the state file
(version 6) only with `ship.bits = 16|32|64` in the params file (bf16, fp32
or the exact doubles). With 64, a stopped and resumed DNN search is the uninterrupted
search node for node (`make shipcheck`): 64 ships the whole matrices, since
with `DNN_CORE_NOSYMV` (the default) Z and Y are symmetric to rounding only and
a mirrored upper triangle is not the parent's state; 16 and 32 ship the upper
triangles. Without it, nodes restored from a
file start their ADMM cold, which costs iterations but nothing else. This
matters for `bench/parallel.py`, whose tasks are state files: on nug30,
`ship.bits = 32` took 15% off the CPU on an identical tree.

**The parent-state store** (state version 9, `docs/parent_state_store.md`).
At m = 40 a shipped state is 10 MB of text, and every file that mentions the
parent carries it. With `QAP_DNN_WARM_DIR=DIR` in the environment the solver
writes each parent's state once to `DIR/<hash>.dnn` and the warm line says
`@<path>` instead (relative to the state file's directory when the store is
under it, absolute otherwise); a reader takes the token from that file, and
a missing file is an error, never a silent cold start. A state file and its
store directory travel together: copy `warm/` along with a checkpoint.
Without the variable the file is self-contained, as before. `make shipcheck`
checks that a stop and resume through the store is the same search node for
node.

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

The workers write their parents' DNN states to a parent-state store
(`--warm-dir`, default `WORKDIR/warm`; `--warm-dir ''` for the inline tokens
of before Nov 9 2026), so a task file carries a reference per parent instead
of 10 MB of text, and the coordinator never parses a token. A checkpoint is
written with references into its own `warm/` directory next to it (the
referenced files are hard-linked or copied there), so `F` and `warm/` are the
unit to keep or move; `--resume F` reads them from there.
`--checkpoint-inline` writes a self-contained checkpoint instead, for a
hand-off to another machine. At every checkpoint the store files no open
node references, and no running worker can still be writing, are deleted.

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
`cut.zone`): nodes, fraction fathomed, CPU. `Fathoms by cause` splits
each depth's fathoms: by QPB, by the DNN bound, by the cut phase, and nodes
popped after the incumbent improved past their bound (`qpb` is what is left
after the other three; a run resumed from a state file before version 7 has no
causes for its first part). Runs without `-x` print exactly what they did.

With DNN strong branching (`sb.depths`, below) the head has one more line and
a fourth table follows, both only where it chose a line (state version 10,
written only then: other state files stay version 9):

```
 strong branch: 7 nodes, 448 child solves, 5220 iterations, 12.2 s

Strong branching by level:
level    nodes	  solves	      iters	    time
 1           7	     448	       5220	   12.20
```

by the depth of the node that branched: the nodes whose line it chose, the
child solves it ran, their ADMM iterations and CPU seconds. That time is in
the level table's and the total, and in none of the DNN tables.

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
4     branching rule: 1 U-sum, 2 U-fathom, 3 strong, 4 lookahead, 5 learned (BUFF), 6 eccentric, 7 eccentric lookahead, 8 DNN marginals, 9 lookahead between the two marginal lines (8 and 9 need -x; lookahead where no DNN bound ran), 10 lookahead with the size score, 11 rule 9's lines with the size score (10 and 11 need branch.la_exp)
150   Frank-Wolfe iterations for the bound
100   Frank-Wolfe iterations when we know we can't fathom
25    Frank-Wolfe iterations for row selection (lookahead children; each starts from the parent's X)
1     bound: 0 QP, 1 QP parametric, 2 QP+GLB, 3 QP improved, 4 GLB, 5 PB, 6 EVB3, 7 BUFF
30    step for parametric bound
30    compute B for the best _ rows of U
```

After the blocks an optional `[settings]` section holds the keyed settings
("Settings" above), one `key = value` per line, `#` starting a comment:

```
[settings]
# the DNN bound from depth 1 through depth 4, with the cut phase in CP form
dnn.depths = 1:4
dnn.giveup = 300
cut.cp = 100:1000
```

The 53 files that have no section parse unchanged. An unknown key, or a
value out of range, is an error. `experiments/wil50_feasibility/fullrun/fullrun.par`
is a complete example; `bench/envpar.py base.par --env OLD=VALUE ... -o new.par`
writes one from a pre-Nov-9 environment recipe.

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
- **Rule 10** is rule 4's lookahead with another score (Nov 23;
  `experiments/branching_lab/ucheap` on `claude/cool-babbage-ilnogj`). The
  candidates are the level's best "compute B for best _ rows of U" rows and
  columns by U sum, every free line when that number is m or more (rows only
  while symmetry is left). Every child of a candidate is bounded by QPB as
  rule 4 bounds it, with the level's "Frank-Wolfe iterations, for row
  selection". A child's shortfall s is the cutoff less that bound. The node
  branches on the line with the least sum, over the children QPB does not
  drop, of exp(L s / G), G the node's own QPB gap and L = `branch.la_exp` (a
  child QPB drops adds nothing): one child far from the cutoff outweighs
  several near it, which rule 4's sum of the children's bounds cannot see. It
  needs no DNN bound.
  - `branch.la_exp` has no default: a params file that uses rule 10 must set
    it, or the run stops at startup. It is in the settings hash only when
    set, so recipes without rule 10 keep their hash. 32 is the value of every
    run measured.
  - Measured with every line a candidate at every lookahead level of the
    nug30 recipe: nug30 3,420 nodes and 2,392 CPU s (rule 9: 4,285 to 4,365 and
    3,314 to 3,421), tho30 8,581 and 6,350 (14,815 and 9,375); the lookahead
    over every line is under 2% of that CPU. Under DNN strong branching at
    depth 1 (`sb.depths = 1:1`): 3,140 and 7,757 nodes. On three depth-3
    subtrees of tho40 every line gains nothing (1,862 nodes against rule 9's
    1,818, one subtree 38% larger): see rule 11. Not measured: a shortlist
    shorter than every line (by U sum), and any tai instance.
  - `PARAM/rule10.par` is the rule alone, for the brute-force check.
- **Rule 11** is rule 9 with rule 10's score: the candidates are the
  `branch.xbar_lines` rows and columns with the largest DNN marginal scores
  (rule 4's candidates where no DNN bound ran or symmetry is left, as for rule
  9), chosen by the exponential score on their children's QPB bounds. It needs
  `branch.la_exp` as rule 10 does.
  - Which of the two depends on the instance, as the width of rule 9's
    shortlist did (above): nug30 3,261 nodes under rule 11 against 3,420 under
    rule 10 with every line and 4,285 to 4,365 under rule 9; three depth-3
    subtrees of tho40 (`branch.xbar_lines = 2`) 1,607 nodes under rule 11,
    1,862 under rule 10 with every line, 1,818 under rule 9; tho30 14,303
    under rule 11 (one line a side: the marginal lines hold nothing there),
    8,581 under rule 10 with every line, 14,815 under rule 9.

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
| tho30, the same with the CP form (`cut.cp = 100:1000`; Nov 8) | 17,786 (DNN: 17,779, 17,056 fathomed, 5,957 s) | 6,994 | 17.3 min (experiments/tho30_solve) |
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
solver's default since Nov 8 2026**: with no `[settings]` section,

```sh
python3 bench/parallel.py qaplib/nug30.dat --jobs 8 --ub 6124 \
    --params PARAM/nug_deep_xbar9_all.par --extra "-x 1:7:1:200:1" \
    --checkpoint nug30.ckpt --checkpoint-every 300
```

gives the record tree (3,716 nodes with rule 9, verified node for node
against the explicit recipe). A `[settings]` key still overrides every
default, and writing the recipe out is harmless. The defaults that changed
on Nov 8, old -> new, for anyone whose A/B control was "unset" (the keys are
Nov 9's names; until then each was an environment variable, `QAP_DNN_BETA`
for `dnn.beta` and so on, `bnbdrv --settings-map` lists them):

| setting | old default | new default | setting | old default | new default |
|---|---|---|---|---|---|
| `dnn.gamma1` | 0 | 1.0 | `DNN_CHOLQR_PASSES` | 2 | 1, then removed Nov 9 |
| `dnn.gamma` | 1.6 | 1.0 | `DNN_R32_NONCERT` | 0 | 1, then removed Nov 9 |
| `dnn.cheb` | 2 | 3 | `dnn.cross` | 0 | 2 |
| `dnn.tight_lam` | 0 | 2 | `dnn.lazycert` | 0 | 1 |
| `dnn.beta` | 0.3 | 0.21 | `DNN_CORE_NOSYMV` | 0 | 1, then removed Nov 9 |
| `dnn.refresh` | 20 | 120 | `DNN_NOSYMW` | 0 | 1, then removed Nov 9 |
| `dnn.warm_y` | 0 | 2 | `DNN_CORE_REFRESH_PARTIAL` | 0 | 1, then removed Nov 9 |
| `dnn.ext` | 0 | 20 | `dnn.beta_sched` | unset (off) | `w=0.36` (`""` turns it off) |
| `dnn.gu_zone` | 0 | 0.003 | `ship.bits` | 0 | 32 |

"Removed Nov 9": the variable and its losing code path are gone
(`docs/settings_plan.md`, step 1), so the winning behaviour is unconditional
and the name is an error. `DNN_GUARD_DIV` (20) and the `DNN_CERTX*`
experiment went the same way. The `dnn.*` keys apply to the search and, with
`dnntest -p`, to `dnntest`, whose cold solves are unaffected by
`dnn.beta_sched = w=...` (a warm-child factor). The cut phase (`cut.*`,
`pool.cp`, `cut.cp_lazycert`) stays opt-in. The previous form of this block,
for the record (an environment that the solver now rejects;
`bench/envpar.py` converts it):

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

The DNN settings, as `[settings]` keys of the params file (the table above
has every one with its default; each has its experiment report; the values
here were the opt-in ones before Nov 8): `dnn.gamma1 = 1.0 dnn.gamma = 1.0`
the Peaceman-Rachford step and `dnn.cheb = 4` the filter degree
(experiments/eigensolver: 1.5x less DNN CPU on nug30); `dnn.tight_lam = 2`
with `dnn.cheb = 3` (degree-4 strength at 85% of the cost);
`dnn.gap_ratio = 0.65` the gap-ratio give-up rule
(experiments/wil50_feasibility/tuning); `dnn.ext = 20` keeps the DNN bound
up to 20 levels below a last-level miss instead of handing its subtree to
QPB alone, and `dnn.gu_zone = 0.003` makes a give-up within 0.3% of the
cutoff wait 100 iterations for confirmation (both from J.-P. Goux's nug30
solver; experiments/wil50_feasibility/tuning/goux). `cut.zone = 0.002` adds
a cut phase after a DNN miss that ends within that fraction of the cutoff:
regional triangle cuts (latest/dnncut.h) from the node's own iterate,
`cut.grid = R:C` the grid whose row and column prefixes are the regions
(wil50 `5:10`; nug30 `5:6` with `cut.side = 1`, the second matrix),
`cut.rounds` (12), `cut.iters` (300), `cut.k` (1000), `cut.delta` (50),
`cut.k5` (0: also this many of the regional 5-point count inequalities a
round, latest/dnncut.h); every value stays certified
(experiments/wil50_feasibility/stronger_bounds). On nug30 the cut phase
costs more than it saves (+34% CPU on the laptop, +42% on the cloud), so the
record setting leaves it off.

Later switches, all used in the nug30 record setting except the last two
items (defaults in parentheses; experiment reports in
`experiments/wil50_feasibility/`):
- `dnn.beta = F`: the ADMM penalty is F·m/3 (0.3).
- `dnn.refresh = N`: a full eigendecomposition every N iterations (20).
  Goux's 0.21 and 120 together took 20% off nug30 (`tuning/goux`); on wil50
  they do not pay.
- `dnn.warm_y = 2`: a child starts from its parent's final Y, conditioned
  on its new assignment, as well as its Z (0: Z only; `warm_y/`). With it, a
  child's `-x ...:W` warm start carries the primal too.
- One CholQR pass in the subspace step (was `DNN_CHOLQR_PASSES`, default 2;
  `roadmap/node_speed.md`). Unconditional since Nov 9.
- `ship.bits = 32`: carry the parent's ADMM state in state files, as
  above (0).
- `dnn.amx = 1`: Intel AMX bf16 kernels for the subspace products, on CPUs
  with AMX-BF16 only (`amx/`).
- The wil50 cut-phase switches `pool.on`, `cut.warmy` and `cut.damp1`
  (`pool_search/`, `pipeline_ab/`). Unused on nug30.
- `dnn.beta_sched = two,f=1.6667,t=60,w=0.7` (unset): a schedule for the ADMM
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
  `dnn.beta_sched = two,f=1.6667,t=60,w=0.25,wc=0.7`. `w=0.25` for the cut
  rounds as well (`wlow`) loses: the rounds take 1.5-1.9x the iterations.
- `dnn.cross = 2` (unset): between two regular certificates (every 10
  iterations), certify once more at the iteration where the running best is
  predicted to cross the cutoff (`speed/nov6_admm/`): fathoms come 0-9
  iterations earlier (4.2 per fathom on 24 warm children, 1.05x their CPU),
  nothing else changes. In the search: 1.035x in DNN iterations and 1.04x in
  DNN CPU (`pipeline_ab/`, Nov 6). `DNN_ADMM_CROSS_FIRST` and
  `DNN_ADMM_WARM_DAMP` are the same study's other probes (measured not worth
  it); all off by default.
- `DNN_ADMM_WARM_BETA=F1:T:F2` (unset): a node's own warm solve runs at F1 x
  the base penalty until the first full eigendecomposition at or after
  iteration T, then at F2 x it, in place of `dnn.beta_sched`'s `w`; cut rounds
  keep `w` or `wc`. Under test Nov 6 as `0.25:60:0.7` (`pipeline_ab/` arm
  `wtwo`): children fathom as at 0.25, and a miss hands its cut phase a state
  closer to 0.7's.
- R's rank term in fp32 except on certificate iterations (was
  `DNN_R32_NONCERT=1`; unconditional since Nov 9; `speed/r32_noncert/`). 4.2% less per wil50 iteration without
  AMX, 2.5% on nug30; about 1% with `dnn.amx = 1` (the `fast` arm above).
- `dnn.lazycert = 1` (unset): Rump's Cholesky test of the certificate only
  at checks that can decide something (the first, any whose polyhedral value
  could clear the cutoff, the last iteration, and the stopping check and the
  stored best at the exit); the same fathom decisions and returned lb by
  construction (`speed/nov6_core/`).
- The subspace tracker's full eigendecomposition back-transforms only the
  columns it reads (was `DNN_CORE_REFRESH_PARTIAL=1`; unconditional since
  Nov 9; bit-identical eigenvectors, about a third of the refresh saved;
  `speed/nov6_core/`).
- No symmetrisation pass after the expand and forward transforms (was
  `DNN_CORE_NOSYMV=1` with `DNN_NOSYMW=1`; unconditional since Nov 9); the
  iterate is symmetric to rounding only and the certificate evaluates the
  symmetric part exactly (`speed/nov6_core/`).
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

### DNN strong branching (sb.*)

`sb.depths = D1:D2` in the params file's `[settings]` turns it on at branching
nodes of those depths (Nov 22; `experiments/branching_lab/root` and
`experiments/branching_lab/port` for the measurements). It chooses the line in
place of the level's rule; everything else about the node is unchanged. Off by
default: without `sb.depths` the search is node for node what it was.

- **What it does at a node.** The candidates are the `sb.k` rows and the `sb.k`
  columns with the largest marginal scores (rule 9's score), or every free line
  through depth `sb.all`. With symmetry left only rows are candidates and a row's
  children are the columns the node would create. Every child of a candidate that
  QPB does not already drop is bounded with `sb.iters` ADMM iterations, started
  as its own first solve would be (warm from the node, under the node's cut pool
  where it hands one on), with no give-up and no cut phase; a child's bound is
  never taken below the one it inherits (QPB's for it, or the node's own DNN
  bound, which holds for every child). A child solve fails when it ends
  unfathomed below QPB's bound for the child, or with no certificate. (Not the
  1% test of a node's own solve: a child solve stops at `sb.iters`, and one that
  is only short of converged ends above QPB's bound.) With `dnn.retry = 1` a
  failed child is solved again as a node's own would be: without the pool, then
  cold, each with `sb.iters`, and the best certified bound counts. The trace
  counts the children retried and those still failed, which score as the node
  itself does. A line's score is the
  sum over its children of shortfall^P, the shortfall being the cutoff less the
  child's bound (0 when it fathoms); the node branches on the line with the
  least, ties to the larger marginal score.
- **The exponent, by depth.** `sb.power = D:P[,D:P...]`: from depth D on, P (a
  bare P: every depth). The
  default `0:8,1:1` is the 8th power at the root and the plain sum below. At a
  root every child is far from fathoming and a row is as bad as its worst
  children: the sum's first choice was the worst tree run at nug30 and at tho30;
  the 8th power's is the best root row run at nug30, production's at tho30, and
  at tho40 row 2 for every power from 2 to 12. Below the root many children
  fathom or nearly and their number counts: on nug30 the 8th power at depths 1
  and 2 gave 3,374 nodes where the sum gave 3,020.
- **Where it needs the DNN bound.** A node needs its own DNN state: `dnn.depths`
  must cover `sb.depths` with `dnn.warm = 1`, or the run stops at startup. For the
  root that means `dnn.depths` from 0 (the grids' recipes start at 1).
- **The gap gate.** `sb.gap = D:G[,D:G...]`: from depth D on, only at nodes whose
  relative gap (incumbent - the node's last DNN bound) / (incumbent - root QPB
  bound) is at least G. The same quantity as the params blocks' "use this rule
  when gap is worse than", with the DNN bound in place of QPB's and tested after
  the node is bounded. For example `sb.depths = 0:2`, `sb.all = 0`,
  `sb.gap = 2:0.06`: every row at the root, the shortlist at every depth-1 node,
  and at the depth-2 nodes at least 6% of the root gap short.
- **What it costs, and where it is counted.** Its CPU is the run's: in the total
  time, in the level table's time at the node's depth, and on its own in the
  RESULTS line `strong branch: N nodes, S child solves, I iterations, T s` and
  the table "Strong branching by level" (state files of version 10, written only
  once strong branching has chosen a line: until then, and in every run without
  it, state files are version 9 as before). A child solve, a retry included, is
  not a node. One node costs (2 x `sb.k` x m - `sb.k`^2) x `sb.iters`
  iterations at dimension m - 1 at most (`experiments/branching_lab/budget` for
  the open instances: about 1% of a run at depth 1, more than the run at depth 3).
- **`QAP_SB_TRACE=file`** (diagnostic): a line per strong-branching node:
  `sb depth D fixes F gap G | chose LINE score S open K | N solves, R reused, Y retried, Z failed, I iterations, T CPU s, largest shortfall X | LINE:marginal:score:open ...`.
- **`QAP_SB_PROGRESS=file`** (run): every child solve is appended to the file and a
  child already in it (the same settings hash, root bound, node, cell, `sb.iters`
  and cutoff) is not solved again, so a node whose process stopped mid-way repeats
  only what is missing. Within a process the child solves are kept in memory
  whether or not the file is set, so the estimator's dives pay for a node once. The
  coordinator does not manage this file: it helps where the file outlives the
  process (a laptop, an instance that restarts), not across a lost instance. The
  root's own pass is kept by the coordinator's root cache as before.
- **The settings hash.** The `sb.*` keys are in the hash only when `sb.depths` is
  set (then all of them, defaults included). A run without it has the hash it had
  before these keys existed, so its state files resume as before.
- **The marginals.** The shortlist's marginal scores are the node's own copy:
  `sb.depths` does not change what rules 8 and 9 or the estimator's `est.xbar`
  see at a node where strong branching does not choose the line.
- `make sb` (`tests/sb_test.py`) checks the statistics, the trace, the progress
  file, `bench/parallel.py`'s totals, the gate, and that a state file is version
  10 only once strong branching ran.

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
