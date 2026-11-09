# Settings cleanup: one home for every solver setting (plan, Nov 9 2026)

This plan is for the laptop session, alongside the Nov 8 review's cleanup
(`docs/review_nov8.md`). Both phases land before the wil50 run. The rule they
implement is in `CONTRIBUTING.md`, "Where settings live".

## Why

The solver takes settings through four channels, and nothing records which
ones a run used.

| channel | where it's read | what's in it |
|---|---|---|
| command line | `driver.C`, about 20 flags | the run: instance, `-u`, `-s`, `-p`, `-r`/`-o`, `-k`/`-t`, `-F`, `-d`; and two policy flags, `-x` (DNN depths, gap, give-up, warm start) and `-w` (FW patience) |
| `.par` file | `util.C:107` `qap_read_params` | tree policy, by blocks keyed on relative gap and depth: rule, FW iterations, which bound, step, NBEST. Positional, no keys; 53 files |
| environment | about 72 `getenv` names in `bnb.C`, `dnn.C`, `dnncut.C`, `state.C`, `branch.C`, `estimate.C` | the DNN bound, the cut phase, state shipping, core switches, estimator, diagnostics |
| compile time | `util.h`, inline constants | weighted branching on with weight 0.5 (`util.h:44-45`); rule 9's two candidate lines (`branch.C`, `kr = kc = 1`) |

The problems:

- **Provenance.** The solver does not print the environment it ran with, and
  the state file does not record it. A resume with a different environment
  behaves differently without a word. A mistyped name is ignored. The
  launchers work around this with `env.json`/`env.txt` beside each run.
- **Defaults move under the launchers.** 34e347c (Nov 8) made the nug30
  record recipe the default. `fullrun/settings.sh` sets 16 of the 18 changed
  variables, but not `QAP_DNN_BETA` (0.3 -> 0.21) or `QAP_DNN_GU_ZONE`
  (0 -> 0.003). The gated `cpcutt` arm and stages A and B ran the old values.
  A binary built from cognac after 34e347c would have run the full run with
  two ungated settings, and nothing would have shown it. Step 0 pins them.
- **Policy is split across channels.** The par file's gap and depth blocks
  exist to vary policy over the tree, but the DNN's depth range is on the
  command line (`-x`), the cut phase's depths are in the environment
  (`QAP_DNN_CUT_MAXDEPTH`), and `dnntest` has its own flags (`-B`, `-L`, `-g`,
  `-G`, `-H`, `-d`) for the same DNN parameters `bnbdrv` reads from the
  environment.
- **Duplicates and leaks.** The sym file and the depth limit can be set in the
  par file and by flag (and `util.C:126-130` lets the par field overwrite the
  `-s` pointer, review D). Rule 9 ignores the par file's NBEST. Review D's
  "env switches read inconsistently": some are latched in function-local
  statics, some re-read every call.

## Where it runs

- **The laptop does all of the work.** Phases 1 and 2 touch `latest/`, the
  launchers and the docs, in the same files as the Nov 8 cleanup.
- **The cloud box doesn't build or edit `latest/` until stage B finishes** (about
  10:30-12:30 UTC Nov 9). Stage B runs a frozen copy of the stage-A binary
  (a667727), so nothing here affects it. After stage B, the cloud box runs the
  wil50 depth-2 replays of the Phase 2 gate, which need the 7.2 GB of stage-A
  parent states that live only there.
- **Equivalence runs compare old and new binaries on the same machine.** The
  laptop's Accelerate and the cloud box's OpenBLAS do not give bitwise-equal
  floats, so a cross-machine diff proves nothing.

## Order

0. **Done Nov 9**: `fullrun/settings.sh` pins `QAP_DNN_BETA=0.3` and
   `QAP_DNN_GU_ZONE=0`, the values the gated arm and stages A and B ran.
1. **Review D's retirements first.** Delete the CERTX experiment and the
   measured-out switches before building the registry, so dead knobs don't get
   migrated.
2. **Phase 1, provenance.**
3. **Phase 2, one home.**
4. **Re-gate the full-run configuration** (the last section).

## Phase 1: provenance (about a day)

**Status, Nov 9 (laptop)**: done, with three deviations. (a) The unknown-name
check covers the solver's own prefixes `QAP_DNN_`, `QAP_EST_`, `QAP_XBAR_` and
`DNN_`, not every `QAP_*` name: the launchers use `QAP_ROLE`, `QAP_BUCKET`,
`QAP_BACKEND`, `QAP_RUN` and a dozen more for themselves, and those reach the
solver's environment. (b) The golden dump (`fullrun/settings.golden`) is taken
from the environment alone; the launcher's `-p`/`-x` show as defaults in it and
a launcher that passes them diffs its own dump. (c) `sampling/stage_a.py`
writes the dump as `settings.txt` beside `env.json` rather than in its place,
since other tools read `env.json`. The hash covers the `experiment` class as
well as `policy` (both change the search); `par.hash`, the params file's
content hash, is a policy entry, so a different par file on `-r` is refused
too. Gate: `make test`, the regression check at 0%, nug30 3,716 nodes, and the
dump under `settings.sh` equal to the dump rebuilt from the cpcutt arm's
`env.txt` plus the two pins.

1. **A settings registry.** Add one table, in a new `latest/settings.{h,C}`.
   Each entry has:
   - the name
   - the type and parse rule
   - the default
   - the valid range
   - a one-line description
   - a class: `policy` (changes the search), `run` (this run's identity or I/O),
     `experiment` (an A/B in flight) or `diagnostic` (output only)

   Every `getenv` in `latest/` reads through it, once at startup, into one
   settings struct that is passed down. That also settles review D's
   latched-versus-re-read item. `dnntest` uses the same registry for the DNN
   entries.
2. **A startup dump.** `bnbdrv` and `dnntest` print a `settings:` block at
   startup. It lists every registry entry with its value and its source:
   `default`, `par`, `env` or `flag`.
   - `bnbdrv --print-settings [-p file.par] [flags]` prints the block and exits.
     It needs no instance.
   - The block is the oracle for every equivalence check below. It is far
     cheaper than tree runs and catches a moved default directly.
3. **Unknown names are errors.** A `QAP_*` or `DNN_*` variable in the environment
   that the registry doesn't know stops the run (exit 2), naming the variable.
   `QAP_ALLOW_UNKNOWN=1` overrides this for bisecting old commits. An error, not
   a warning: a typo on 1,000 cores has to stop the run, not cost a week.
4. **Settings in the state file (QAP-STATE 8).**
   - Record a hash of the resolved `policy` entries, plus the entries that differ
     from their defaults, as text.
   - On `-r`, compare them. On a mismatch, print both sets and exit, unless
     `--allow-settings-change` is given.
   - Task files carry it, so `parallel.py` workers and the coordinator's
     workers are checked without changes of their own. Version 7 files still
     read; they just have no record.
5. **Launchers check a golden dump.** `fullrun_check` greps the binary for
   names. Replace that with a diff: the full-run configuration's
   `--print-settings` output against a committed `fullrun/settings.golden`.
   - Any difference, including a changed default, fails the check.
   - The sampling launcher (`sampling/common.py`) stores the dump in each run
     directory in place of `env.json`.
6. **Review C5 in the same pass.** Fix the option parsing:
   - a flag with no argument is an error
   - `atoi`/`atof` garbage is an error (`-u 3O00` must not become 0)
   - an unknown option is an error
   - `-F` with `-r`, and `-t` without `-o`, are errors

   It is the same failure class: a run silently not doing what was asked.

**Gate (bit for bit; Phase 1 must not change any search):**
- `make test` passes.
- `bench.py` standard suite gives a node ratio of exactly 1.000.
- nug30 with rule 9 gives its record tree of 3,716 nodes.
- `--print-settings` under `settings.sh` equals the dump rebuilt from
  `pipeline_ab/runs_cpt/cpcutt/env.txt` plus `QAP_DNN_BETA=0.3` and
  `QAP_DNN_GU_ZONE=0`. This proves the registry's defaults and parse rules.

## Phase 2: one home (2-3 days, most of it launchers and equivalence runs)

**Status, Nov 9 (laptop)**: done, with these deviations and notes.
- **Syntax.** `[settings]` is one section after the positional blocks,
  `key = value` per line, `#` comments; `[block k]` is not implemented (no
  per-block policy exists yet; the reader rejects such a block as "not
  supported yet"). The 53 existing files parse unchanged.
- **The flags stay.** `-x` and `-w` remain as overrides of their keys (source
  `flag`); an empty `-x` field now leaves that key to the file or its
  default instead of resetting it. The par file's positional sym and
  max-depth fields are left as they were (item 3's warning is not done).
- **Keys** (`bnbdrv --settings-map` lists the old names): `dnn.depths`,
  `dnn.maxgap`, `dnn.giveup`, `dnn.warm`, `dnn.ext`, `dnn.gamma1`, `dnn.gamma`,
  `dnn.cheb`, `dnn.tight_lam`, `dnn.gap_ratio`, `dnn.gu_zone`, `dnn.momentum`,
  `dnn.beta`, `dnn.refresh`, `dnn.warm_y`, `dnn.beta_sched`, `dnn.cross`,
  `dnn.lazycert`, `dnn.amx`, `ship.bits`, 25 `cut.*`, `pool.on`, `pool.cp`,
  `branch.weighted`, `branch.weight`, `branch.xbar_lines`, `fw.patience`,
  `est.xbar`, `est.xbar_mix`, `par.hash`. The `DNN_*` experiment switches
  that the production arms use (`DNN_BETA_SCHED`, `DNN_ADMM_CROSS`,
  `DNN_CORE_LAZYCERT`, `DNN_AMX`) graduated to keys; `DNN_R32`,
  `DNN_ADMM_CROSS_FIRST`, `DNN_FUSE`, `DNN_ADMM_WARM_BETA` and
  `DNN_ADMM_WARM_DAMP` stay in the environment as experiments.
- **The old name is an error** whose message names the key, not a silent
  alias, so no launcher can keep exporting it unnoticed. `bench/envpar.py`
  converts a recipe (and `--from-env --shell` converts a bash launcher's
  exports in place); it reads the map from the binary.
- **dnntest** takes `-p file.par`; its own flags override the file, its
  DNNParams defaults are the registry's defaults in that process
  (`settings_override_default`), and the dump shows which.
- **The table is generated**: `bnbdrv --settings-table`, included in
  `RUNNING.md` between markers, checked by `tests/settings_doc.py` in
  `make test` (`--update` rewrites it). `--settings-table` and
  `--settings-map` print nothing else and do not read the environment.
- **Launchers converted**: `tests/parallel_test.py` and `tests/ship_test.py`
  (`tests/par/`), `infra/phase1/recipes/nug30.sh` (`recipes/nug30.par`,
  `EXTRA` and `ENV` empty), `infra/phase0/bench.sh` (`dnntest.par`),
  `fullrun/settings.sh` + `fullrun.par` (item 7; `QAP_DNN_POOL_CHECK=20`, a
  diagnostic, stays exported so the configuration is the measured one to the
  letter), the sampling launchers (`common.PARAMS` is `fullrun.par`; each
  stage-A run writes its own `params.par` from the base plus `--env`
  settings; stage B derives `giveup.par` and `converge.par`; `validate.py`
  its `params.par`), `pipeline_ab/run_pair.sh` (each arm's exports become
  `runs/ARM/arm.par` through envpar; the certx arm no longer runs),
  `tail_probe/tail_probe.py` (the arm dicts, through envpar, to `arm.par`),
  `coord` (`--env` is for experiment switches and diagnostics; policy is the
  `--params` file). Historical experiment scripts are untouched.
- **Gate, laptop (M2 Pro, Accelerate)**: all with the binary of this commit, Nov 9
  (`make test` with the settings-table check; `make shipcheck`; `make dnncheck`;
  bench standard fe14ac0 vs this tree: 19 instances, nodes exactly 1.000, time
  0.965; nug30 rule 9: 3,716 nodes, 3,574 DNN-fathomed; dumps old (environment
  and flags) vs new (par file) value for value identical through the map for
  the nug30 recipe, the tho30 CP recipe, the ship test and the full run; wil50
  depth-1 tasks r.0, r.10, r.14 (and r.20, which is r.14 again: the root has 15
  children) replayed from the root with fe14ac0 + `settings.sh`'s environment +
  `-x 1:4:1:300:1` against this binary + `fullrun.par`: identical trace lines
  (bound, status, 300 iterations each, children; pid and CPU aside) and
  identical children (49, 29, 29). **tho30 with the CP cut phase through depth
  2 gives 17,470 nodes (16,740 DNN-fathomed), not the 17,786 of the Nov 8
  log**, and the old binary fe14ac0 under the environment recipe gives the same
  17,470 on this machine, so the difference predates Phase 2 (commits between
  9d6e216 and fe14ac0: the Nov 8 default changes, the retirements, the
  tidying) and is not a settings fault; `experiments/tho30_solve/cutd2_cp.par`
  is the recipe as a par file. Worth a bisect before the next tho30 claim.
- **Gate, cloud box (after stage B)**: the 8 depth-2 replays with
  `sampling/validate.py` are still to run there, old binary (a667727 or
  fe14ac0) against this one, on the same machine.
  Run Nov 9, a667727 against c378cad: `experiments/wil50_feasibility/sampling/results/phase2_gate_wil50.md` (every DNN result identical; the QPB values differ, by 18ebd58 alone).


1. **Keyed settings in the par file.** Add an optional section after the
   positional blocks, so the 53 existing files parse unchanged:

   ```
   [settings]
   dnn.depths = 1:4
   dnn.maxgap = 1
   dnn.giveup = 300
   dnn.warm = 1
   dnn.ext = 20
   cut.grid = 5:10
   cut.side = 1
   cut.cp = 100:1000
   ...
   [block 1]
   branch.xbar_lines = 1
   ```

   - `[settings]` keys apply everywhere. `[block k]` keys override inside the
     k-th gap/depth block, for entries the registry marks as per-block.
   - An unknown key is an error, and so is a value out of range.
   - The syntax is a proposal; settle it in the implementation. The
     requirements are: backward compatible, keyed, validated, per-block where
     the policy varies over the tree.
2. **Migrate every surviving `policy` entry to a par key.** Use the families
   `dnn.*`, `cut.*`, `pool.*`, `ship.*`, `branch.*`, `est.*`.
   - The environment name is removed for every migrated entry. The par file is
     already the A/B mechanism for policy (`CONTRIBUTING.md`: a new par file
     passed to the candidate only).
   - The environment keeps only `experiment` entries that haven't graduated,
     and `diagnostic` ones (`QAP_DNN_TRACE`, `QAP_XBAR_DEBUG`). Diagnostics are
     printed in the dump but left out of the state-file hash.
3. **Policy flags become keys.**
   - `-x` maps to `dnn.depths`, `dnn.maxgap`, `dnn.giveup` and `dnn.warm`;
     `-w` maps to `fw.patience`. The flags stay as overrides (source `flag` in
     the dump) until the launchers below are converted, then go.
   - `-s` and `-d` are run-scoped and stay. The par file's positional sym and
     max-depth fields become deprecated: a warning if set, an error if they
     conflict with a flag. This fixes `util.C:126-130`.
4. **Compile-time policy becomes keys, with today's values as defaults:**
   - `branch.weighted` and `branch.weight` (`util.h:44-45`).
   - `branch.xbar_lines`, rule 9's number of candidate lines per side. Default
     1, which is today's two lines. Rule 9 keeps ignoring NBEST, so nothing
     changes. The branching variant (strong-branching a shortlist by DNN child
     bounds) comes later as a new rule on these keys.
5. **`dnntest` reads the same keys** (`dnntest -p file.par`). Its `-B`, `-L`, `-g`,
   `-G`, `-H` and `-d` become aliases printed as `flag`, so a cold node solve
   and the search can run from one par file.
6. **Changing a default is a search change.** It goes in its own commit, with a
   bench re-baseline, and every committed golden dump is updated in the same
   commit. This is the rule 34e347c would have needed.
7. **`fullrun/fullrun.par`** is built from `nug_deep_xbar9_all.par`, the
   `settings.sh` exports and `-x 1:4:1:300:1`. `settings.sh` then keeps only
   the instance, the sym file, the upper bound and the thread variables.
8. **The reference table is generated.** `bnbdrv --print-settings --all` prints
   the registry as a Markdown table. `RUNNING.md` includes it, and `make test`
   fails if the two differ, so the docs can't drift.

**Launchers to convert:**
- `bench/bench.py` and `bench/parallel.py` (`--extra`)
- `coord/` (qapctl and the worker supervisor) and `infra/phase1/` (recipes
  nug17 and nug30)
- `experiments/wil50_feasibility/sampling/common.py`, `fullrun/settings.sh`
- `pipeline_ab/run_pair.sh`, `tail_probe/`

Historical experiment scripts are not converted: their logs record the commit
they ran, and they reproduce on that commit.

**Gate (bit for bit):**
- **Dumps:** identical between the old configuration (environment and flags)
  and the new (par file) for:
  - the nug30 record recipe
  - the tho30 CP cut-phase run (`tho30_cutd2_cp_j8_m2pro.log`)
  - the full-run configuration
- **Trees:**
  - `make test`
  - `bench.py` standard, exactly 1.000
  - nug30 with rule 9, 3,716 nodes
  - tho30 with the CP cut phase through depth 2, 17,786 nodes
- **wil50, laptop:** depth-1 stage-A tasks `r.0`, `r.10`, `r.14` and `r.20`,
  replayed from the root with old and new binaries. They must give the same
  trace lines: the bound to 1e-9, the same iterations, the same line and
  children.
- **wil50, cloud box after stage B:** 8 depth-2 tasks replayed from their kept
  parent states, old binary against new, with `sampling/validate.py`:
  - 2 that fathomed in the first solve
  - 2 that fathomed in the CP phase
  - 2 that stalled
  - 2 that hit the iteration limit

## Re-gate the full run

After Phase 2, the full run's binary is new and its configuration is
`fullrun.par`.
- The golden dump must equal the settings stage A ran, and the wil50 replays
  must pass.
- If both hold, stage A's frontier (14,192 nodes) and stage B's class fractions
  remain valid inputs to the run and the estimate, because the settings are
  proven identical, not just believed to be.

## Out of scope, and related

- **Tuning.** Whether wil50 should adopt the record's beta 0.21 or the give-up
  confirmation zone of 0.003 is an A/B on the node-10 subtree, decided
  separately.
- **18ebd58 (review B1, the QPB h term)** changes QPB bounds. `update_STG` is
  called from both `qpbnd.C` and `qpbnd2.C`, and `qpbnd2` is what rules 4 and
  9 use to score children. So it can change QPB-scored line choices as well as
  QPB fathoms. It postdates stage A's binary. It is an algorithm change gated on
  its own; the settings work doesn't settle it. Before the run, decide whether
  the full run's binary carries it and, if so, A/B it on wil50 nodes.
- **The branching variant** lands after Phase 2, as a new rule reading
  `branch.*` keys.
