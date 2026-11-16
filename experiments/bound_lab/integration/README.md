# bound_lab integration: the day's solver patches in one tree (Nov 16 2026)

**IN PROGRESS: the verdict below is filled in as the runs finish. Sections marked (pending) are not done.**

## Verdict

(pending)

## What is in, what is out

Commits on `claude/cool-babbage-ilnogj-integration` (base 26894b2, whose `latest/` equals `cognac`'s
e1596e22: `git diff origin/cognac -- latest` is empty):

- **b86352e** (commit A): `latest/`, `RUNNING.md`'s settings table and five settings goldens. `all.patch` here is its
  `latest/` part against `cognac`.
- **e514979**: `experiments/dnn_cuts/dnncuts.C`, the driver's side of two of the patches (`-phase`, `-regions`,
  the cpfast flags). Kept apart because the brief asked for `latest/` alone in commit A; nothing in `latest/`
  depends on it, and the checks below pass with or without it.

| source | keys | in? |
|---|---|---|
| `tho40_search/stop/combined.patch` (handover + cpfast) | `cut.handover`, `cut.cp_wc_gap`, `cut.cp_wc_deficit`, `cut.cp_extrap` | in (applied clean) |
| `regions/regions.patch` | `cut.regions` | in (one conflict) |
| `admm/admm.patch` | `dnn.dscale`, `dnn.dscale0`, `dnn.dscale_until` | in (one conflict). Three keys, not the two the brief named: `dnn.dscale0` is in the patch |
| `origin/claude/cool-babbage-ilnogj-nug30` at 48a624c, `git diff origin/cognac -- latest` | `cut.giveup_depth` | in (applied clean) |
| `origin/claude/cool-babbage-ilnogj-tai50b` at 5337856, `tai50b_size/est_resume.patch` | none | **in**: a fix, see below |
| `origin/claude/cool-babbage-ilnogj-tho40dives` at 022973d, `tho40_dives/dive_trace.patch` | none (`DIVE_TRACE`, `DIVE_PATH`, `DIVE_FULL_DEPTH` read with `getenv`) | **out**: instrumentation, see below |
| `swapcuts.patch`, `cpfast_experiments.patch`, `cpfast_wc_test.patch`, `handover.patch`, `cpfast.patch` | | out, as the brief says (the last two are the parents of `combined.patch`) |

**`est_resume.patch`: in.** Ten lines in `driver.C` and `estimate.C`. Without it `bnbdrv -e -r F` estimates with
root bound 0, because the run that resumes never bounds the root: the dives' child weights and the parameter
level are then computed from a wrong relative gap. With it the root bound is F's. Measured here
(`logs/est_resume.log`, nug15, a state written at depth 2): unpatched prints `root bound: 0`, merged prints
`root bound: 992.015`, which is the state file's `root_bound` line. It touches nothing unless `-e` and `-r` are
both given; `make estimate` (cold `-e`) gives the same numbers in both builds. It has no key because it has no
"off" worth keeping. It is the `driver.C` and `estimate.C` parts of `all.patch` and nothing else touches those
two files, so it can be dropped or taken alone.

**`dive_trace.patch`: out.** It is that study's instrument: three environment variables read with `getenv`,
outside the settings table (CONTRIBUTING.md, "Where settings live": every variable goes in the registry),
a file opened and closed per node, and two modes that change what a dive does (`DIVE_FULL_DEPTH` pushes every
child, `DIVE_PATH` forces the child on a given permutation). It also sets `importance = 1.0` at every pop,
unconditionally. That line is needed only by its own all-children mode, but it is a change with no switch, in
the estimator's weighting: today a dive ends at its first node that does not branch, so the stale value is
never used, and I believe the line changes nothing, but I did not test that and it should not ride in on a
trace patch. If the dives study wants a permanent "dive along a permutation", it should come back as a
registered `QAP_EST_*` setting.

## Conflicts

Applied in the order combined, regions, admm, giveup_depth, est_resume, with `git apply --3way`.

1. **`latest/bnb.C`, regions on top of combined.** Both touch the line `if (cp.grid_rows == 0) cp.prefix = true;`:
   combined inserts `cut.handover`'s parse just above it, regions changes it to
   `if (cp.grid_rows == 0 && cp.regions.empty()) cp.prefix = true;`. Kept both: the handover block, then the
   regions form of the line (`bnb.C:978-991`).
2. **`latest/dnn.C`, admm on top of combined.** Both add code at the head of the ADMM loop: combined the
   `cp_extrap` state and its checkpoint/jump block, admm `dnn.dscale_until` and its switch back to the plain
   iteration. Kept both: the two sets of declarations before the loop, then inside it the dscale switch and
   `w.sc = sc_on`, then the extrapolation block (`dnn.C:1930-1975`). The order does not matter: extrapolation
   needs the cut hook (`ex_E` is 0 without `yc`, `dnn.C:1931`) and scaling is refused with it
   (`!par.ycuts`, `dnn.C:1748`), so the two are never on in the same call.

`settings.C`, `dnncut.C`, `dnncut.h`, `dnn.h`, `dnncuts.C`, `RUNNING.md`: merged without conflict. `RUNNING.md`'s
table row from `regions.patch` was then replaced by `tests/settings_doc.py --update` (all nine rows).

## Review findings

Read before the long runs, on the merged tree (line numbers are b86352e's). "Fixed": nothing. I found no plain
bug, so `latest/` is the authors' code plus the two conflict resolutions. What I would change is listed for
the owner instead.

### Defaults and initialisation

1. **Every new field has a default on every path I could find.** `DNNParams` (`dnn.C:251`), `DNNCutParams`
   (`dnncut.C:33-35`), the projector's hand-over members (`dnncut.C:715-716`), `Work::sc = false`
   (`dnn.C:537`; `Work` is a local of `dnn_bound`, `dnn.C:1670`, so nothing carries over between calls).
   `cut.handover`'s W is 0 (off) in the constructor and 200 when the key is given without it (`bnb.C:981`).

### State files

2. **The format does not change.** `state.C` is untouched, files are still `QAP-STATE 9`, and the pool encoding
   is the same: a pool carries its regions as index sets, so a pool made under `cut.regions` is self-contained.
3. **The settings hash changes, so an old state file is refused.** The table goes from 71 to 80 entries and the
   hash covers every entry. Measured (`logs/state_compat.log`, nug16a with `tests/par/ship.par`): a state
   written by the unpatched build is refused by the merged one (exit 1, "the state file was written with
   settings ..."), and with `--allow-settings-change` it finishes with the same totals as the unpatched build
   (32 / 30 / 0 nodes, 3,150 FW iterations). **For a run in flight this means all workers move to the new
   binary together, with `--allow-settings-change` on the first resume, or not at all.** The same is true of
   any commit that adds a key; it is the main thing to know before cherry-picking.
4. **`cut.regions` is a path, and only the path is in the record.** The hash covers the text of the value;
   `par.hash` covers the params file, not the regions file. A regions file edited between a stop and a resume,
   or two workers with different files at the same path, pass the `-r` check. The path must also exist on
   every worker (relative paths resolve from the solver's directory). The file is read at the first cut phase
   (`bnb.C:893`, a function-local static), so a wrong path stops the run with exit 2 at the first near miss
   rather than at startup (`logs/bad_values.log`), and never if no phase runs. Suggested before any
   distributed run with regions: read it at startup and add a `cut.regions.hash` entry next to `par.hash`.

### Paths where a key at its default still changes something

5. **Arithmetic and control flow: none found by reading, none in the identity runs** (below). The places that
   were rewritten rather than added to: `hh_block`, `rows_D`, `col_rot`, `rows_rot` (`dnn.C:663-765`) now pick
   the reflector by `w.sc`, the same vectors and the same operations when it is false; the warm-phase test
   `par.cut_round` became `par.cut_round == 1` (`dnn.C:1872`), and the only writers are `dnncut.C:1468` and
   `:2314` (1) and the two new keys (2); `stop()` returns early as before unless `ho_w > 0` (`dnncut.C:1270`).
6. **Output with `-v` does change.** A verbose non-CP solve prints a new line every 50 iterations
   (`split it N arrow share ...`, `dnn.C:2380`) whether or not `dnn.dscale` is set, and a verbose CP phase
   prints `cp phase start: ...` (`dnncut.C:2325`). `bnbdrv` never runs verbose, so traces and logs of a search
   are unaffected; `dnntest -v` and `dnncuts` output gain lines. Measured: `logs/dnntest_verbose.diff`
   (nug15, 300 iterations: every iteration line equal, the differences are the settings lines and the six
   `split it` lines). A script that counts lines of `dnntest -v` output would notice; `make dnncheck` does not.
7. **`dnn_bound` now looks up three settings by name on every call** (`dnn.C:1746-1748`, `:1938`), where the
   function's other lookups are `static`. A table scan per DNN solve: nothing measurable, a matter of style.

### Keys that interact

8. **Three rules can stop a CP phase, in this order** (`CutProjector::stop`, `dnncut.C:1268-1303`): the
   `cut.giveup` pairs first; then, if `cut.handover` is set, the hand-over rule **instead of** the stall rule
   (it returns before the stall test); else the stall rule. Consequences:
   - `cut.giveup_depth` works by not loading the pairs at shallower nodes (`bnb.C:968`), so above that depth
     only the hand-over rule (or the stall rule) can stop a phase. That is the combination nug30's `g2` and
     tho40's `zallhw` would make together; nobody has run it.
   - With both `cut.giveup` and `cut.handover`, a give-up pair can stop a phase at iteration 100 while the
     hand-over rule cannot act before max(W, 200) (`dnncut.C:1283`). The give-up wins early by construction.
   - With `cut.handover` set, `cut.stall` is ignored without a message.
   - A hand-over stop is reported as `stalled 1 gaveup 0` in the DNN trace, the same as a stall-rule stop.
     An analysis that counts stalls counts hand-overs.
   - The hand-over rule's "not before 200" and "not moved by 400: stop" (`dnncut.C:1283-1285`) are constants
     in the code, not fields of the key.
9. **`cut.handover`, `cut.giveup` and the three `cut.cp_*` keys do nothing in the rounds form** (no `cut.cp`):
   the `cut.cp_*` keys are read only inside the `cut.cp` branch (`bnb.C:1015-1018`) and the stop rules live in
   the CP projector. No warning is printed, unlike the rounds-only keys under `cut.cp` (`bnb.C:1023`).
10. **`dnn.dscale` is silently off in every CP-ADMM solve** and with the accelerations (`dnn.C:1748`). On an
    instance with a cut phase the first solves and the pooled children (a cut term in the objective, no cut
    hook) run scaled and the phases do not. Its author measured instances without a cut phase. The only run
    here with scaling, pools and phases together is the all-keys soundness run below, which says the bounds
    stay valid, not that the combination helps.
11. **`cut.cp_extrap` and `dnn.dscale` are never on in the same call** (finding 10 and `dnn.C:1931`), which is
    why the `dnn.C` conflict could keep both sides as they were.
12. `cut.cp_wc_gap` and `cut.cp_wc_deficit` together: either one keeps `w` (`dnncut.C:2318-2324`).

### Validation of values

13. **A malformed `cut.cp_extrap` is accepted** (`sscanf` with no check, `bnb.C:1017`). Measured
    (`logs/bad_values.log`): `abc` runs with extrapolation off, `40:x` runs with E = 40 and the default
    factor, guard and mode. `cut.handover`, `cut.giveup_depth` and `dnn.dscale` reject bad values with exit 2.
    `cut.grid` and `cut.cp` have always parsed this way, so this is the house style for the older text keys,
    but RUNNING.md promises that a value that does not parse stops the run.
14. **`cut.regions` indices out of range are dropped without a message in `bnbdrv`** (`dnncut.C:612`), while
    `dnncuts -regions` rejects them. A 1-based file would run. Any region gives valid cuts, so this is about
    getting the regions that were meant, not about soundness. A line longer than 16 KB would be split into two
    regions (`bnb.C:899`); no file here comes near that.
15. **Six of the new entries carry an "old environment name" that never existed** (`DNN_DSCALE`, `DNN_DSCALE0`,
    `DNN_DSCALE_UNTIL`, `QAP_DNN_CUT_REGIONS`, `QAP_DNN_CUT_GIVEUP_DEPTH`, `QAP_DNN_CUT_HANDOVER`:
    `settings.C:47-49, 61, 81, 83`), so `bnbdrv --settings-map` lists them as retired names; the
    cpfast keys pass `NULL`. Harmless: setting one in the environment is an error that names the key.
16. `dnn.dscale0 = 0` is in range and means 1 (`dnn.C:1768`).

### Seen in passing, not from the merge

17. A phase that finds no regions at a node returns before doing anything and is traced as
    `status error ... iters 0` (`dnncut.C:2263`; once in the all-keys run on tai15b, depth 1). The node goes on
    with its first bound and is counted as a cut-phase node. The path is old; region files make it more
    reachable than grids do. I did not establish whether the cause there was an empty region list or no
    violated triangle.
18. Five goldens under `experiments/estimates/sko42/tune/` and `sko49/tune/` fail `fullrun_check` with the
    **unpatched** binary too (they say `settings 65`; the table has had 71 since before today). I left them
    alone (`logs/goldens.log`).

### What I did not review

- `experiments/dnn_cuts/dnncuts.C` beyond a skim of the diff.
- The algebra of the scaled transform (`Work::set_scaling` and the scaled Y/Z step, `dnn.C:546-592`,
  `:2059-2108`). What I relied on: every certificate is computed with `w.sc = false` on a multiplier in
  original units (`dnn.C:2224`, `:2496`), so a mistake there could cost iterations but not validity; the
  author's comparison against a dense prototype; and the all-keys run.
- Whether the hand-over rule's cost model or cpfast's thresholds are right. That is the authors' studies.

(checks, identity and reproductions: pending)
