# bound_lab integration: the day's solver patches in one tree (Nov 16 2026)

One solver tree with all of today's solver patches, checked against the tree before them. Commit to take:
**b86352e** (`latest/`, RUNNING.md's settings table, five settings goldens); `all.patch` is its `latest/` part
against `cognac` e1596e22. [m] = measured here, with the log named.

## Verdict

- **Builds clean**: yes, from a clean checkout (exit 0; the same three old compiler warnings as before the merge).
- **The repo's checks pass**: `make test`, `make nodetrace`, `make shipcheck`, `make dnncheck` pass.
  **`make continuationcheck` fails, and fails identically before the merge**
  (`FAIL split Y/Z differs from uninterrupted full-spectrum state`); I did not find out why.
  RUNNING.md's table and five goldens were regenerated for the nine new keys and nothing else.
- **Switch-off identity is exact**: with no new key set, the unpatched and merged builds give the same node
  trace and the same DNN trace, line for line, on six runs: nug20 in full (1,138,754 trace lines), nug20 with the
  DNN bound, a 404-node tho40 slice under production settings, a tho40 slice with the cut phase, a tai30b
  slice, and a 232-node wil50 slice under production settings. No difference of any kind was found.
- **Each patch still does what its author measured**: the hand-over subtree `38:32` takes 410 iterations (410);
  tai30b with region cuts is solved in 363 nodes (363), the same trace; `dnn.dscale` with `dscale_until = 100`
  takes 13,965 iterations on the tai30b subtree (13,965), every node equal; nug30 under `cut.giveup_depth = 2`
  takes 77,534 iterations (77,535) over the same 1,646 nodes (43 of them differ from the author's in the
  fourth decimal of a bound, as two parallel runs of one build do; see the reproductions).
- **The one thing to know before taking it: the settings hash changes** (71 to 80 entries), so the merged
  binary refuses a state file written before it unless `--allow-settings-change` is given [m]. With the flag
  the resumed run's totals are unchanged. A run in flight has to switch binaries on all workers at once.
- No bug fixed, because I found none that was plainly one. 18 review findings are listed below; the ones I
  would act on before a long run are 4 (`cut.regions` is a path whose contents are not in the settings
  record), 8 (how the three stop rules combine), 13 and 14 (two keys accept bad values silently).

About 1 h 50 min of wall time on a 4-core container; nothing outside it was touched.

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

## The repo's checks

Built from nothing in a clean checkout of e514979 (`git worktree add --detach /home/user/wt-clean`),
`clean_checks.sh`, logs `logs/clean_*.log`. The machine: 4 cores, Xeon with AMX, g++ 13.3.0, OpenBLAS from the
`scipy-openblas32` wheel installed today (the container was new; the authors' builds used the same wheel name,
I do not know that it was the same version). The `bnbdrv` from the clean checkout has the same sha256
(cec0f39a8ae8...) as the one the identity runs and reproductions used; the unpatched one is e6a5558d99c8...
(worktree `/home/user/wt-base` at 26894b2), and rebuilds to the same bytes.

| check | command | outcome [m] |
|---|---|---|
| build | `make -C latest -j2`, `make -C latest dnntest`, `make -C latest dnncuts` | exit 0. 3 compiler warnings (`console.C:283`, `console.C:441`, `driver.C:543`), the same 3 as the unpatched tree, none in a file the merge touches except `driver.C`, where the line is old (`system("hostname")`) |
| `make test` | brute force, resume, parallel, coordinator, settings table, regression, estimator | **exit 0.** "all parallel checks passed", "all coordinator scenarios passed" (the AWS scenario skipped: no boto3 here), "settings table: RUNNING.md agrees", "regression check passed: 7 instances match bench/reference/ (tolerance 0.0%)", "estimator: all checks passed" |
| `make nodetrace` | | exit 0, "node trace: all checks passed" |
| `make shipcheck` | | exit 0: ship64 and ship64store "every node identical to the straight run" (184 DNN nodes, 13,231 iterations) |
| `make dnncheck` | | exit 0, "OK" (prototype agreement 8.6e-08, nug20 nodes with known optima, cancellation, eliminations) |
| `make continuationcheck` | | **exit 2: `FAIL split Y/Z differs from uninterrupted full-spectrum state`. The unpatched tree fails with the same line** (`logs/base_continuationcheck.log`), so it is not from the merge. Not chased: I do not know whether it fails on the owner's machine or only with this container's BLAS. `make test` does not run this target |
| settings goldens | `fullrun_check` per settings file, `logs/goldens.log` | the five that passed with the unpatched binary fail with the merged one, as expected, and pass after regeneration; diff in `logs/goldens.diff` |

The unpatched tree's `make test` (`logs/make_test_base.log`), run at the same time as the merged tree's first one
and while the identity runs started, failed once in the estimator test, on CPU time alone
(`FAIL QAP_EST_CHILD_QPB=0: 1200 dives, nodes 1.031, CPU 1.167 of the run's`: the tolerance is 15% and the box
was loaded). Three reruns pass (`logs/base_estimate_rerun.log`). The node ratio, 1.031, is the same in both
trees every time.

The coordinator test's `noub` scenario (no incumbent given) reports different node counts from run to run in
both trees (101,528 to 113,419); it is marked ok by the test and depends on when workers find the incumbent.

### The goldens

`RUNNING.md`'s table: nine rows added (`tests/settings_doc.py --update`). Each of the five goldens changes in
the same way, shown here for `experiments/estimates/tho40/settings.golden` (all in `logs/goldens.diff`):

```
-settings 71 hash d8f32b9fa3087c58
+settings 80 hash 49171709a99d1b81
+setting dnn.dscale 0 default policy
+setting dnn.dscale0 1 default policy
+setting dnn.dscale_until 0 default policy
+setting cut.regions "" default policy
+setting cut.giveup_depth 0 default policy
+setting cut.handover "" default policy
+setting cut.cp_extrap "" default policy
+setting cut.cp_wc_gap 0 default policy
+setting cut.cp_wc_deficit 0 default policy
```

Regenerated: `experiments/estimates/{tho40,sko42,sko49,tai50b}/settings.golden` and
`experiments/wil50_feasibility/fullrun/settings.golden`. Not touched: the five stale ones under `tune/` (finding 18).

## Switch-off identity

Same command, same machine, the unpatched build (26894b2, `latest/` = `cognac` e1596e22) against the merged
one (b86352e), no new key set. `ident.sh NAME base|new ARGS`, `ident_queue.sh`; `ident_cmp.py` compares the
node traces (`QAP_NODE_TRACE`) and the DNN traces (`QAP_DNN_TRACE`) line by line in order, with only the CPU
fields and the pid removed, so every bound to four decimals, every iteration count, status, pool size,
branching line and child count must agree. All [m]; output in `logs/ident_cmp.log`.

| run | settings | node-trace lines | DNN-trace lines | result |
|---|---|---|---|---|
| `nug20`: full search, no incumbent given | `latest/PARAM/nug20.par` (QPB only) | 1,138,754 | none | **identical** (1,138,620 / 911,370 / 103 nodes in both) |
| `nug20_dnn`: full search, `-x 2:6:1:300:1 -u 2570` | the same with the DNN bound | 140 | 133 | **identical** |
| `tho40_a`: `-F 26:10,8:0,2:32,18:39,35:24` | `experiments/estimates/tho40/tho40.par` (production; no cut phase) | 404 | 404 | **identical** |
| `tho40_t3_z030`: `-F 26:10,8:0,38:32` | `tho40_search/par/z030.par` (grid cut phase in CP form, pools: keys from before the merge) | 38 | 38 | **identical** |
| `tai30b_F3`: `-F 0:3,1:7,2:10` | `admm/bb_tai30b_base.par` (tai50b's settings; no grid) | 301 | 233 | **identical** |
| `wil50_a220`: `-F 21:22,31:23,40:15` | `experiments/wil50_feasibility/fullrun/fullrun.par` (production; cut phase, pools), `QAP_DNN_POOL_CHECK=20` | 232 | 232 | **identical** (232 / 227 / 0 nodes, 7 cut phases, 17,191 DNN iterations in both) |

- The tho40 production slice is `rule2_ab`'s subtree a; the cut-phase slice is `tho40_search`'s `t3`. I added the
  second because tho40's production settings have no cut phase and most of the merge is in the cut phase.
  The nug20-with-DNN and tai30b pairs are extras (cheap).
- **The three slices that have an author's trace also equal it**, line for line (`logs/ident_vs_authors.log`:
  `rule2_ab/runs/a_base`, `tho40_search/runs/t3_z030`, `admm/logs/bb_tai30b_F3_base`), so this container's
  build reproduces the authors' builds on those commands.
- CPU of each pair, indicative (shared box): within 1.2% of each other (`logs/ident_cpu.log`), e.g. tho40_a
  694 s and 693 s. The `w.sc` tests in the transform cost nothing visible.
- Node traces of the two nug20 full searches are 96 MB each and are not committed; their line count and digest
  are in `logs/ident_cmp.log`, and `ident_queue.sh` regenerates them in 2.5 minutes.
- The wil50 slice is the depth-3 subtree `a220` of `experiments/estimates/wil50/subtrees/`. The owner's run of it
  (another machine and BLAS, through the coordinator) has the same 232 DNN nodes, 227 fathomed and 7 phases, and
  18,852 iterations where this container takes 17,191; no line of the two DNN traces agrees to the fourth
  decimal, which is the known dependence on the BLAS, not something the merge could cause (both builds here
  give 17,191).

## Reproductions, one per patch

Merged build, `repro.sh`, one process at a time at low priority beside the identity runs (the nug30 whole tree
separately, four workers). Iterations and nodes are the metric. All [m]; comparisons in `logs/repro_cmp.log` and
`logs/repro_handover.log` (`repro_cmp.py` matches the two DNN traces node by node).

| patch | run | author | here | |
|---|---|---|---|---|
| `cut.handover` + `wc` (combined) | parent `26:10,8:0` bounded under `zallhw`, then the warm subtree `38:32` from its state | parent: phase stopped by the rule at 2,550 iterations, 237,004.5614. `38:32`: 110 + 300 = **410 iterations**, fathomed at 240,515.9671, 1 node | the same numbers, to the fourth decimal and the pool sizes (`kept 1428 dropped 145`, `active 1735`) | **exact** |
| `cut.regions` | tai30b, `regions/bnb/tai30b_cut.par`, `-u 637117114`, in one piece | 637,117,113 proved in **363 nodes** (DNN 299 nodes, 250 fathomed; cut phase 49 nodes, 22 fathomed), in three resumed pieces | 363 nodes, 299 / 250, 49 / 22; the node trace's 375 lines equal the author's three pieces concatenated | **exact** |
| `dnn.dscale = 0.5`, `dnn.dscale_until = 100` | tai30b `-F 0:3,1:7,2:10`, `admm/bb_tai30b_ds0.5u100.par` | **195 DNN nodes, 184 fathomed, 13,965 iterations** (243 nodes) | 195, 184, 13,965 (243); all 195 nodes equal in status, iterations and bound | **exact** |
| `cut.giveup_depth = 2` | the author's command: the whole nug30 tree under `nug30_phase`'s `g2.par`, `bench/parallel.py`, 4 workers, `--ub 6124` | 1,646 nodes; 1,645 DNN nodes, 1,583 fathomed; 85 phases, 53 stopped by the give-up; **77,535 iterations** (61,900 + 15,635). The nine depth-1 phases run to the stall rule: 630-1,240 iterations, 6,570 in all | 1,646; 1,645, 1,583; 85, 53; **77,534** (61,899 + 15,635). Depths 1 and 2 (234 nodes): every node equal in every field. Below: 43 of 1,411 nodes differ from the fourth decimal (first: depth 3, 6093.1458 here, 6093.1454 there) | **same tree, 1 iteration in 77,535** |
| | the same settings in one process, `-u 6124 -d 2` (run first, as the cheap check) | (the author's depth-1 and depth-2 lines) | depth 1 exact; depth 2: 20 of 225 nodes differ from the fourth decimal, 14 more iterations in 30,258 | see below |
| `est_resume` | `bnbdrv -e -r` on a nug15 state | (no number in the branch) | root bound 992.015 from the state file, 0 without the patch (`logs/est_resume.log`) | works as described |

- **nug30: I put the differences down to how a parent's state reaches a child, not to the merge.** With `ship.bits = 32` a
  child bounded by another process starts from its parent's state rounded to 32 bits, and one bounded by the
  same process starts from it at full precision. My one-process `-d 2` run therefore differed from the
  author's four-worker run at depth 2 (20 nodes, fourth decimal). Run the author's way, depths 1 and 2 agree in
  every field, and the 43 nodes that still differ are all at depth 3 or below, where (as I read
  `bench/parallel.py`; not tested) which process bounds a node depends on timing. For scale, the author's own two runs of the unpatched baseline arm (`runs/base`,
  `runs/base2`, same command) differ from each other in 8 of 3,687 nodes in the same way
  (`logs/repro_cmp.log`). 43 of 1,645 is a larger share than 8 of 3,687. I did not go further
  than that. Node count, fathom counts, phase count and give-up count are equal.
- The author's `cut.giveup_depth` build was the nug30 branch's own; the other three authors built their patch
  alone. So these four rows also say that each patch behaves the same with the others compiled in.

## Checks beyond the brief

- **Every new key on at once, small instances with known optima** (`allkeys/`, `logs/allkeys.log`) [m]: tai12b and
  tai15b with region files, nug15 and nug16b with a grid (nug16b with a region file as well), pools,
  `dnn.dscale = 0.5` with `dnn.dscale_until`, `cut.handover`, `cut.giveup` with `cut.giveup_depth = 2`,
  `cut.cp_extrap`, both `cut.cp_wc_*` keys, `wc=2.88`, `QAP_DNN_POOL_CHECK=20`, incumbent = optimum + 1 (+ 2 on
  nug). All four find the optimum and prove it (36, 60, 31, 34 nodes; 8 to 12 cut phases each, most stopped at
  iteration 200 by the hand-over rule, two by the give-up). This is a soundness check of the combination, on
  instances where a wrong fathom on the optimal path would lose the optimum. It says nothing about speed.
- **Stop and resume with every new key on** (`logs/shiptest_newkeys.log`) [m]: `tests/ship_test.py` with the keys
  passed by `--set` (nug20, CP phase, hand-over, give-up from depth 2, extrapolation, a regions file on top of
  the grid, scaling until 100): "ship64: 12 slices, 184 DNN nodes, 183 warm-started, 10,155 DNN iterations:
  every node identical to the straight run, and the 3 level tables", and the same through the parent-state
  store. So a search under the new keys resumed from state files is the same search.
- **A state file across the two builds** (`logs/state_compat.log`): finding 3.
- **Malformed values** (`logs/bad_values.log`): findings 4, 13.
- **`dnntest -v` output** (`logs/dnntest_verbose.diff`): finding 6.

## What I did not verify

- sko49 (the brief allowed wil50 or sko49; I ran wil50).
- `make continuationcheck` fails before and after the merge; I did not find out why.
- The coordinator's AWS scenario (no boto3/moto in this container) and anything on EC2.
- The four reproductions are one cheap run each. I did not rerun any author's larger tables (the `18:32`,
  `30:32`, `2:32` subtrees, the path checks under the combined arm, tai50b dives with scaling, the whole nug30
  tree under `g2` or `g5`).
- Combinations nobody has measured for speed: `cut.giveup_depth` with `cut.handover`; `dnn.dscale` on an
  instance with a cut phase; `cut.regions` together with `cut.grid` (run only in the all-keys and ship tests).
- Why the nug30 whole-tree run differs from the author's in 43 nodes' fourth decimals beyond "parallel runs
  under 32-bit shipping are not reproducible to the last digit" (above).
- `dive_trace.patch`'s unconditional `importance = 1.0` (left out; see "What is in").
- `dnncuts` was built (`make -C latest dnncuts`, exit 0) but not run.
- The other branches were fetched once, at 48a624c (nug30), 022973d (tho40dives) and 5337856 (tai50b); neither
  of the last two had a README yet, so my reading of their patches is from the code. Later commits there are
  not in this merge.

## Files

- `all.patch`: `git diff origin/cognac b86352e -- latest` (863 lines, eight files). It applies to `cognac`'s
  `latest/` with `git apply`. Cherry-picking b86352e also brings RUNNING.md's table and the five goldens,
  which `make test` and `fullrun_check` need.
- `ident.sh`, `ident_queue.sh`, `ident_cmp.py`, `ident/`: the identity runs. `repro.sh`, `repro_cmp.py`, `repro/`:
  the reproductions (state files, 44 MB each, are not committed). `clean_checks.sh`: the clean-checkout checks.
  `allkeys/`: the all-keys runs. `logs/`: everything quoted above.
- The scripts name this container's paths (`/home/user/qap`, `/home/user/wt-base`, `/home/user/wt-clean`).
