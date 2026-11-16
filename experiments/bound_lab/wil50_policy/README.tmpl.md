# wil50 / sko49: tho40's phase policy on complete subtrees (bound_lab `wil50_policy`, Nov 16 2026)

Question (the coordinator's brief): do the three changes to how the cut phase is *used* that halved
tho40's warm subtrees (`../tho40_search/README.md`, "The stop rule": a wide zone, the hand-over rule
`cut.handover`, early entry) make wil50's and sko49's subtrees cheaper? **[m]** measured here (the run
under `runs/` is named), **[e]** estimated. Iterations are the metric; CPU is this box's (4 cores, four
solver processes at a time, `lscpu.txt`).

## Answer

**A modest gain, not tho40's 2x: about 1.1-1.35x on wil50 and 1.2-1.35x on sko49, from two of the three
changes. The wide zone changes nothing, and the hand-over rule has to be refitted to be *less* patient
far from the cutoff than tho40's, or it loses most of the gain.**

| | subtree (cold `-F`) | base | recommended arm | factor in iterations (CPU) | the tree |
|---|---|---|---|---|---|
| wil50 | a220 `21:22,31:23,40:15` | 232 nodes, 17,191 it | `hoe2`: 12,699 it | **1.35x** (1.36x) [m] | the same 232 nodes, node for node |
| wil50 | b83 `8:9,21:11,45:40,34:0` (a depth-4 child of a251) | 181 nodes, 13,296 it | `hoe2`: 12,230 it | **1.09x** (1.07x) [m] | 181 nodes in both |
| wil50 | a251, a325 (laptop traces, not rerun: 35 and 100 CPU-hours an arm here) | 387,735 and 1,066,452 it | | **1.17x and 1.16x** [e, `extrap.py`] | |
| sko49 | s188 `20:6,40:4,41:16` | 47 nodes, 4,162 it | `hoe2w`: 3,446 it | **1.21x** (1.23x) [m] | the same 47 nodes |
| sko49 | s232 `23:0,40:8,41:10` | 271 nodes, 34,557 it | `hoe2w`: 25,763 it | **1.34x** (1.37x) [m]; 1.56x on the 271 nodes both arms bounded | 404 nodes: 133 more, cheap ones (3,627 it) |

- **wil50: about 1.15-1.2x [e], measured 1.09x and 1.35x.** Both measured factors are inside the 1.4x
  that single subtrees vary by, but they are not noise of that kind: the trees are the same and the
  whole difference sits on the few nodes that run a phase (a220: 7 nodes, 10,649 -> 6,111 it, 1.74x;
  the 225 nodes that fathom in their first solve: 6,542 against 6,588). What limits it is the
  composition: on the large subtrees 46-48% of the iterations are first solves that fathom, which no
  phase policy touches. a220 overstates the gain (44% of its iterations are missed phases, one of
  them a 3,400-iteration phase that ended 0.1 short); a251 and a325 have 25%.
- **sko49: about 1.2-1.35x [m, two subtrees]**, with `wc=1.4` in the arm. Attribution: on s188 the
  arm without `wc` (`hoe2`) is 1.04x and `wc=1.4` alone 1.05x, the two together 1.21x; on s232
  `wc=1.4` alone is 1.08x (1.17x on matched nodes) and `hoe2` alone was stopped at the session's budget after 73 of its nodes: on the 73 nodes it
  shares with the others, base took 10,623 iterations, `wc14` 9,166, `hoe2` 8,769 (1.21x) and
  `hoe2w` 7,726 (1.37x) [m, partial: `runs/s232_hoe2`, incomplete]. So on sko49 each half gives
  about 1.05-1.2x and they add.
- **Recommended blocks** (on top of each instance's current par file; `cut.handover` needs
  `../tho40_search/stop/handover.patch`, which this branch has applied in `latest/`):

  ```
  # wil50 (arm hoe2)                      # sko49 (arm hoe2w)
  dnn.giveup_warm = 100                   dnn.giveup_warm = 100
  cut.handover = 32:0:25                  cut.handover = 27:0:5
                                          dnn.beta_sched = two,f=1.6667,t=60,w=0.25,wc=1.4
  ```
  Not `wc=1.4` on wil50: it made both wil50 subtrees worse than `hoe2` (below).
- **What did not carry over from tho40.**
  1. *The wide zone is already there.* Under the base zones (0.0075: 366 short; 0.018: 421 short)
     every unfathomed node of every subtree here, and of the laptop's a220, a251 and a325, already
     runs a phase (no "no phase" node in any trace). The zone arm is the base arm, bit for bit, at
     depth 3 and below; it could only differ at depths 1-2, outside these subtrees. Not run.
  2. *Patience far from the cutoff is not repaid.* With constants in tho40's spirit (`32:0:250`, arm
     `hoe`) a220's root phase ran 2,660 iterations instead of 980 for 14.7 units more (105.3 -> 90.6
     short), and its children's first solves cost the same (5,405 against 5,314 on the 179 matched
     ones): 1,680 iterations lost, 1.19x instead of 1.35x. The phase's tail is flat here (1 unit per
     100 iterations after 1,500), where tho40's far phases were still closing 20 units per 100.
  3. *The gain is near the cutoff*: phases that the stall rule lets run 1,600-3,400 iterations to end
     0.1-6 short, whose 46 children then fathom in 25-30 iterations each (1,150-1,250 in all). The
     rule stops them at 1,000-1,150.
  4. *Early entry cuts both ways.* It saves 190 iterations at each give-up (first solve 110
     instead of 300), but a warm node that would have fathomed in its first solve between iterations
     110 and 300 now gives up and runs a phase: on b83, 7 such nodes cost 2,701 iterations where
     they cost 2,000 (and the phases fathom them). That is why b83 is 1.09x and not more.
- **`wc` (the coordinator's added arm).** On the far roots `wc=1.4` reaches the base phase's end value
  1.27x sooner on wil50 and 1.44x sooner on sko49; 2.8 is no better and 5.6 worse (table below). In
  the search, alone: sko49 s188 1.05x, s232 1.08x with a larger tree (359 nodes: two phases that
  fathomed after 2,050-2,460 iterations in base stalled 0.7 and 3.5 short); 1.17x on the matched
  nodes. On wil50 on top of `hoe2` it loses: a220 13,451 against 12,699 (277 nodes against 232), b83
  13,645 against 12,230 (268 nodes against 181; 0.97x of base). Near
  misses fathom less often at the higher penalty, as cpfast found on tho40; `cut.cp_wc_gap` (cpfast's
  patch, not in this build) is the thing to try before using `wc` on wil50.
- **Sanity: PASS.** No run printed `New objective value` (wil50 48,816, sko49 23,386 stand) and every
  finished run ends with its subtree fathomed (columns "new obj" and "complete" below). The unpatched
  keys reproduce the documented trees: a220 232 DNN nodes (1, 47, 184), s188 47 nodes and 4,162
  iterations to the digit, s232 271 nodes (1, 46, 180, 44).

## How it was run

- Solver: `latest/bnbdrv` of this branch = `claude/cool-babbage-ilnogj` + `handover.patch` (committed
  here; unset, the key changes nothing; the settings goldens differ by the one new unset key and are
  not regenerated). `pip install scipy-openblas32 scipy`, `make -C latest -j4`.
- `run.sh INST NAME ARM SPEC`: one `bnbdrv -T 1`-style process (one thread) per subtree and arm, cold
  from `-F`, with the instance's symmetry file and incumbent, `QAP_DNN_TRACE` + `QAP_DNN_TRACE_HIST` +
  `QAP_NODE_TRACE` and the full run's `QAP_DNN_POOL_CHECK=20`. `queue.sh` keeps four running from
  `jobs.txt` (`queue.log` has the order). Each `runs/NAME_ARM/` has `cmd.txt`, `log.txt`, `err.txt`,
  `dnn.txt` (one line per DNN node, with each phase's history), `nodes.txt`.
- Arms (`par/INST_ARM.par` = the base file plus the lines at its end):

  | arm | wil50 | sko49 |
  |---|---|---|
  | `base` | `experiments/wil50_feasibility/fullrun/fullrun.par` | `experiments/estimates/sko49/sko49.par` (Cw) |
  | zone | = base (see above) | = base |
  | early | `dnn.giveup_warm = 100` (both instances already have `dnn.gu_zone = 0`): inside `hoe*`, not run alone | the same |
  | `hoe` | early + `cut.handover = 32:0:250` (the first fit, patient far away) | early + `27:0:16` |
  | `hoe2` | early + `cut.handover = 32:0:25` (refit after a220 `hoe`) | early + `27:0:5` |
  | `hoe2w` | `hoe2` + `wc=1.4` | `hoe2` + `wc=1.4` |
  | `wc14`, `wc28`, `wc56` | base with `wc` = 1.4, 2.8, 5.6 in `dnn.beta_sched` | the same |

- **The stall rule "on and off" with the hand-over rule:** in the patch the hand-over rule *replaces*
  the stall rule (`CutProjector::stop()`), so every `ho*` arm is "stall off"; the primal gate
  (`cut.pgate`) is a rounds-form key that the CP form ignores, and neither instance sets `cut.giveup`.
  The base arm is "stall on, no hand-over". A both-rules arm would need a code change; not done.
- Warm roots (the parent's state, as the tho40 study): not built. `experiments/estimates/README.md`
  ("Cold against warm") measured a220 warm and pooled at the same 232 nodes and 4% fewer iterations.
- Subtrees. a220 and sko49's two are the brief's. a251 (387,735 iterations on the laptop: about 35
  CPU-hours an arm on this box, against 17 core-hours in all) could not be run; **b83**, a depth-4
  child of a251 (`21:11,34:0,45:40,8:9`: 82.6 short and 180 nodes below it in the laptop's a251), was
  run cold instead as the second wil50 subtree.

## Every run [m]

`python3 summarize.py`; "stalled" is "handed over" in the `ho*` arms; `a220r`, `s232r` are the
single-phase `wc` probes (`-d 3`). x base is base / arm. `s232_hoe2` reads RUNNING: it was stopped by
hand at 07:00 (its `cmd.txt`), incomplete, and is used only for the 73-node comparison in the answer.

@@ python3 summarize.py

Matched nodes (`python3 summarize.py --matched NAME ARM1 ARM2`: the nodes both arms bounded, by how
they ended; where the trees have the same node count but fewer "both", the branching facility
differed under one node, whose 45-46 children then carry other names):

```
@@ for p in "a220 base hoe" "a220 base hoe2" "a220 hoe2 hoe2w" "b83 base hoe2" "b83 hoe2 hoe2w" "s188 base hoe2w" "s232 base wc14" "s232 base hoe2w" "s232 base hoe2" "s232 hoe2 hoe2w" "s188 base hoe2"; do set -- $p; [ -f runs/$1_$3/done ] && python3 summarize.py --matched $p; done
```

## What the phases did [m]

Every node that ran a phase (`python3 summarize.py --phases RUN`): first solve (iterations, shortfall),
the phase (iterations, how it ended, gap left), and the DNN nodes and iterations below it.

a220 base:
```
@@ python3 summarize.py --phases a220_base
```
a220 `hoe` (patient far away: the root's 2,660 iterations):
```
@@ python3 summarize.py --phases a220_hoe
```
a220 `hoe2` (recommended):
```
@@ python3 summarize.py --phases a220_hoe2
```
b83 base and `hoe2`:
```
@@ python3 summarize.py --phases b83_base; python3 summarize.py --phases b83_hoe2
```
s232 base and `hoe2w`:
```
@@ python3 summarize.py --phases s232_base; python3 summarize.py --phases s232_hoe2w
```

- Under `hoe2` no phase runs past 1,140 iterations on wil50; base's ran 1,540, 1,620 and 3,400 to end
  4.5, 4.5 and 0.1 short. The children of the three cost 1,125-1,357 either way.
- On s232 the rule also stops three phases that in base *fathomed* after 1,800, 2,050 and 2,460
  iterations (at 2.0, 5.0 and 7.9 short, after 1,140, 770 and 690): three lost fathoms, each replaced
  by 44-45 children that cost 1,190-1,230 iterations. That is the 133 extra nodes, and it is a gain.
- Phases that fathom quickly are unchanged or shorter (they start from a first solve stopped at 110).

## The `wc` sweep on single phases [m]

`python3 wcsweep.py`: the root of a220 and of s232, cold, first solve 300 iterations, then the phase
under the stall rule; gap to the cutoff by phase iteration (the trace's history), and the iteration at
which each `wc` first reaches the gaps the base phase had at 300, 600 and at its end. Far misses only:
a pooled near miss needs its parent's state, so the near side was measured in the search arms instead
(`wc14`, `hoe2w` above).

@@ python3 wcsweep.py

## The hand-over cost, refitted [m data, e model]

`python3 fitcost.py 50 48816 ../../estimates/wil50/subtrees/{a251,a325,a220}/trace.txt`: the 364
missed phases of the laptop's complete wil50 subtrees, by the gap they left, and what their children
then cost per child (m = 50 - depth children):

```
@@ python3 fitcost.py 50 48816 ../../estimates/wil50/subtrees/a251/trace.txt ../../estimates/wil50/subtrees/a325/trace.txt ../../estimates/wil50/subtrees/a220/trace.txt | sed -n 1,13p
```

- 32 iterations a child at the cutoff, doubling about every 22 units from 25 short on: far steeper
  than a quadratic. In the key's unit (g = 0.1% of the cutoff, 48.8 units) the first fit was
  `32:0:250`. sko49 (`fitcost.py 49 23386` on Cs232 and this study's base runs: 13 missed phases)
  is 27-33 a child up to 20 short, 45 at 31, 130 at 50-60, 725 at 121: `27:0:16`.
- **That curve is not what a longer phase buys.** It compares *different nodes* that the same stall
  rule stopped at different gaps: a node left 100 short is a harder node than one left 85 short. Run
  as the rule (`hoe`), a220's root bought 14.7 units for 1,680 iterations and its children were no
  cheaper. `hoe2`'s constants (wil50 `32:0:25`, sko49 `27:0:5`) keep the measured cost at the cutoff
  (N0) and cut the far term so that a far phase stops about where the stall rule stops it (a220's
  root at 1,120, s232's at 1,040-1,160, b83's at 1,010). They are set from that one observation,
  not fitted: the far term is the open constant.
- Offline replay (`replay.py`, after tho40's `replay_rules.py`) on the base runs' histories, written
  before the search arms: `replay_a220_base.out` (7 phases) and `replay_sko49_base.out` (22). Its
  children column uses the curve above and so shares its flaw; its phase-iteration column said what
  the runs then showed on wil50 (8,499 -> 4,829 phase iterations under `32:0:25`; measured 5,151). On
  sko49 it predicted no gain from any stop rule (x1.00-1.02) because it charges a lost fathom at the
  curve's price; the run found the lost fathoms cheap.

## wil50's larger subtrees [e]

`python3 extrap.py`: the laptop's complete subtrees by where their iterations go, and what the two
effects measured on a220 would save there (early: each warm give-up's first solve ends at 110; near
hand-over: a missed pooled phase ends by 1,050 iterations, its children's cost unchanged). It does not
count early entry's cost on slow first-solve fathoms (b83), so read it as an upper estimate.

@@ python3 extrap.py

By this model a220 itself is 1.30x (measured 1.35x on this box's own base). For the tree: about
**1.15x, 1.1-1.2x [e]**, i.e. wil50's 7,500-18,600 core-hours less 10-15%.

## Not done, and what I would do next

- a251 complete (too large for this session); sko49 s232 `hoe2` to its end; warm roots; an early-only arm and a hand-over-only arm
  on a220 (the matched nodes above separate the two effects, an arm of each would do it cleanly).
- `dnn.giveup_warm` between 100 and 300 (150, 200): b83 says 100 is too early for the nodes that
  fathom slowly in their first solve.
- The far term of `cut.handover` on a deep subtree. a220 and b83 are shallow; in a251 79% of the
  iterations sit under two depth-4 nodes whose pooled phases stalled at 600 iterations, 108 and 143
  short. Whether a longer (or `wc=1.4`) phase there pays is the largest open question for wil50,
  and needs a251's `0:0` or `47:0` child run both ways (2-20 CPU-hours an arm).
- `wc=1.4` with `cut.cp_wc_gap` (near misses at the old penalty), on sko49 first.
- sko49's estimate rests on probes, not complete subtrees, in its tail: the factor here is for the
  188-232 short class only.

## Files

| | |
|---|---|
| `run.sh`, `queue.sh`, `jobs.txt`, `queue.log` | the runner, the four-at-a-time queue, the jobs in the order run |
| `par/` | the arms |
| `runs/NAME_ARM/` | `cmd.txt` (command, start, end, load), `log.txt`, `err.txt`, `dnn.txt`, `nodes.txt` |
| `tracelib.py`, `summarize.py` | the trace reader; the table, `--phases`, `--matched` |
| `fitcost.py`, `replay.py`, `replay_*.out` | the hand-over cost by gap; the offline replay and its outputs |
| `wcsweep.py`, `extrap.py` | the `wc` table; the estimate on the laptop's subtrees |
| `lscpu.txt` | the box |
| `mkreadme.sh`, `README.tmpl.md` | this file's assembly |
