# tai50b: the size of the tree with the cluster-region cut phase (bound_lab session `tai50b_size`, Nov 16 2026)

Question (from the bound lab): how big is tai50b's tree with the cluster-region cut phase of
`../regions/README.md`, in nodes and CPU-seconds, and with what settings? **[m]** measured here (log
named), **[e]** estimated. Incumbent everywhere: 458,821,518 = the best known value + 1 (cutoff
458,821,517.05), so every run has to prove that nothing better than the best known solution exists.

## Answer

@ANSWER@

## 1. Setup, machine, calibration

- `git apply ../regions/regions.patch` (adds `cut.regions`), plus `est_resume.patch` (mine, 8 lines in
  `latest/estimate.C` and `latest/driver.C`): `bnbdrv -e -r STATE` takes the root bound from the state
  file, so dives can start from the open nodes of a state file and nothing above them is recomputed.
  Without it the resumed estimator scores children against a root bound of 0. Neither patch is committed
  to `latest/`; both are files here and in `../regions/`.
- Machine [m, `logs/machine.txt`]: 4 cores, Intel Xeon @ 2.30 GHz (family 6 model 173, AVX-512, AMX), KVM.
  One `-T 1` process per core (`taskset`), OpenBLAS single-threaded, all four cores busy nearly all the time.
- Calibration [m, `logs/calib_m48.txt`]: `dnntest -p tai50b.par -x 7:40,40:39 -i 100` (m = 48, plain DNN,
  cold) **263 ms per iteration** with the other three cores busy; with `dnn.amx = 1` 210 ms
  (`logs/calib_m48_amx1.txt`; not used, so that CPU converts to machines without AMX).
- In the search the iterations cost more than that plain figure, because a pooled first solve carries its
  parent's cuts and a phase projects onto its pool [m, all traces under `runs/`]:

  | depth (m) | first solve, s per iteration | phase, s per iteration |
  |---|---|---|
  | 0 (50) | 0.52 | 0.69 |
  | 1 (49) | 0.51 | 0.59 |
  | 2 (48) | 0.43 | 0.52 |
  | 3-5 (47-45) | 0.28-0.36 | 0.35-0.44 |

  To convert: scale by your machine's `dnntest` figure for the command above against 263 ms. [e] The
  tho40 study measured this box at 0.70 of a c7a core's CPU per matched node (`../tho40_search/README.md`).

## 2. The build, checked on tai30b [m]

`bnbdrv -p par/tai30b_cut.par -u 637117114 ../qaplib/tai30b.dat` (the regions README's run, one piece):
**363 nodes, 1,432 CPU-s, optimum 637,117,113 found and proved**; depth 1: 29 of 30 fathomed (8 by the first
solve, 21 by the phase), then one path [`logs/tai30b_check.log`, `.trace`]. The same node count as the
regions session. The pipeline used below (root bounded with the phase and saved, dives resumed from the
state) was rehearsed on tai30b first [`runs/t30_root`, `runs/t30_dives`].

## 3. The root [m]

`par/root_long.par` = `tai50b.par` with `dnn.depths = 0:48`, `dnn.giveup = 100` and the cut block, the phase
with `cut.stall = 0` and `cut.cp_iters = 6000`; `-d 0 -o root.state` saves the 50 children with the root's
state and pool (`run_root.sh`). `par/root_short.par`: the same with 1,200 phase iterations.

| root | first solve | phase iterations | certified | under the best known value | CPU-s | log |
|---|---|---|---|---|---|---|
| short | 470 it, 417,845,879 | 1,200 | 445,306,477 | 13,515,040 (2.95%) | 1,077 | `runs/root_short/dnn.txt` |
| **long** | 470 it, 417,845,879 | 6,000 | **453,552,676** | **5,268,841 (1.15%)** | 4,341 | `runs/root_long/dnn.txt` |

- **453,552,676 is a new best root bound for tai50b** (QAPLIB's 431,090,700; the regions session's
  449,469,933 at 1,900 iterations).
- **It is not the plateau, and I cannot give the plateau.** With `cut.stall = 0` the solver records no
  phase history (my mistake: the history is kept only for the stall rule), so the root's curve is three
  points: 445.31M at 1,200 phase iterations, 449.47M at 1,600 (regions session, another first solve),
  453.55M at 6,000. The average gain over the last 4,400 iterations is 93k per 100. The on-path child,
  which starts from the root's state and inherits its bound, gains 0.46M in its first 30 phase iterations
  and 2.05M in 690 (`runs/L1_0/dnn.txt`), part of it from the fix, so the root was still climbing. The
  on-path values of section 5 bound the plateau from the other side only loosely.
- A longer root pays for itself: the same depth-1 children start 3M to 7M closer from the long root
  (section 4), and 3,300 more root CPU-s is less than five depth-1 nodes.

## 4. Depth 1

The root branches on row 31 (50 children `31:x`). `par/A.par` (the settings block above, phase limit
2,400). Every node starts warm from the root's state with its pool (`run_state.sh`, `-r` on a piece of the
root's state file, `-d 1`).

**From the long root: the on-path child `31:7` and a random 17 of the other 49** (`pick_d1.py`, seed
20261116; 19 were drawn, the last two of one core's five were not run for time: `logs/d1_sample.txt`)
[m, `runs/L1_*/dnn.txt`, `logs/depth1_long.md`]:

@D1LONG@

- **11 of the 17 fathom at depth 1** (2 in the first solve, 9 in the phase after 130 to 1,990 iterations),
  6 branch; with `31:7` that is 7 of 18 branching. Mean 709 CPU-s per depth-1 node.
- Of the 6 that branch, **two hit the 2,400-iteration phase limit 12,000 and 347,520 short** (`31:5`,
  `31:11`: `31:5` was gaining 30k per 80 iterations and would have fathomed about 50 iterations later), two
  stalled 1.1M and 1.5M short (`31:15`, `31:4`), and two stalled 2.7M and 2.9M short (`31:43`, `31:32`).
  The last two hold solutions within 2M of the best known value (section 5) and cannot fathom at depth 1;
  the other four could, with a longer phase.

**From the short root, the first six children in file order** (not a random sample; stopped when the long
root arrived) [m, `runs/s1_*/dnn.txt`]:

@D1SHORT@

The two nodes in both sets: `31:7` ends 6.19M short from the short root and 3.22M from the long one; `31:4`
4.47M and 1.48M. Under the short root's branched children every sampled depth-2 node still fathomed (18 of
18 dives, mean 313 iterations and 152 CPU-s, `runs/reh`, `logs/dives_reh.txt`), against 23 to 110 CPU-s under
the long root (section 6).

## 5. Why the tree is small, and where it is not: the path of the best known solution and its twins

**The spine** [m, `runs/spine_D/dnn.txt`, `spine.sh`]: the node on the path of `QAPSOLN/tai50b.sln` at each
depth, each warm from its parent with its pool. None may fathom; the question is how far below
458,821,517 the bound stays, because a node holding a solution `delta` above the best known value cannot
fathom while the gap on the path is larger than `delta`.

@SPINE@

**The landscape** [m, `landscape.py`, `logs/landscape.txt`]: no 2-swap of the best known solution improves
it; 15 are within 100k of it, 52 within 300k, 150 within 1M, 250 within 2M (the nearest: +1,534, +9,549,
+23,497). For each depth-1 child `31:x`, the best solution found with p(31) = x (swap from the best known,
then 2-swap descent): `31:7` +0, `31:0` +147k, `31:32` +277k, `31:43` +1.96M, `31:42` +2.51M, `31:21` +3.73M,
`31:15` +4.04M, `31:38` +4.58M, and 5M to 41M for the rest. Against the on-path gap of 3.2M at depth 1:
**`31:0` and `31:32` are twins of the on-path child** (they cannot fathom until the gap is under 0.15M to
0.28M, depth 6 to 7 on the spine), `31:43` and `31:42` are half-twins (until about depth 3), and everything
else can fathom at depth 1 or just below. That is what the sample shows: `31:32` and `31:43` stall 2.9M and
2.7M short, and the four other branched nodes are near misses of the phase limit and the stall rule.
The spine's gap goes on falling below depth 8: 68k at depth 9, 23k, 17k, 11k and 5k at depths 10 to 13.

**The twin, probed** [m, `runs/twin32_D/dnn.txt`]: the path of the best solution found with `31:32`
(`logs/perm_31_32.txt`), depths 2 to 4:

@TWIN@

It tracks the spine (2.30M, 1.81M, 0.46M short at depths 2, 3, 4 against 2.65M, 2.19M, 1.03M), so I count
`31:0` and `31:32` as subtrees of the on-path child's size [e].

## 6. The dives [m for each dive, e for what they add up to]

`bnbdrv -e -r STATE` (`run_dive.sh`; `estimate.par`: start level = the depth of STATE's nodes, exponent 1.5;
child sampling `est.xbar = 1`, `est.xbar_mix = 0.3`; one seed per process). A random draw from all depth-2
nodes almost never meets the spine, so the dives are **stratified by it**: one set per spine level from the
off-path children of the spine node above (`LA_offD`, D = 2..8), one set from all children of the spine node at
depth 8 (`LA_tail8`), and the sets under the sampled off-path depth-1 nodes (`LB`: children of `31:4`,
`31:15`, `31:32`, `31:43`; `LB2`: of `31:5`, `31:11`; `LB_near`: of `31:4`, `31:15` alone; `LB_43`: of `31:43` alone).
The spine nodes themselves are measured (section 5).

@DIVES@

Every dive and every node bound of every set: `logs/dives_SET.txt` (one line per dive:
`depth:outcome/iterations/final shortfall in M`), `runs/SET/sSEED/{dives.txt,dnn.txt,est.out}`.

Stopping, process by process (a process stopped by the clock loses its dive in flight, which under-samples
long dives; a process that finishes its planned count does not):

@STOPS@

Two first attempts at levels 4 and 5 that the clock ended after 2 and 6 dives (`runs/trunc/`), and 7 dives
from all 49 children of `31:7` before the sets were stratified (`runs/LA0`), are not used in the estimate; their
node bounds are in the path check.

What the dives show: off the spine nothing goes deep. Of the @NDIVES@ dives, @NB1@ met a node that branched. Four
did so off the spine, at depths 4 and 5, each at a near miss that ran to the 2,400-iteration phase limit 0.05M
to 0.20M short and whose child then fathomed. The other three are `LA_tail8` dives that drew the spine's own
node at depth 9 and then an off-path child, which fathomed. The longest fathoming phases are at depths 2 and 6 (up to 1,700 and 980 iterations). So
the cost below depth 1 is: every child of every open node bounded once (20 to 100 CPU-s each, more at
depths 4 to 6), plus the open nodes themselves (300 to 600 CPU-s on the spine, about 1,000 for a near miss).

## 7. The estimate [e] (`estimate.py`, `logs/estimate.txt`)

@ESTIMATE@

## 8. Validity [m] (`validity.py`, `logs/validity.txt`)

@VALIDITY@

- tai30b from an incumbent one above its optimum: optimum found and proved, 363 nodes (section 2).
- Not checked here: the cut code itself (the regions session's brute-force set, `../regions/validity.py`,
  was not rerun), and nothing was checked below depth 8 on the path.

## 9. Settings: what the runs say [m unless marked]

- **A long root** (6,000 phase iterations or more, `cut.stall = 0`, `cut.cp_iters`): section 3. To see its
  curve, give it the history another way: `cut.stall` must be on for `QAP_DNN_TRACE_HIST` to record anything.
- **`dnn.giveup = 100`**: a pooled child that will not fathom in its first solve hands over to the phase
  after 110 iterations; the phase then gains from its first iterations (0.2M to 0.5M in the first 30).
  First solves that are going to fathom are not cut short (they ran 22 to 930 iterations).
- **The phase limit of 2,400 (`cut.rounds = 12` x `cut.iters = 200`) and the stall rule (`cut.stall = 3`) are
  about break-even as they stand.** At depth 1, `31:5` ended at the limit 12,000 short (gaining 30k per 80
  iterations: 50 more iterations would have fathomed it; instead 49 children) and `31:11` 348k short; `31:4`
  and `31:15` stalled 1.5M and 1.1M short while gaining 0.05M to 0.07M per 100, so about 2,000 more iterations
  each (1,100 CPU-s), against 49 children at about 30 CPU-s (1,500 CPU-s). At depths 4 and 5 the near misses
  end at the limit 0.05M to 0.20M short gaining 0.012M per 100, and their 45 children cost 700 to 2,000 CPU-s
  together. A limit of 4,000 (`cut.rounds = 20`) would save the `31:5` kind and cost little elsewhere [e: not rerun].
  The stall rule did its job on the nodes that cannot fathom (the spine, `31:32`, `31:43`: 600 to 800 iterations).
- Not tried: a smaller region set (the 16-region `tai50b.small.txt`; a phase iteration with the 61 regions
  costs 1.15 to 1.3 times a first-solve iteration), `pool.cp = 1`, `dnn.amx = 1` (0.80 of the plain
  iteration's CPU on this box), `dnn.gu_zone = 0`.

## 9b. After the report: the coordinator's messages, and one A/B [m]

Three messages from the coordinating session reached this session only after the report above was written
(they were queued at 01:39, 02:31 and 04:13 and read at 05:01). Every number above was measured without
them. What they change:

**The phase penalty, `dnn.beta_sched = w=0.36,wc=2.88`** (`../cpfast`), tested here on four depth-1 nodes
already bounded under `par/A.par`, same root state and pool (`par/A_wc.par`, `runs/wc_X/dnn.txt`):

| node | `A.par` (no `wc`): phase iterations, end, CPU-s | `wc=2.88`: phase iterations, end, CPU-s |
|---|---|---|
| `31:27` | 681, fathomed, 457 | **533, fathomed, 331** |
| `31:5` | 2,400, at the limit 12,000 short, branched, 1,504 | **1,323, fathomed, 763** |
| `31:43` (holds a solution 1.96M above the best known) | 800, stalled 2.68M short, 552 | 1,310, stalled **1.95M** short, 762 |
| `31:7` (on the path) | 690, stalled 3.22M short, 451 | 1,000, stalled **2.51M** short, 532 |

Fathoming phases take 1.3x to 1.8x fewer iterations and the nodes that cannot fathom end 0.7M closer, so
their children start closer. **Use `wc=2.88`; the estimates above were made without it and should be high for
that reason** [e: four nodes at one depth; the root, the web levels and `wc` = 1.0 and 5.76 were not run].

**The other session's samples** (`../regions/README.md`, "Follow-up (Nov 16)", read after the report): the same
picture of the top (one path through the 3-point cluster's locations 31, 33, 32, then a fringe of
same-cluster children from depth 4), and its figure for "every survivor's children die at once", 50 to 200
core-hours and 10,000 to 60,000 nodes, is this README's range. Its central figure (2,000 to 7,000 c7a
core-hours, order 1e6 nodes) assumes the fringe compounds for three or so levels at about 3 survivors per
survivor. What was measured here on that question:

- Off the path at depths 4 and 5, 4 of 13 dives met a node that branched (each at the phase limit, 0.05M to
  0.20M short); the child then drawn, by the DNN marginals, which favour the likely survivors, fathomed in
  all 4 (45 to 205 iterations in three, 681 in one). At depths 6, 7 and 8, 0 of 66 off-path nodes branched.
- The web model counts the survivors level by level from the solutions they hold: 43 and 65 open nodes at
  depths 4 and 5, then 7, 5, 6 at depths 6 to 8 (section 7): it compounds for two levels and stops, because
  the gap on the path falls from 0.63M at depth 5 to 0.20M at depth 6.
- The roots differ: its children start from a 3,000-iteration root (on-path nodes 4.82M, 3.88M, 3.10M short
  at depths 1 to 3), these from the 6,000-iteration one (3.22M, 2.65M, 2.19M), and it notes itself that a long
  phase on the path changes its children's fate (8 of 10 fathom at once instead of 2 of 10).

So I keep the answer, with this added to "What could make it wrong": four dives is all the direct evidence
against compounding at depths 4 and 5, and the two sessions' numbers differ by 20x on exactly that. Its
other settings (`dnn.gu_zone = 0`, `dnn.giveup_warm = 100`, `cut.stall = 0` with `cut.cp_iters = 3000`) were
not run here. Its report that `bnbdrv -e -r` printed a wrong level table is what `est_resume.patch` fixes
(without it the resumed estimator has no root bound).

## 10. Not done

@NOTDONE@

## Files

| | |
|---|---|
| `par/A.par` | the settings of every node below the root; `par/root_long.par`, `par/root_short.par` the roots; `par/tai30b_cut.par`, `par/t30_root.par` the tai30b check and rehearsal |
| `est_resume.patch` | `bnbdrv -e -r`: the root bound from the state file |
| `run_root.sh`, `run_state.sh`, `run_dive.sh`, `queue.sh`, `spine.sh`, `core*.sh`, `after_l1.sh` | the runners (root with state; a state file's nodes down to a depth; dives from a state file; a core's queue; the path of a permutation; what each core ran) |
| `splitstate.py`, `pick_d1.py` | state files: split, merge (renumbering the parents' warm lines), filter by depth, fix or path; the depth-1 sample |
| `nodes.py`, `dives.py`, `estimate.py`, `validity.py`, `landscape.py`, `web.py`, `webcheck.py`, `webcost.py`, `mkreport.sh`, `mkreadme.py`, `README.tmpl.md`, `answer.md`, `estimate_text.md`, `notdone.md` | tables, the estimates E1 and E2, the path check, the near-optimal solutions, the web model with its check and its cost (E3), this file and its hand-written parts |
| `runs/NAME/` | `cmd.txt`, `log.txt` or `est.out`, `dnn.txt` (one line per node bound, with the phase's history), `nodes.txt` (node trace, search runs), `dives.txt`. State files and the parents' DNN states (123 MB each) are not kept: `run_root.sh` and `run_state.sh` rewrite them |
| `logs/` | machine, calibration, tai30b check, the tables (`*.md`), `dives_SET.txt`, `pool_SET.txt` (per level, `../../tai50b_feasibility/estimate/pool.py`), `estimate.txt`, `validity.txt`, `landscape.txt`, the depth-1 sample |
