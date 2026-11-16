# tai50b: the size of the tree with the cluster-region cut phase (bound_lab session `tai50b_size`, Nov 16 2026)

Question (from the bound lab): how big is tai50b's tree with the cluster-region cut phase of
`../regions/README.md`, in nodes and CPU-seconds, and with what settings? **[m]** measured here (log
named), **[e]** estimated. Incumbent everywhere: 458,821,518 = the best known value + 1 (cutoff
458,821,517.05), so every run has to prove that nothing better than the best known solution exists.

## Answer

**tai50b with the cluster-region cut phase is a run of about 100 core-hours on this machine: roughly 5,000 to
13,000 nodes and 2e5 to 6e5 CPU-seconds. Hundreds of core-hours at the outside, not thousands.** [e]

| estimate | nodes | CPU-s | core-hours here | what it is |
|---|---|---|---|---|
| E1 | 4,573 (2,473 to 7,135) | 228,101 | **63** (35 to 105) | depths 0-1 measured, 18 of the 50 depth-1 nodes; below them 275 Knuth dives, expanded by the sampling design alone |
| E2 | 12,023 (5,786 to 19,493) | 420,631 | **117** (75 to 164) | the same runs, with the two depth-1 children that hold a near-optimal solution counted as copies of the on-path child's subtree |
| E3 | 13,189 | 601,790 | **167** | a model of which nodes cannot fathom (those holding a solution within the gap measured on the best known solution's path), priced with measured node costs; 60 core-hours at half that gap |
| far case | 66,159 | 3,615,595 | 1,004 | E3 at twice the measured gap, which the nodes bounded here contradict (1 of 199 nodes that far out branched) |

- Intervals are 5th to 95th bootstrap percentiles and are too narrow: they do not know that the levels that
  carry the cost have few dives (see "What could make it wrong"). **My range is 60 to 250 core-hours**, with E2
  as the number I would plan on. On a c7a core [e]: about 1.4 times these figures, if this box's 0.70 of a c7a
  per matched tho40 node carries over (section 1).
- The old estimate (no cut phase) was a floor of 1e11 CPU-s. This one is 5 to 6 orders of magnitude lower.
- **Where the cost is** (E2 by depth, section 7): depth 5 36%, depth 2 20%, depth 6 17%, depth 4 12%, depth 1
  10%, the root 1%; everything from depth 7 down 3%. The DNN iterations are the whole cost (0.3 to 0.7 s each).
- **The shape of the tree.** The root ends 5.27M (1.15%) under the best known value. 11 of 17 sampled
  depth-1 children fathom there. Under the ones that branch, all 93 depth-2 node bounds of the dives (79 distinct nodes) fathomed.
  What does not fathom is the path of the best known solution and the few nodes that hold a solution nearly
  as good: the bound on that path is 3.2M short at depth 1, 0.63M at depth 5, 0.11M at depth 8, 5k at depth
  13, and a node holding a solution within about that gap has to be branched. The landscape is flat (150
  2-swaps of the best known solution are within 1M of it), so this "web" is a few dozen open nodes wide at
  depths 4 and 5 and thin again by depth 7.

**Settings** (`par/A.par`: `experiments/estimates/tai50b/tai50b.par` with these changed or added; the root
separately with `cut.stall = 0` and `cut.cp_iters = 6000`, `par/root_long.par`):

```
dnn.depths = 0:48        # the root too: it carries the cut pool and the state its children start from
dnn.giveup = 100         # a pooled first solve hands over to the phase early
pool.on = 1
cut.regions = experiments/bound_lab/regions/regions/tai50b.comb.txt     # 61 location-side regions (absolute path in the .par)
cut.side = 0
cut.cp = 100:1000
cut.warmy = 1
cut.zone = 1             # a phase at every DNN miss
cut.rounds = 12          # x cut.iters = 2,400 phase iterations at most; I would raise it to 20 (section 9) [e]
cut.iters = 200
cut.freeze = 0
cut.stall = 3
```

Add `dnn.beta_sched = w=0.36,wc=2.88` (the coordinator's finding, checked on four depth-1 nodes after the
estimates were made, section 9b): fathoming phases took 1.3x to 1.8x fewer iterations. The estimates were made
without it, so they should be high on that account.

Needs `../regions/regions.patch`. These are the settings every number here was measured under, chosen once from
the tho40 block and the regions README, not tuned: no A/B was run (section 9 says what the runs suggest).

**What it rests on** [m]

- The root: **453,552,676 certified** after 470 + 6,000 iterations (4,341 CPU-s), a new best root bound, and
  not yet the plateau (section 3).
- 18 of the 50 depth-1 nodes, warm from the root with its pool: the on-path child and a random 17 of the other
  49 (section 4).
- The path of the best known solution bounded node by node to depth 13, and the path of a second near-optimal
  solution (under `31:32`) to depth 4, which behaves the same (section 5).
- 275 dives from saved states, stratified by the path's levels, 359 node bounds in all with section 9b's (section 6); a model of
  the open nodes checked against the 235 distinct nodes bounded (section 7).
- Validity: 20 node bounds on the path of the best known solution, none fathomed, none above 458,821,517;
  tai30b solved from an incumbent one above its optimum in 363 nodes, as in the regions session (sections 2, 8).
- 49,000 CPU-s of node bounds in total, on 4 cores.

**What could make it wrong**

1. **The web is thinly sampled.** Levels 4 and 5 carry about half of E2 and have 8 and 5 dives, each stopped
   by its wall limit with a long dive in flight; the dive estimator weights a rare child by 100 or more, so one
   dive moves those levels by tens of thousands of CPU-s. E1 and E2 can be off by 2x either way from this alone.
2. **The twins are assumed, not measured, below depth 4.** E2 gives `31:0` (never bounded) and `31:32` the
   on-path child's subtree. E3 assumes a twin's gap is the spine's, that 2-swap descent finds the best solution
   inside a node, and that every node branches on the spine's line. More near-optimal solutions than the
   descent finds would widen the web; the nodes bounded so far agree with the rule at scale 1 (section 7).
3. **32 depth-1 nodes were not bounded.** By the landscape only two of them (`31:0`, `31:42`) are near, and E2
   and E3 count them; a far node that stalls anyway costs 49 quick children (4 of the 15 sampled far nodes did).
4. **Below depth 13 on the path nothing was run.** The estimates extrapolate the path to depth 45 at 13,000
   CPU-s per path (E2 counts three such paths: 9% of it), on the pattern that the bound there is within 5k and its children fathom in seconds.
   tai30b's path behaved that way. If the open nodes near the bottom multiply (many solutions within a few
   thousand of the best known: the nearest 2-swap is 1,534 away), the bottom of the tree grows; those nodes
   are small (m under 37) and cheap, so I expect nodes more than CPU.
5. **The cut code's validity is the regions session's**, plus the path checks here. A wrong cut that fathoms
   off the path would make every number here too small and the run's answer wrong. tai35b and tai40b solved
   end to end are the test, not yet run.
6. If 458,821,517 is not optimal, the search finds the better value and the tree is smaller.
7. CPU-seconds are this machine's (263 ms per plain iteration at m = 48); nodes and iterations carry over.
8. **The other session's samples** put the tree at 2,000 to 7,000 c7a core-hours if the fringe of near-optimal
   nodes compounds for about three levels; the dives and the model here say it compounds for two levels and
   stops, on thin evidence (section 9b).
9. In the other direction: a longer root, a longer phase limit and any tuning would shrink it (section 9).

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

| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |
|---|---|---|---|---|---|---|---|---|
| `31:7` | yes | 110 | 5,269,086 | 690 | 3,222,232 | stalled | branched | 451 |
| `31:4` |  | 110 | 4,358,673 | 1260 | 1,477,383 | stalled | branched | 776 |
| `31:27` |  | 110 | 2,849,001 | 681 | -678 | fathomed | fathomed (phase) | 457 |
| `31:20` |  | 70 | -87,228 |  |  | - | fathomed (first solve) | 38 |
| `31:46` |  | 110 | 4,023,989 | 1830 | -4,074 | fathomed | fathomed (phase) | 1104 |
| `31:5` |  | 110 | 4,011,441 | 2400 | 12,000 | limit | branched | 1504 |
| `31:36` |  | 110 | 3,969,462 | 1658 | -773 | fathomed | fathomed (phase) | 991 |
| `31:11` |  | 110 | 3,820,529 | 2400 | 347,520 | limit | branched | 1479 |
| `31:25` |  | 110 | 2,780,829 | 522 | -658 | fathomed | fathomed (phase) | 390 |
| `31:6` |  | 930 | -99 |  |  | - | fathomed (first solve) | 491 |
| `31:47` |  | 110 | 4,180,529 | 1643 | -12 | fathomed | fathomed (phase) | 998 |
| `31:30` |  | 110 | 3,038,990 | 887 | -2,683 | fathomed | fathomed (phase) | 589 |
| `31:15` |  | 110 | 3,823,147 | 1440 | 1,106,730 | stalled | branched | 856 |
| `31:43` |  | 110 | 4,501,903 | 800 | 2,679,372 | stalled | branched | 552 |
| `31:45` |  | 110 | 2,147,814 | 266 | -2,699 | fathomed | fathomed (phase) | 219 |
| `31:10` |  | 110 | 1,822,066 | 130 | -27,040 | fathomed | fathomed (phase) | 130 |
| `31:32` |  | 110 | 4,998,939 | 790 | 2,939,928 | stalled | branched | 501 |
| `31:38` |  | 110 | 3,330,029 | 1990 | -257 | fathomed | fathomed (phase) | 1233 |

- **11 of the 17 fathom at depth 1** (2 in the first solve, 9 in the phase after 130 to 1,990 iterations),
  6 branch; with `31:7` that is 7 of 18 branching. Mean 709 CPU-s per depth-1 node.
- Of the 6 that branch, **two hit the 2,400-iteration phase limit 12,000 and 347,520 short** (`31:5`,
  `31:11`: `31:5` was gaining 30k per 80 iterations and would have fathomed about 50 iterations later), two
  stalled 1.1M and 1.5M short (`31:15`, `31:4`), and two stalled 2.7M and 2.9M short (`31:43`, `31:32`).
  The last two hold solutions within 2M of the best known value (section 5) and cannot fathom at depth 1;
  the other four could, with a longer phase.

**From the short root, the first six children in file order** (not a random sample; stopped when the long
root arrived) [m, `runs/s1_*/dnn.txt`]:

| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |
|---|---|---|---|---|---|---|---|---|
| `31:0` |  | 110 | 12,051,521 | 1160 | 5,296,925 | stalled | branched | 721 |
| `31:3` |  | 110 | 9,456,713 | 2089 | -4,271 | fathomed | fathomed (phase) | 1291 |
| `31:1` |  | 110 | 5,061,353 | 282 | -4,827 | fathomed | fathomed (phase) | 215 |
| `31:4` |  | 110 | 11,400,049 | 1340 | 4,470,036 | stalled | branched | 853 |
| `31:7` | yes | 110 | 12,561,721 | 990 | 6,192,043 | stalled | branched | 623 |
| `31:2` |  | 110 | 11,262,032 | 1340 | 4,482,149 | stalled | branched | 871 |

The two nodes in both sets: `31:7` ends 6.19M short from the short root and 3.22M from the long one; `31:4`
4.47M and 1.48M. Under the short root's branched children every sampled depth-2 node still fathomed (18 of
18 dives, mean 313 iterations and 152 CPU-s, `runs/reh`, `logs/dives_reh.txt`), against 23 to 110 CPU-s under
the long root (section 6).

## 5. Why the tree is small, and where it is not: the path of the best known solution and its twins

**The spine** [m, `runs/spine_D/dnn.txt`, `spine.sh`]: the node on the path of `QAPSOLN/tai50b.sln` at each
depth, each warm from its parent with its pool. None may fathom; the question is how far below
458,821,517 the bound stays, because a node holding a solution `delta` above the best known value cannot
fathom while the gap on the path is larger than `delta`.

| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |
|---|---|---|---|---|---|---|---|---|
| `31:7` | yes | 110 | 5,269,086 | 690 | 3,222,232 | stalled | branched | 451 |
| `31:7,33:32` | yes | 110 | 3,525,708 | 600 | 2,652,428 | stalled | branched | 330 |
| `31:7,32:0,33:32` | yes | 110 | 3,199,122 | 650 | 2,185,798 | stalled | branched | 319 |
| `31:7,32:0,33:32,45:34` | yes | 110 | 3,569,687 | 1240 | 1,029,797 | stalled | branched | 550 |
| `31:7,32:0,33:32,40:39,45:34` | yes | 110 | 1,643,719 | 960 | 626,571 | stalled | branched | 370 |
| `30:2,31:7,32:0,33:32,40:39,45:34` | yes | 230 | 733,184 | 1670 | 199,906 | stalled | branched | 573 |
| `29:31,30:2,31:7,32:0,33:32,40:39,45:34` | yes | 240 | 556,039 | 1630 | 112,268 | stalled | branched | 485 |
| `21:20,29:31,30:2,31:7,32:0,33:32,40:39,45:34` | yes | 260 | 356,033 | 1180 | 107,349 | stalled | branched | 330 |
| `21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34` | yes | 230 | 245,752 | 1430 | 67,660 | stalled | branched | 331 |
| `21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34,49:41` | yes | 260 | 283,350 | 1790 | 22,865 | stalled | branched | 388 |
| `21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41` | yes | 230 | 173,421 | 1590 | 17,189 | stalled | branched | 338 |
| `21:20,22:33,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41` | yes | 250 | 35,128 | 1040 | 11,084 | stalled | branched | 226 |
| `1:23,21:20,22:33,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41` | yes | 230 | 57,575 | 1800 | 5,428 | stalled | branched | 466 |

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

| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |
|---|---|---|---|---|---|---|---|---|
| `31:32,33:7` |  | 110 | 3,332,829 | 640 | 2,296,390 | stalled | branched | 346 |
| `31:32,32:0,33:7` |  | 110 | 2,736,089 | 660 | 1,812,881 | stalled | branched | 322 |
| `31:32,32:0,33:7,45:34` |  | 110 | 3,150,868 | 2020 | 458,995 | stalled | branched | 864 |

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

| set | start nodes (depth) | processes | dives | node bounds: fathomed by the first solve / by the phase / branched | deepest dive | mean CPU-s per node bound | estimate of the set's subtrees: nodes | CPU-s |
|---|---|---|---|---|---|---|---|---|
| `LA_off2` | 48 (2) | 1 | 15 | 11 / 4 / 0 | 2 | 78 | 49 ± 1 | 3,729 ± 1,305 |
| `LA_off3` | 47 (3) | 1 | 25 | 23 / 2 / 0 | 3 | 47 | 47 ± 1 | 2,184 ± 571 |
| `LA_off4` | 46 (4) | 1 | 8 | 4 / 4 / 2 | 5 | 267 | 1,218 ± 825 | 43,272 ± 26,579 |
| `LA_off5` | 45 (5) | 1 | 5 | 3 / 2 / 2 | 6 | 308 | 1,348 ± 1,280 | 34,402 ± 21,498 |
| `LA_off6` | 44 (6) | 2 | 19 | 15 / 4 / 0 | 6 | 72 | 43 ± 1 | 3,216 ± 1,060 |
| `LA_off7` | 43 (7) | 1 | 26 | 24 / 2 / 0 | 7 | 32 | 43 ± 1 | 1,277 ± 421 |
| `LA_off8` | 42 (8) | 2 | 21 | 20 / 1 / 0 | 8 | 9 | 42 ± 1 | 362 ± 132 |
| `LA_tail8` | 42 (9) | 3 | 85 | 81 / 4 / 3 | 10 | 23 | 97 ± 32 | 1,305 ± 472 |
| `LB` | 196 (2) | 3 | 8 | 5 / 3 / 0 | 2 | 195 | 193 ± 4 | 37,105 ± 17,520 |
| `LB2` | 98 (2) | 1 | 25 | 25 / 0 / 0 | 2 | 23 | 99 ± 2 | 2,287 ± 208 |
| `LB_43` | 49 (2) | 1 | 8 | 5 / 3 / 0 | 2 | 220 | 51 ± 2 | 11,300 ± 4,878 |
| `LB_near` | 98 (2) | 1 | 30 | 30 / 0 / 0 | 2 | 35 | 97 ± 1 | 3,445 ± 521 |

Every dive and every node bound of every set: `logs/dives_SET.txt` (one line per dive:
`depth:outcome/iterations/final shortfall in M`), `runs/SET/sSEED/{dives.txt,dnn.txt,est.out}`.

Stopping, process by process (a process stopped by the clock loses its dive in flight, which under-samples
long dives; a process that finishes its planned count does not):

| set | seed | planned dives | finished | ended |
|---|---|---|---|---|
| `LA_off2` | 212 | by the clock | 15 | its wall limit |
| `LA_off3` | 213 | by the clock | 25 | its wall limit |
| `LA_off4` | 224 | 12 | 8 | its wall limit |
| `LA_off5` | 235 | 12 | 5 | its wall limit |
| `LA_off6` | 216 | 15 | 15 | its count |
| `LA_off6` | 226 | 12 | 4 | its wall limit |
| `LA_off7` | 217 | by the clock | 26 | its wall limit |
| `LA_off8` | 218 | 15 | 9 | stopped by hand |
| `LA_off8` | 228 | 12 | 12 | its count |
| `LA_tail8` | 229 | by the clock | 42 | stopped by hand |
| `LA_tail8` | 239 | by the clock | 23 | stopped by hand |
| `LA_tail8` | 249 | by the clock | 20 | stopped by hand |
| `LB` | 302 | by the clock | 6 | stopped by hand |
| `LB` | 303 | by the clock | 1 | stopped by hand |
| `LB` | 304 | by the clock | 1 | stopped by hand |
| `LB2` | 322 | by the clock | 25 | its wall limit |
| `LB_43` | 312 | 20 | 8 | stopped by hand |
| `LB_near` | 311 | 30 | 30 | its count |

Two first attempts at levels 4 and 5 that the clock ended after 2 and 6 dives (`runs/trunc/`), and 7 dives
from all 49 children of `31:7` before the sets were stratified (`runs/LA0`), are not used in the estimate; their
node bounds are in the path check.

What the dives show: off the spine nothing goes deep. Of the 275 dives, 7 met a node that branched. Four
did so off the spine, at depths 4 and 5, each at a near miss that ran to the 2,400-iteration phase limit 0.05M
to 0.20M short and whose child then fathomed. The other three are `LA_tail8` dives that drew the spine's own
node at depth 9 and then an off-path child, which fathomed. The longest fathoming phases are at depths 2 and 6 (up to 1,700 and 980 iterations). So
the cost below depth 1 is: every child of every open node bounded once (20 to 100 CPU-s each, more at
depths 4 to 6), plus the open nodes themselves (300 to 600 CPU-s on the spine, about 1,000 for a near miss).

## 7. The estimate [e] (`estimate.py`, `logs/estimate.txt`)

Three estimates, from the same runs:

- **E1, the sampling design alone**: depths 0 and 1 measured (the root, the on-path child, 17 random others
  expanded 49/17), everything below from the dives: the on-path child's subtree by spine level, the sampled
  branched nodes' subtrees from `LB*`.
- **E2, post-stratified by the landscape**: the two depth-1 children that hold a solution within 1M of the best
  known value (`31:0`, not sampled, and `31:32`) are counted as copies of the on-path child's subtree (the twin
  probe of section 5 supports it for `31:32`); `31:43` and `31:42` get `31:43`'s dives; the other 45 are
  expanded 45/15 from the 15 sampled. E2 is above E1 because the dives under `31:32` (few, and blind to its
  path) saw none of what the spine-stratified dives see under `31:7`.
- **E3, the web model priced with measured costs** (`web.py`, `webcost.py`, `logs/web.txt`, `logs/webcost.txt`,
  section 5's rule: a node is open iff it holds a solution within g(depth) of the best known value). The rule
  against the 235 distinct nodes bounded here (`webcheck.py`, `logs/webcheck.txt`):

  | delta / g(depth) | nodes | branched | of those: at the phase limit or stalled within 1.5M (near misses) | fathomed | mean iterations of the fathomed |
  |---|---|---|---|---|---|
  | 0 to 0.5 | 13 | 13 | 7 | 0 | - |
  | 0.5 to 1 | 7 | 4 | 3 | 3 | 1024 |
  | 1 to 1.5 | 5 | 2 | 2 | 3 | 1066 |
  | 1.5 to 2 | 11 | 2 | 2 | 9 | 170 |
  | 2 to 3 | 28 | 0 | 0 | 28 | 302 |
  | 3 to 5 | 36 | 0 | 0 | 36 | 286 |
  | 5 to inf | 135 | 0 | 0 | 135 | 109 |

  Below half the spine's gap every node branched; between 0.5 and 1.5 about half do, most of them near misses
  at the phase limit; above twice the gap one node in 199 did (a depth-1 near miss). So gap scale 1.0 is the
  model to read, 0.5 a floor, and 2.0 is contradicted by these nodes (it is in the table as the far case).

```
cost of an open node by depth (the spine's) [m]: 1: 451, 2: 330, 3: 319, 4: 550, 5: 370, 6: 573, 7: 485, 8: 330, 9: 331, 10: 388, 11: 338, 12: 226, 13: 466
mean cost of a fathomed node by depth [m] (node bounds): 1: 604 (11), 2: 70 (93), 3: 47 (25), 4: 99 (6), 5: 40 (5), 6: 74 (21), 7: 32 (26), 8: 9 (21)
below an open depth-8 node (LA_tail8, 85 dives, + the spine below depth 8) [e]: 948 nodes, 14436 CPU-s
children of the depth-1 near misses: 12 parents x 49 x 30 CPU-s (55 node bounds) [e]

gap scale 1.0: **13189 nodes, 601790 CPU-s = 167 core-hours**

| depth | nodes | open | CPU-s |
|---|---|---|---|
| 0 | 1 | 1 | 4341 |
| 1 | 50 | 5 | 29416 |
| 2 | 245 | 14 | 20700 |
| 3 | 672 | 11 | 34732 |
| 4 | 517 | 43 | 70500 |
| 5 | 1978 | 65 | 100950 |
| 6 | 2925 | 7 | 221324 |
| 7 | 308 | 5 | 12015 |
| 8 | 215 | 6 | 3811 |
| 2, under the 12 near misses | 588 | 0 | 17384 |
| 9 and below | 5690 | | 86618 |

gap scale 0.5: **5285 nodes, 215873 CPU-s = 60 core-hours**

| depth | nodes | open | CPU-s |
|---|---|---|---|
| 0 | 1 | 1 | 4341 |
| 1 | 50 | 3 | 29722 |
| 2 | 147 | 6 | 11794 |
| 3 | 288 | 6 | 15235 |
| 4 | 282 | 10 | 32378 |
| 5 | 460 | 10 | 21791 |
| 6 | 450 | 2 | 34510 |
| 7 | 88 | 2 | 3692 |
| 8 | 86 | 3 | 1717 |
| 2, under the 12 near misses | 588 | 0 | 17384 |
| 9 and below | 2845 | | 43309 |

gap scale 2.0: **66159 nodes, 3615595 CPU-s = 1004 core-hours**

| depth | nodes | open | CPU-s |
|---|---|---|---|
| 0 | 1 | 1 | 4341 |
| 1 | 50 | 20 | 27121 |
| 2 | 980 | 116 | 98444 |
| 3 | 5568 | 124 | 296717 |
| 4 | 5828 | 177 | 655751 |
| 5 | 8142 | 483 | 486533 |
| 6 | 21735 | 38 | 1637615 |
| 7 | 1672 | 17 | 60621 |
| 8 | 731 | 22 | 13470 |
| 2, under the 12 near misses | 588 | 0 | 17384 |
| 9 and below | 20864 | | 317599 |
```

`estimate.py`'s output (E1, E2, the on-path child's subtree T_A, and E2 by depth):

```
root: 6470 iterations, 4,341 CPU-s, certified 453,552,676
depth 1 [m]: 31:7 branched (451 CPU-s); B sample 17 of 49: 11 fathomed, 6 branched; mean CPU-s 724
spine [m]: depth: shortfall {2: 2652428, 3: 2185798, 4: 1029797, 5: 626571, 6: 199906, 7: 112268, 8: 107349, 9: 67660, 10: 22865, 11: 17189, 12: 11084, 13: 5428}
dives per set: {'LA_off2': 15, 'LA_off3': 25, 'LA_off4': 8, 'LA_off5': 5, 'LA_off6': 19, 'LA_off7': 26, 'LA_off8': 21, 'LA_tail8': 85, 'LB': 8, 'LB2': 25, 'LB_43': 8, 'LB_near': 30} | sets without dives: []
deepest dives: {'LA_off2': 2, 'LA_off3': 3, 'LA_off4': 5, 'LA_off5': 6, 'LA_off6': 6, 'LA_off7': 7, 'LA_off8': 8, 'LA_tail8': 10, 'LB': 2, 'LB2': 2, 'LB_43': 2, 'LB_near': 2}
the spine below depth 8 [e] (spine_below8): 851 nodes, 13,132 CPU-s

T_A (below 31:7) [e]: 3,745 nodes, 105,838 CPU-s; by set (nodes, CPU-s): {'LA_off2': (49, 3729), 'LA_off3': (47, 2184), 'LA_off4': (1218, 43272), 'LA_off5': (1348, 34402), 'LA_off6': (43, 3216), 'LA_off7': (43, 1277), 'LA_off8': (42, 362), 'LA_tail8': (97, 1305)}
T_i of the branched B parents [e] (nodes, CPU-s): {'31:4': (54, 1588), '31:5': (47, 1125), '31:11': (53, 1162), '31:15': (43, 1686), '31:43': (73, 22887), '31:32': (0, 0)}

E1 nodes: 4,573 nodes  (bootstrap 5%-95%: 2,473 to 7,135)
E1 CPU: 228,101 CPU-s  (bootstrap 5%-95%: 127,719 to 376,809)  = 63 core-hours (35 to 105)
E2 nodes: 12,023 nodes  (bootstrap 5%-95%: 5,786 to 19,493)
E2 CPU: 420,631 CPU-s  (bootstrap 5%-95%: 268,303 to 590,126)  = 117 core-hours (75 to 164)
T_A nodes: 3,745 nodes  (bootstrap 5%-95%: 1,672 to 6,226)
T_A CPU: 105,838 CPU-s  (bootstrap 5%-95%: 56,384 to 161,910)  = 29 core-hours (16 to 45)

E2 by depth [e] (depths 0, 1 and the spine measured):

| depth | nodes | CPU-s | share of CPU |
|---|---|---|---|
| 0 | 1 | 4,341 | 1.1% |
| 1 | 50 | 36,321 | 9.5% |
| 2 | 654 | 74,635 | 19.6% |
| 3 | 141 | 7,510 | 2.0% |
| 4 | 135 | 46,558 | 12.2% |
| 5 | 3,651 | 136,438 | 35.7% |
| 6 | 4,043 | 64,626 | 16.9% |
| 7 | 129 | 5,288 | 1.4% |
| 8 | 130 | 2,057 | 0.5% |
| 9 | 122 | 2,802 | 0.7% |
| 10 | 168 | 1,112 | 0.3% |
| all | 9,224 | 381,688 | |
```

## 8. Validity [m] (`validity.py`, `logs/validity.txt`)

```
359 traced node bounds in 53 traces; 20 on the path of the best known solution
  depth 0 -: certified 453,552,676.2916  (5,268,841 under)  branched  [runs/root_long/dnn.txt]
  depth 0 -: certified 445,306,477.0661  (13,515,040 under)  branched  [runs/root_short/dnn.txt]
  depth 1 31:7: certified 455,599,284.8629  (3,222,232 under)  branched  [runs/L1_0/dnn.txt]
  depth 1 31:7: certified 452,629,473.6287  (6,192,043 under)  branched  [runs/s1_1/dnn.txt]
  depth 1 31:7: certified 456,310,021.6095  (2,511,495 under)  branched  [runs/wc_7/dnn.txt]
  depth 2 31:7,33:32: certified 456,169,088.9166  (2,652,428 under)  branched  [runs/spine_2/dnn.txt]
  depth 3 31:7,32:0,33:32: certified 456,635,718.7546  (2,185,798 under)  branched  [runs/spine_3/dnn.txt]
  depth 4 31:7,32:0,33:32,45:34: certified 457,791,720.1686  (1,029,797 under)  branched  [runs/spine_4/dnn.txt]
  depth 5 31:7,32:0,33:32,40:39,45:34: certified 458,194,945.8105  (626,571 under)  branched  [runs/spine_5/dnn.txt]
  depth 6 30:2,31:7,32:0,33:32,40:39,45:34: certified 458,621,611.2378  (199,906 under)  branched  [runs/spine_6/dnn.txt]
  depth 7 29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,709,249.2532  (112,268 under)  branched  [runs/spine_7/dnn.txt]
  depth 8 21:20,29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,714,168.5303  (107,348 under)  branched  [runs/spine_8/dnn.txt]
  depth 9 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,753,856.6484  (67,660 under)  branched  [runs/LA_tail8/s239/dnn.txt]
  depth 9 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,753,856.6484  (67,660 under)  branched  [runs/LA_tail8/s249/dnn.txt]
  depth 9 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,753,856.6484  (67,660 under)  branched  [runs/LA_tail8/s249/dnn.txt]
  depth 9 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34: certified 458,753,856.6484  (67,660 under)  branched  [runs/spine_9/dnn.txt]
  depth 10 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,45:34,49:41: certified 458,798,651.6812  (22,865 under)  branched  [runs/spine_10/dnn.txt]
  depth 11 21:20,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41: certified 458,804,328.4075  (17,189 under)  branched  [runs/spine_11/dnn.txt]
  depth 12 21:20,22:33,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41: certified 458,810,432.5647  (11,084 under)  branched  [runs/spine_12/dnn.txt]
  depth 13 1:23,21:20,22:33,27:10,29:31,30:2,31:7,32:0,33:32,40:39,44:1,45:34,49:41: certified 458,816,089.3467  (5,428 under)  branched  [runs/spine_13/dnn.txt]
violations: 0
largest certified value at an unfathomed node anywhere: 458,816,089.3467
largest certified value at a fathomed node: 461,716,998.5014 (a fathom stops at the first certificate above the cutoff)
```

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

- **The root's plateau** (section 3): 6,470 iterations, still climbing, no curve.
- **A separate random sample of 40 or more depth-2 nodes, and depth 3 below them**: the dives' start nodes are
  that sample instead (the table of section 6 counts them by depth; `logs/dives_SET.txt` has each one).
- **Thirty complete dives through the web.** The dives number well over 200, but most start and end at
  one fathoming node; the ones that matter (a branch at depths 4 to 6, section 6) are a handful, and levels
  4, 5, 6 and 8 have few dives. E1 and E2's intervals understate that.
- **Any tuning A/B**: one settings arm (`par/A.par`) throughout. The phase limit, the stall rule, the region
  set, `pool.cp`, `dnn.amx` and the root's length are all open (section 9).
- **The 32 unsampled depth-1 nodes**, among them `31:0` (a twin by the landscape) and `31:42`.
- **The spine below depth 13 and any twin below depth 4**; the path check stops there too.
- **End-to-end validity on solved instances other than tai30b** (tai35b, tai40b with their own region
  files), which the regions session listed as the next check and which should come before a tai50b run.
- **The run itself.** At this size the next step is not a better estimate but the search: `bnbdrv` from
  the saved long root with `par/A.par` (phase limit raised), under the coordinator or `bench/parallel.py`.

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
