# nug30: does today's cut phase beat the record recipe? (bound_lab agent `nug30_phase`, Nov 16 2026)

Every arm below is the WHOLE nug30 search (`-u 6124`, `SYM/nug30.sym`, rule 9, DNN at depths 1-7),
4 workers of `bench/parallel.py` on one 4-core box (`lscpu.txt`: Xeon 2.1 GHz with AMX), one arm at a
time, nothing else running. **[m]** measured (the run's directory under `runs/` is named), **[e]** estimated.

## Answer

**Yes, by about a quarter, but not with the phase as the brief specified it, and not where tho40 gains.**

| | nodes | DNN iterations | CPU s | wall s |
|---|---|---|---|---|
| record recipe (`base`) [m] | 3,688 | 110,926 | 3,030 | 771 |
| best cut arm (`g5`) [m] | 1,674 | 72,147 | 2,197 | 564 |
| ratio | 2.20x fewer | 1.54x fewer | **1.38x less (-27%)** | **1.37x less (-27%)** |
| the same without `dnn.gu_zone = 0` (`g2`) [m] | 1,646 | 77,535 | 2,309 | 596 |
| ratio | 2.24x fewer | 1.43x fewer | 1.31x less (-24%) | 1.29x less (-23%) |

- A rerun of the baseline (`base2`) gave the same tree and iterations, 3,083 CPU s and 784 s wall: the
  CPU noise is under 2%.
- With `dnn.amx = 1` on both sides (this box has AMX; the record recipe leaves it off): `basea` 2,505 CPU s
  and 638 s wall, `g2a` 2,024 and 528: the same trees, 1.43x fewer iterations, **1.24x less CPU (-19%)**.
- **The brief's four zones do not win, except the narrowest, barely** [m]: `cut.zone = 0.003` 7% less CPU
  (1.10x fewer iterations), 0.01 **39% more** CPU, 0.03 19% more, 0.06 the same search as 0.03 (@@Z060@@).
  The Nov 1 verdict (+34% / +42%) still describes a phase that runs to its end at every give-up.
- **What wins is the pool, not the fathoming phase.** `g2` is the brief's block with `cut.zone = 0.03`
  plus an early stop for phases below depth 1:
  `cut.giveup = 100:5,300:2` and `cut.giveup_depth = 2` (a new key, three lines of `latest/bnb.C`
  and `settings.C`, default 0 = unchanged; see "The code change"). `g5` is `g2` with `dnn.gu_zone = 0`:
  a near miss goes to its phase at the give-up (200 iterations) instead of waiting 100 more for
  confirmation (median 320 in `g2`). In `g2`:
  - The nine depth-1 nodes each run a phase to the stall rule (630-1,240 iterations, 6,570 in all; none
    can fathom, they end 20-81 short). Their pools go to the 225 depth-2 nodes: 186 fathom in the first
    solve where the baseline fathoms 164, in 17,261 iterations where the baseline takes 23,770.
  - Below depth 1 a phase stops at iteration 100 unless it is within 5 of the cutoff (at 300, within 2),
    and hands its pool on. 23 of 76 fathom (median 100-280 iterations); the 53 stopped ones cost 100
    iterations each and their pooled children take about a quarter fewer iterations.
  - Phases are 20% of `g2`'s iterations (15,635 of 77,535) and 22% of `g5`'s (15,631 of 72,147).
- **Sanity: PASS.** Every arm ends at 6124 with no solution below it and no open node. The `g2` arm
  started from incumbent 6126 (`path_g2`) finds 6124 itself; all 28 nodes on the path to it branch, and no
  certified value on the path passes 6124.0000 (`pathcheck.out`).

## The arms [m]

`python3 summarize.py` (CPU s: the solver's CPU summed over the workers, as its RESULTS block; wall: `time.txt`):

@@TABLE@@

- `base`: `latest/PARAM/nug_deep_xbar9_all.par` with the record settings written out (`par/base.par`;
  the same as the solver's defaults with `-x 1:7:1:200:1`, `ship.bits = 32`). 3,688 nodes here; the
  record run on the M2 Pro had 3,705 (BLAS rounding at the cutoff).
- `z003`, `z010`, `z030`, `z060`: `base` plus the brief's block, `cut.zone` = 0.003, 0.01, 0.03, 0.06:
  a phase when the first solve ends within 18.4, 61.2, 183.7, 367.4 units of the cutoff (6122.05).
- `sko`: nug30's blocks and rule 9 with sko42's whole `[settings]` block (`cut.zone = 0.010`, beta 0.3,
  refresh 40, the two-phase beta schedule, give-up 300, `dnn.gu_zone = 0`, `dnn.amx = 1`, `dnn.depths = 1:40`).
- `d1`: `z030` with `cut.maxdepth = 1`: the nine depth-1 phases only; their pools are inherited down the tree.
- `g1`: `z030` with `cut.giveup = 100:5,300:2` at every depth (the depth-1 phases stop at 100 iterations too).
- `g2`: `g1` with `cut.giveup_depth = 2`. `g4`: `g2` with `cut.stall = 5`. `g5`: `g2` with `dnn.gu_zone = 0`.
- `basea`, `g2a`: `base` and `g2` with `dnn.amx = 1`. `g2s`: sko42's `dnn.*` block with `g2`'s phase.
- `path_g2`: `g2` from incumbent 6126 (the path check).

In the nodes column the tree shrinks 2-3.6x under every cut arm; the CPU does not follow, because the
nodes that disappear are the cheapest in the search.

## Why the zones lose and the pool wins

**1. On nug30 a give-up is cheap to branch on** [m, `base`]. Below depth 1 the baseline gives up at 131
nodes. Under the median one lie 400-530 iterations: 28 children that fathom in 19 iterations each from
the parent's state.

| depth | give-ups | first-solve shortfall: median (range) | iterations below each: median / mean | nodes below (mean) |
|---|---|---|---|---|
| 1 | 9 | 115.7 (61.6-121.8) | 15,161 / 12,125 | 409 |
| 2 | 61 | 23.3 (0.6-79.7) | 532 / 1,399 | 57 |
| 3 | 49 | 15.5 (1.7-52.5) | 464 / 859 | 37 |
| 4 | 17 | 8.9 (1.0-23.7) | 403 / 742 | 32 |
| 5 | 4 | 3.7 (0.9-4.4) | 403 / 546 | 25 |

On tho40 the same node costs thousands of iterations below it (`../tho40_search/README.md`), which is
what its phases had to beat. Here a phase has about 450 iterations to fathom in, at 1.0-1.07x the cost
per iteration (`python3 costs.py ARM`).

**2. The phase climbs about 3 units in its first 100 iterations and 10 in 600** [m, the phase
histories in `runs/z030/dnn.txt`]. So it fathoms inside 450 iterations only from about 8-10 units short
(0.15% of the cutoff), where tho40's fathomed from 1%:

| first-solve shortfall | `z003` phases: fathomed, median iterations | `z030` (pooled starts): fathomed, median iterations |
|---|---|---|
| 0-5 | 20 of 20, 83 | 17 of 17, 160 |
| 5-10 | 25 of 25, 222 | 9 of 9, 510 |
| 10-20 | 24 of 24, 438 | 9 of 15, 1,900 (misses: 1,420, 5 short) |
| 20-40 | | 0 of 18 (790 iterations, 18 short) |
| 40-61 | | 0 of 5 (600, 35 short) |
| 61-184 | | 0 of 10 (650-680, 53-73 short) |

**3. Phase by phase against the baseline's cost under the same node** (`python3 matched.py ARM`; a node
is matched by its set of fixes; phase iterations weighted by their measured cost):

| first-solve shortfall | `z010`: phases (fathomed), arm / base | `z003`: phases (fathomed), arm / base | `g2`: phases (fathomed), arm / base |
|---|---|---|---|
| 0-5 | 14 (14), 0.30 | 20 (20), 0.27 | 10 (10), 0.30 |
| 5-10 | 15 (15), 0.82 | 25 (25), 0.64 | 14 (5), 0.57 |
| 10-15 | 12 (10), 1.00 | 17 (17), 0.92 | 4 (0), 0.52 |
| 15-20 | 5 (3), 1.10 | 7 (7), 1.53 | 10 (0), 0.66 |
| 20-30 | 14 (11), 1.87 | | 9 (0), 0.71 |
| 30-40 | 10 (5), 2.49 | | 9 (0), 0.57 |
| 40-61 | 17 (0), 1.07 | | 5 (0), 0.67 |
| 61-184 (depth 1 and its far children) | | | 10 (0), 0.75 |
| all | 87 (58), 1.25 | 69 (69), 0.73 | 71 (15), 0.69 |

- `z010`: every phase from 20 to 40 short fathoms or nearly does, and costs 2-2.5x what branching would
  have. That is the +39%.
- `z003` keeps the phases that pay and nothing else: 1.10x fewer iterations, and no pool ever forms
  (all 69 phases fathom). `python3 replay.py` on `z010`'s unpooled phases puts the best plain zone at
  0.002-0.003 and its gain at about 3% of the iterations [e]; measured, 0.003 gave 9%.
- `g2`: a stopped phase's 100 iterations buy pooled children that cost 0.5-0.7 of the baseline's.

**4. Which zone wins, then.** Of the brief's four, 0.003 (-7% CPU). 0.01 is the worst place to be: wide
enough to run long phases on depth-2 and depth-3 nodes 20-60 short, too narrow (61.2 units) to reach
any depth-1 node (the nearest is 61.6 short), so no pool comes from above. 0.03 reaches depth 1 and gets the pools,
then spends 35,113 iterations on 39 depth-2 phases (18 fathom) where `g2` spends 4,627:

| | `base` | `z010` | `z030` | `d1` | `g1` | `g2` |
|---|---|---|---|---|---|---|
| depth-1 phase iterations | 0 | 0 | 6,570 | 6,570 | 900 | 6,570 |
| depth 2: fathomed in the first solve, of 225 | 164 | 164 | 186 | 186 | 169 | 186 |
| depth 2: first-solve iterations | 23,770 | 23,770 | 17,261 | 17,261 | 21,621 | 17,261 |
| depth 2: phases (fathomed), iterations | | 58 (40), 51,568 | 39 (18), 35,113 | 0 | 56 (14), 6,068 | 39 (8), 4,627 |
| nodes at depth 3 and below | 3,453 | 836 | 794 | 2,225 | 1,845 | 1,411 |
| iterations at depth 3 and below (first + phases) | 85,356 | 57,076 | 47,470 | 57,920 | 56,926 | 47,277 |
| total iterations | 110,926 | 134,214 | 108,214 | 83,551 | 87,315 | 77,535 |
| CPU s | 3,030 | 4,224 | 3,592 | 2,434 | 2,513 | 2,309 |

- `d1` alone (nine phases, 8% of its iterations) is 1.33x fewer iterations and 20% less CPU: most of the
  gain. The pooled first solves below it take 19-27% fewer iterations than the same nodes in the baseline
  at every depth (`python3 costs.py d1`) and fathom a few more.
- `g1` shows the depth-1 phases need their length: stopped at 100 iterations their pools fathom 169 of
  225, not 186.
- `z030` has 1.03x fewer iterations than the baseline and 19% more CPU: its iterations sit higher in the
  tree, where each costs more (44 ms in a depth-1 phase, 33-36 at depth 2, 19-24 at depths 4-5).
- `sko`: sko42's `dnn.*` block does not change the picture. Against `z010` it has 18% fewer iterations
  (76 phases, 58,515 iterations, where `z010` has 87 and 81,467; one depth-1 node, 0:0, falls inside its zone
  after 300 first-solve iterations) and each costs 0.84 of `z010`'s (`dnn.amx = 1`): 2,899 CPU s, which is
  16% MORE than the AMX baseline `basea` (2,505).
  @@G2S@@

**5. Refinements of `g2`** [m]:
- `g5`, `dnn.gu_zone = 0`: **72,147 iterations, 2,197 CPU s, 564 s wall**, the best arm. The record
  recipe's 100-iteration confirmation of a near give-up is what the phase now does better: the first
  solves that precede a phase below depth 1 take 17,930 iterations instead of 23,390, the phases and
  the tree are almost the same (1,674 nodes against 1,646; 185 of 225 fathom at depth 2 against 186).
- `g4`, `cut.stall = 5` (the depth-1 phases run 1,090-2,280 iterations, 11,680 in all, instead of 6,570):
  190 of 225 fathom at depth 2 instead of 186, 75,437 iterations, but 2,420 CPU s and 636 s wall: the
  extra depth-1 iterations (44 ms each) cost more than they save. The stall rule's length is about right.

## Sanity

- `./sanity.sh`: every run of incumbent 6124 has no `New objective value` line, ends with
  `best objective 6124` and a final state with no open node. Since the baseline proves that no
  solution below 6124 exists, no cut-arm fathom contradicts it.
- `python3 pathcheck.py path_g2 path_g2 g2 z003 z010 z030 d1 g1 sko` (output in `pathcheck.out`):
  `path_g2` is `g2` from incumbent 6126 (cutoff 6124.05), so the optimum is not given and nodes holding
  it must not fathom. It reports one `New objective value: 6124`; the permutation in its final state
  costs 6124 on `qaplib/nug30.dat`; the 28 nodes whose fixes agree with it (depths 0-27) all branch
  (the last is enumerated). On the path the pooled values are: depth 1 first solve 6000.23, phase
  6041.98; depth 2 6065.60, 6067.95; depth 3 6096.23, 6099.86; depth 4 6111.83, 6112.77; depth 5 6122.70,
  phase **6124.0000**; depths 6-26 6124.0000 in the first solve and the phase. In the incumbent-6124 arms
  no value on that path passes 6124 either. PASS.
- `path_g2` took 2,065 nodes, 111,153 iterations, 2,866 CPU s: with the optimum unknown and an incumbent
  2 above it, `g2` still costs less than the baseline costs with the optimum given.

## The code change

`cut.giveup_depth = D` (`latest/settings.C`, `latest/bnb.C`; `RUNNING.md`'s settings table has the row):
`cut.giveup` applies only at nodes of depth D or more. Default 0, so nothing changes unless it is set.
Without it the settings cannot say "long phases at depth 1, short ones below": `cut.giveup` stops by the
gap alone, and the depth-1 gaps (44-99 at iteration 100) overlap the far depth-2 ones. `g1` and `d1` are
the two arms that need no new key (-17% and -20% CPU). The `g2`-family runs used a build of this branch
in a worktree (`/home/user/wt-n30`), the earlier arms the build before the change; the solver is
otherwise the same source (commit in each `cmd.txt`).

## Settings

The winning arm, `par/g5.par`: the record recipe with `dnn.gu_zone = 0` (the recipe has 0.003) plus

```
pool.on = 1
cut.grid = 5:6
cut.side = 1
cut.cp = 100:1000
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 1000
cut.zone = 0.03
cut.giveup = 100:5,300:2
cut.giveup_depth = 2
```

`qaplib/nug30.dat`'s second matrix is the 5 x 6 grid (checked: rectilinear, row-major, 6 columns; 3% zeros
against the first matrix's 35%), so `cut.side = 1`. `cut.zone = 0.03` covers every give-up in the tree
(the farthest is 121.8 short). `par/g2.par` is the same with `dnn.gu_zone = 0.003`. Without the new key: `par/d1.par` (`cut.zone = 0.03`, `cut.maxdepth = 1`).

## What this does and does not show

- One instance, and the stop rule (100:5, 300:2) was read off this instance's own phase histories and
  run once. nug30 is solved; what transfers is the mechanism: **where children fathom in about 20
  iterations, phases should be long only high in the tree, where one pool serves hundreds of nodes,
  and short everywhere else.** tho40's finding (a wide zone, pools from high phases) is the same
  mechanism with the opposite price for the deep phases.
- For tho40 and sko42 the measured thing worth carrying over is the split of the stop rule by depth,
  and the question it raises: tho40's `w18` spent 66% of its iterations in phases; how many of those
  would a `cut.giveup` below some depth have saved?
- Not tried: a depth-1 phase longer than the stall rule allows other than `g4`; `cut.giveup_depth = 3`
  (in `z030` the depth-2 phases cost 35,113 iterations, so only with their own, longer, stop rule);
  rule 9's line from the phase's iterate (`cut.xbar = 1`); the M2 Pro, where the record was set
  (Accelerate, a different cost per iteration by depth).

## Files

| | |
|---|---|
| `par/*.par` | the arms |
| `run.sh`, `queue*.sh`, `queue.log` | one whole search per arm (`./run.sh ARM`; `PAR=`, `UB=`, `SOLVER=` override); the order run and each arm's wall and user time |
| `runs/ARM/` | `cmd.txt` (command, commit, times, load), `log.txt` (the solver's report), `nodes.txt` (node trace), `dnn.txt` (DNN trace with each phase's history), `time.txt`; `final.state` is not committed |
| `summarize.py` | the arms table; `--depth ARM` per depth; `--phases` phases by shortfall |
| `matched.py`, `replay.py`, `costs.py` | phases against the baseline's subtrees; zone and cap replay; ms per iteration and the pooled children |
| `pathcheck.py`, `pathcheck.out`, `sanity.sh` | the sanity checks |
| `lscpu.txt` | the machine |
