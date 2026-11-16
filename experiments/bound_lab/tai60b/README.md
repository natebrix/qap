# tai60b and tai80b with the cluster-region cuts (bound_lab agent "tai60b", Nov 16 2026)

Question: do the cluster-region cuts that transformed tai50b (`../regions/README.md`) do the same for
the two larger open instances of the family? Is either within reach of a few thousand cloud
core-hours, and does one root run give a new best lower bound?

**[m]** measured here (the run is named; `runs/TAG/`), **[e]** estimated. Every bound quoted is a
certified value of the solver (a dual certificate at each check). Iterations and nodes are the
metrics; CPU seconds are from a 4-core box that ran 3 to 4 solver processes at once.

## Answer

1. **Both roots give new best lower bounds** (against Kim, Kojima, Toh's plain DNN bounds, the
   published ones given in the brief):

   | instance | best known (checked from `QAPSOLN/*.sln`) | published bound | gap | this root, certified | gap | share of the published gap closed |
   |---|---|---|---|---|---|---|
   | tai60b | 608,215,054 | 592,829,061 | 2.53% | **605,315,696** at 470 + 3,000 iterations [m `rootL`] | **0.48%** | **81%** |
   | tai80b | 818,415,043 | 786,900,529 | 3.85% | **@@R80@@** at about 470 [e] + @@R80IT@@ iterations, still rising [m `root80`] | **@@R80GAP@@** | **@@R80SHARE@@** |

   tai60b's root passes the published bound 250 iterations into the phase (592,871,409 at
   470 + 250 [m `ab_c6`]); tai80b's at phase iteration 250 (787,039,986 [m `root80`]). Neither has
   plateaued: tai60b gains 124k over its last 250 iterations (207k, 167k, 143k, 124k per 250 from
   2,000 on; by that trend about 0.8M more, a plateau near 606.1M [e]); tai80b was stopped by the
   clock while gaining about @@R80RATE@@ per 100.

2. **tai60b's tree is closer to tai30b's than to tai50b's at the top, and then it is a path with a
   fringe**, like tai50b's:
   - **Depths 1 and 2: every off-path child is fathomable at once, except the mirror.** The root
     branches on location 44, then 43: the isolated 2-point cluster. Of 12 random depth-1 children
     10 fathom in the first (pooled, warm) solve in 10-91 iterations, one fathoms after 482 phase
     iterations, and one does not fathom: `44:58`, the child that swaps the two facilities of that
     cluster (its best completion costs 2,786 more than the best known solution; no bound can
     separate that before the bottom of the tree). Under the on-path node 12 of 12 random depth-2
     children fathom in the first solve (10-60 iterations) [m `d1p1`, `d1_44_39`, `d1_44_58`,
     `d2p1_under_path`].
   - **The fringe appears with the first location in a large cluster** (depth 4, location 1 of the
     14-point cluster; depth 3 is a column branch on a facility). There 5 of 5 random children
     fathom in the first solve, and of 9 children whose facility belongs to the same cluster in the
     best known solution 5 fathom in the first solve and 4 do not; 3 of those 4 were given a phase
     and all three end it unfathomed (0.23M, 0.25M and 0.40M short at the 1,000-iteration limit)
     [m `d4p1_under_path`, `d4_1_39`, `d4_1_59`, `d4_1_29`].
   - **A survivor's children mostly die at once, but not all, and the ones that do not are again
     same-cluster ones**: under the depth-4 survivor `1:39`, 8 of 10 random children fathom in the
     first solve (27-90 iterations); the other two (0.24M and 0.62M short) carry a facility of their
     location's cluster, and the one given a phase ends it 0.17M short at the limit
     [m `d5p1_under_1_39`, `d5_1_39_22_48`]. tai50b's survivor had 1 of 10. So the fringe does
     compound for at least two levels.
   - **Which children survive is predictable from the near-optimal solutions**: a child survives
     when the best completion a 2-swap descent finds for it costs less than best known + the
     slack of the bound at that depth (the on-path node's shortfall). That rule gets 26 of the 27
     depth-1 and depth-4 first solves above right (`nearopt.py`, section 4). The best known
     solution has 8 pair swaps within 10,000 of it, 40 within 100,000 and 213 within 1M, all but
     two inside one cluster [m `nearopt.py swaps`]: that is the fringe, and it compounds only as
     far as the bound's slack exceeds those excesses.
   - **The on-path nodes close faster than tai50b's**: 2.72M short at depth 1, 1.58M, 1.07M, 0.73M,
     0.61M, 0.46M, 0.41M, 0.37M, 0.34M, 0.22M, 0.14M, 0.13M at depth 12, then 0.10M, 0.10M, 0.09M,
     0.08M, 0.08M, 0.06M at depths 13-18 (0.45% of the best known value at depth 1, 0.02% at depth
     12, 0.01% at depth 18; tai50b: 1.05%, 0.85%, 0.67% at depths 1-3) [m chain `d*_path`]. That is
     far closer than tai50b and far less close than tai30b (0.01% at depth 3, exact from depth 15).

3. **Estimate for tai60b: 500-2,500 core-hours on this box (5e6-4e7 ADMM iterations, 6e4-4e5
   nodes; about 700-3,500 c7a core-hours), with a range of 400 to more than 10,000** [e]. It is a
   model, not a dive: the tree of all prefixes whose best completion lies within the measured
   slack of the bound, costed with the measured iterations per node class (section 5, with what it
   rests on). About 440 of those core-hours are depths 0-16, where the slack is measured; the rest
   is the extrapolation. So: **tai60b looks within reach of a few thousand core-hours, but not
   safely: the number to get next is the bound's slack at depths 17-30 and off the path**, which
   decides between 400 and 10,000. Knuth dives from the saved root would give it. tai80b's tree was
   not sampled (its root alone costs 8 s per iteration here).

## 1. Build, regions, cost per iteration

**Build.** A clean worktree of this branch's base (26894b2) with `../regions/regions.patch`
(`git apply`, no conflict, no fix-up to the patch) and `make -C latest -j4 && make -C latest dnntest
dnncuts`. Two things about the container, not the patch: without the `scipy-openblas32` wheel the
makefile links the reference BLAS (7x slower, RUNNING.md), so `pip install scipy-openblas32` and a
rebuild came first; `gen_regions.py -stats` needs `scipy`. One diagnostic addition of mine,
`live.patch` (10 lines in `latest/dnncut.C`): `QAP_CP_LIVE=FILE` appends the phase's certified
running best at every check, so that a long root leaves a record if it is killed (`runs/*/live.txt`).
It does not touch the numbers.

**Best known values** (`structure.py`): `sum A[i,k] B[p(i),p(k)]`, A the first matrix, gives
608,215,054 for `tai60b.sln` and 818,415,043 for `tai80b.sln`, the listed values [m].

**Structure** (`gen_regions.py -stats`, `structure.py`) [m]:

| instance | 2-D MDS stress | share of the first axis | natural clusters (sizes) | merge height jump | mean distance within / between | share of sum(D) between |
|---|---|---|---|---|---|---|
| tai50b (for comparison) | 0.0003 | 0.971 | 5 (3, 4, 16, 15, 12) | 100 -> 457 | 71 / 1070 | 0.979 |
| tai60b | 0.0003 | 0.918 | 4 (14, 2, 8, 36); at the next levels the 36 split into 20, 12, 4 | 341 -> 791 | 195 / 1188 | 0.892 |
| tai80b | 0.0005 | 0.781 | 7 (10, 10, 21, 10, 6, 20, 3) | 137 -> 248 | 51 / 640 | 0.983 |

tai60b is the collinear kind, less cleanly than tai50b: one 36-point blob (diameter 546 against 1,188
between clusters) holds 11% of the distance sum inside the "natural" clusters, and the 2-point cluster
{43, 44} lies off the axis. tai80b is two-dimensional (as tai35b and tai40b, where the same recipe
closed 66-71%), with 7 tight clusters carrying 98% of the distance.

**Regions chosen** (location side, `cut.side = 0`; files in `regions/`), on the tai50b pattern
`line:1,top:K,linegap:4:0.03`:
- tai60b: **K = 6** (`tai60b.comb6.txt`, 86 regions: 59 prefixes along the principal axis, the 31
  unions of the six clusters 14, 2, 8, 20, 12, 4, and the between-cluster half-plane cuts in four
  directions). K = 6 and not the "natural" 4 because the 36-point blob is not tight: its three parts
  are as far apart as tai50b's clusters are wide. Checked against K = 4 (`tai60b.comb4.txt`, 64
  regions) at the root [m `rootM`, `c4M`, same settings]: equal to 250 phase iterations (592.87M
  against 592.94M), then K = 6 ahead: 601.51M against 601.37M at 700, **602,989,748 against
  602,786,441 at 1,000**, at 0.4% more CPU per iteration.
- tai80b: **K = 7**, the natural clusters (`tai80b.comb7.txt`, 138 regions). Not compared with
  anything: no time at 8 s per iteration. With clusters spread in two dimensions `line:2` or `line:4`
  may do better (not tried here or on tai35b/40b).

**Phase penalty**: `dnn.beta_sched = w=0.36,wc=2.88` from the start. At the tai60b root [m `ab_c6`,
`ab_c6_nowc`, 470 + 250 iterations]: 592,871,409 with `wc`, 587,781,768 without; the run without
reaches at 250 what the run with reaches at about 115.

**Cost** [m]:

| | m = 50 (tai50b, `../regions/size50`) | m = 60 | m = 80 |
|---|---|---|---|
| s per iteration, first solve (cold root) | 0.53 | 1.36 (`rootL`: 470 iterations, 647 CPU-s) | about 7 (`root80`: 3,316 CPU-s for about 470 [e]) |
| s per iteration, phase | 0.72 loaded, 0.4-0.5 idle | 1.53 (`rootL`: 3,000 in 4,649 CPU-s), 1.1-1.3 at depth 1 | 7.9-8.3 (`root80`) |
| peak memory | | 2.5 GB | 7.3 GB in the phase |
| root state file (`ship.bits = 64`) | 134 MB | 277 MB | not written (run stopped) |

So 3x and 16x tai50b, as expected from (m-1)^6. Cores: every solver process had one core
(`OPENBLAS_NUM_THREADS=1`; a tai80b run started with 2 BLAS threads still used one). The container
has 4 cores and a 14.3 GB memory limit: **the first tai80b root was killed for memory** after 51
minutes when three tai60b processes ran beside it (`runs/root80_oomkilled`); after that at most two
tai60b processes ran beside tai80b, and one core idled. `dnn.amx` was not tried (`par/ab_c6_amx.par`
prepared; the CPU has AMX-BF16).

## 2. tai60b root [m]

`root.sh NAME`: `bnbdrv -p par/NAME.par -u 608215054 -d 0 -o runs/NAME/root.state` (the settings of
`../regions/size50/par/rootL.par` with `wc=2.88`, `cut.stall = 0`, `cut.cp_iters` as below). The
first solve stops as "hopeless" at 470 iterations at **581,315,106** (26.9M short): that is the plain
bound at the give-up, not converged; Kim, Kojima and Toh's 592,829,061 is the converged plain DNN
bound, which I did not rerun.

| phase iterations | certified (`rootL`, comb6, wc=2.88) | short of best known | share of the published gap closed |
|---|---|---|---|
| 0 (470 plain) | 581,315,106 | 26,899,948 | (below the published bound) |
| 250 | 592,871,409 | 15,343,645 | 0.3% |
| 500 | 599,190,046 | 9,025,008 | 41% |
| 1,000 | 602,989,748 (= `rootM`) | 5,225,306 | 66% |
| 1,500 | 604,106,475 | 4,108,579 | 73% |
| 2,000 | 604,674,349 | 3,540,705 | 77% |
| 2,500 | 605,048,662 | 3,166,392 | 79% |
| 3,000 | **605,315,696** | **2,899,358** | **81%** |

5,296 CPU-s in all. Not a plateau (gains per 250: 207k, 167k, 143k, 124k); the on-path depth-1 child
reaches 605.49M and the depth-2 one 606.64M, so a longer root has at most about 1M to give.
`rootM` (1,000 phase iterations, 5.23M short) is the root every tree sample below starts from: it was
ready 50 minutes earlier. A dive runner should start from the longer one.

## 3. tai60b tree samples [m]

Protocol of `../regions/size50`: `kid.sh` (one node warm from its parent's state with the pool, first
solve then phase, `par/kidW.par`: `cut.cp_iters = 1000`, `cut.stall = 3`), `pass1.sh` + `cutset.py`
(first solves only of a random subset, `par/kid0.par`), `chain.sh` (the on-path node of each level,
incumbent at best known + 1, cutoff 608,215,054.05). `F:L` as in the node trace: entry F (a location,
the first matrix) of the permutation holds facility L; the on-path child of a branch on row r is
`r:p[r]` with p from `tai60b.sln` (`onpath.txt`). `collect.py` prints the table of all runs.

**The path** (`runs/d*path`; first solves are 110 iterations each):

| depth | child (m) | first solve: certified, short | phase end: certified | short | % of best known | phase iterations (end) | CPU-s | s/it |
|---|---|---|---|---|---|---|---|---|
| 1 | `44:45` (59) | 603,840,108, 4,374,946 | 605,492,773 | 2,722,281 | 0.448% | 780 (stalled) | 1,185 | 1.33 |
| 2 | `43:58` (58) | 606,006,543, 2,208,511 | 606,638,236 | 1,576,818 | 0.259% | 770 (stalled) | 968 | 1.10 |
| 3 | `45:26` (57) | 606,838,380, 1,376,674 | 607,149,697 | 1,065,357 | 0.175% | 600 (stalled) | 678 | 0.95 |
| 4 | `1:33` (56) | 607,144,783, 1,070,271 | 607,481,185 | 733,869 | 0.121% | 830 (stalled) | 803 | 0.85 |
| 5 | `6:24` (55) | 607,330,940, 884,114 | 607,605,946 | 609,108 | 0.100% | 820 (stalled) | 707 | 0.76 |
| 6 | `7:19` (54) | 607,184,172, 1,030,882 | 607,752,706 | 462,348 | 0.076% | 1000 (limit) | 805 | 0.72 |
| 7 | `0:32` (53) | 607,757,734, 457,320 | 607,805,120 | 409,934 | 0.067% | 600 (stalled) | 503 | 0.71 |
| 8 | `14:54` (52) | 607,767,805, 447,249 | 607,841,134 | 373,920 | 0.061% | 600 (stalled) | 463 | 0.65 |
| 9 | `24:2` (51) | 607,818,783, 396,271 | 607,878,851 | 336,203 | 0.055% | 600 (stalled) | 401 | 0.56 |
| 10 | `28:50` (50) | 607,774,646, 440,408 | 607,995,616 | 219,438 | 0.036% | 1000 (limit) | 603 | 0.54 |
| 11 | `18:43` (49) | 607,914,004, 301,050 | 608,071,159 | 143,895 | 0.024% | 1000 (limit) | 543 | 0.49 |
| 12 | `23:18` (48) | 607,998,639, 216,415 | 608,087,137 | 127,917 | 0.021% | 1000 (limit) | 515 | 0.46 |
| 13 | `38:41` (47) | 608,033,437, 181,617 | 608,111,140 | 103,914 | 0.017% | 1000 (limit) | 468 | 0.42 |
| 14 | `15:12` (46) | 607,879,688, 335,366 | 608,113,739 | 101,315 | 0.017% | 1000 (limit) | 478 | 0.42 |
| 15 | `31:10` (45) | 608,097,002, 118,052 | 608,121,522 | 93,532 | 0.015% | 600 (stalled) | 314 | 0.44 |
| 16 | `20:21` (44) | 608,084,974, 130,080 | 608,130,965 | 84,089 | 0.014% | 1000 (limit) | 426 | 0.38 |
| 17 | `16:17` (43) | 608,086,740, 128,314 | 608,136,872 | 78,182 | 0.013% | 900 (stalled) | 357 | 0.35 |
| 18 | `19:13` (42) | 608,113,127, 101,927 | 608,154,280 | 60,774 | 0.010% | 920 (stalled) | 383 | 0.37 |

Shortfalls are against the best known value 608,215,054 (the cutoff is that + 0.05). "limit": the phase
reached `cut.cp_iters = 1000` still gaining; "stalled": stopped by `cut.stall = 3`.

The branching order is the search's own (rule 9): the 2-point cluster's locations 44 and 43, a column
branch on facility 26 (its children are the 58 free locations; on-path 45), then locations of the
14-point cluster (1, 6, 7, 0), then others.

**Depth 1** (`d1p1`: 12 random children + the on-path one, first solves; then phases):

| child | first solve | phase | outcome |
|---|---|---|---|
| `44:2`, `44:5`, `44:9`, `44:22`, `44:23`, `44:25`, `44:40`, `44:42`, `44:54`, `44:59` | fathomed in 10-91 iterations (13-100 CPU-s) | | 10 of 12 fathom at once |
| `44:39` | 1,958,421 short at 110 | fathomed after 482 iterations | fathomed (592 iterations, 741 CPU-s) |
| `44:58` (the swap of the cluster's two facilities; excess 2,786) | 3,846,034 short | stalled 2,293,578 short after 820 | **branches**: the mirror of the path |
| `44:45` (on-path) | 4,374,946 short | stalled 2,722,281 short after 780 | branches |

**Depth 2** (`d2p1_under_path`: 12 random children `43:x` of `44:45`): **12 of 12 fathom in the first
solve**, 10-60 iterations, 13-60 CPU-s. On-path `43:58`: 2.21M short, phase stalls 1.58M short.

**Depth 3**: a column branch (facility 26 over 58 locations); only the on-path child was run.

**Depth 4** (`d4p1_under_path`: children `1:x` of `44:45,43:58,45:26`; location 1 is in the 14-point
cluster, whose facilities in the best known solution are 0, 19, 23, 24, 29, 32, 33, 35, 36, 39, 49,
55, 57, 59):

| child | kind | excess of its best completion (`nearopt.py`) | first solve | phase |
|---|---|---|---|---|
| `1:8`, `1:13`, `1:17`, `1:50`, `1:53` | random (5 of 5), other clusters | 13.6M-34.3M (`1:53`, `1:8`, `1:13`, `1:17`, `1:50`: 13.6M, 15.9M, 30.3M, 33.0M, 42.9M) | fathomed, 23-50 iterations | |
| `1:23`, `1:0` | same cluster | 2.86M, 2.47M | fathomed, 56 and 70 iterations | |
| `1:32`, `1:35`, `1:36` | same cluster | 1.57M, 1.33M, 1.23M | fathomed, 87-117 iterations (by 2k-14k) | |
| `1:29` | same cluster | 0.96M | 705,728 short | 247,676 short at the 1,000 limit, gaining 25k per 100: **branches** |
| `1:24` | same cluster | 0.51M | 1,124,476 short | not run |
| `1:59` | same cluster | 0.51M | 798,843 short | 231,441 short at the 1,000 limit: **branches** |
| `1:39` | same cluster | 0.25M | 922,686 short | 396,970 short at the 1,000 limit: **branches** |
| `1:33` (on-path) | | 0 | 1,070,271 short | stalled 733,869 short after 830 |

So at depth 4 the 43 other-cluster children fathom at once [e from 5 of 5], and of the 13 same-cluster
ones about 8 fathom in the first solve and **4-5 survive** (the on-path one included) [e from 9 first
solves and 3 phases, all three unfathomed at the limit]. tai50b had about 15 of 46 surviving at its
depth 4. The phases of these survivors end 0.74M to 1.2M below the upper bounds on their subtrees'
optima (shortfall + excess), against 0.73M for the on-path node: off the path the slack is up to 1.6
times larger.

**A survivor's children** (`d5p1_under_1_39`: first solves of 10 random children of `1:39`, which
branches on location 22, in the 12-point cluster {14..25}, whose facilities in the best known solution
are 2, 9, 12, 13, 17, 18, 21, 25, 31, 43, 48, 54): 8 fathom in the first solve (27-90 iterations,
20-62 CPU-s), among them one same-cluster facility (`22:17`); `22:9` ends 237,006 short and `22:48`
624,649 short, both same-cluster. Phases: `22:48` 174,189 short at the 1,000 limit: branches
(`d5_1_39_22_48`); `22:9` @@D5B@@.

## 4. Why the fringe is what it is: near-optimal solutions (`nearopt.py`) [m for the counts]

- Pair swaps of the best known permutation: 3 cost less than 1,000 more, 8 less than 10,000, 18 less
  than 50,000, 40 less than 100,000, 77 less than 250,000, 123 less than 500,000, 213 less than 1M;
  211 of those 213 exchange two locations of the same cluster. The 8-point cluster {52..59} (diameter
  14) alone has 22 of the 40 cheapest: almost any arrangement of its eight facilities is within 0.1M.
- A node whose subtree holds a solution with excess e cannot fathom before the bound's slack at its
  depth falls below e. Measured at depth 1 and 4: children with excess up to 0.96M did not fathom in
  the first solve (first-solve slack of the on-path node there: 4.37M and 1.07M), children from 1.23M
  up did at depth 4; at depth 1 the first solve fathomed children with excess 6.28M and more (the one
  at 6.28M, `44:40`, by 22k after 91 iterations) and missed one at 6.75M (`44:39`, 1.96M short, then
  fathomed by its phase). Of the 27 sampled first solves at depths 1 and 4, the rule "fathoms in the
  first solve iff excess > the on-path node's first-solve shortfall" is wrong on that one only: the
  slack off the path is somewhat larger than on it. The same shows in the phases: the three depth-4
  survivors end 0.74M-1.2M below their subtrees' upper bounds, the on-path node 0.73M below its optimum.
- This is the tai50b agent's "same-cluster children survive" seen from the solutions' side, and it
  says when the fringe stops: when the slack is below the swap excesses, i.e. gradually from about
  0.5M (123 swaps) down to 10k (8 swaps).

## 5. Tree-size estimate for tai60b (`model.py`) [e]

**Model.** Fix a branching order of the rows. A depth-d prefix *branches* when the best completion
found for it (2-swap best-improvement descent from its parent's completion) has excess below s(d),
the phase's slack at depth d; it *fathoms in a phase* when its excess is between s(d) and 1.6 s(d)
(the first-solve slack: on-path first-solve shortfalls are 1.3-2.2 times the phase-end ones); otherwise
it *fathoms in the first solve*. s(d) = the on-path node's shortfall at the end of its phase, measured
for d = 1..16 (or 18), continued below by a ratio r per level. Costs from section 3: a branching node
110 + 800 iterations, a phase fathom 110 + 400, a first-solve fathom 45; 1.3 s per iteration at
m = 59 scaled by (m/59)^6 (measured: 0.93-0.95 s at m = 56-57, 0.72 at m = 54, 0.46 at m = 48).
Order: the search's measured rows down the path, then the remaining locations cluster by cluster
(`path`) or at random (`rand`).

s(d) measured to depth 16 (2.72M ... 0.084M, the table in section 3); `logs/model16_*.txt` have the
per-depth profiles:

| order | s(d) for d = 1..16 | ratio r below depth 16 | branching nodes | nodes bounded | ADMM iterations | core-hours on this box |
|---|---|---|---|---|---|---|
| path | measured | 0.8 | 987 | 43,057 | 3.4e6 | 475 |
| path | measured | 0.9 | 4,090 | 137,417 | 1.1e7 | 661 |
| random (seed 1) | measured | 0.9 | 1,805 | 60,347 | 4.9e6 | 492 |
| path | measured | 0.95 | 22,888 | 635,972 | 5.3e7 | 1,378 |
| path | measured x 1.5 | 0.9 | 11,574 | 407,171 | 3.4e7 | 2,362 |
| path | measured x 1.5 | 0.95 | 76,303 | 2.2e6 | 1.8e8 | 4,925 |
| path | measured | 1.0 (the slack never falls below 0.084M) | 1.7e6 | 1.8e7 | 2.3e9 | 7,918 |
| path, with depths 17-18 measured too (0.078M, 0.061M; rows 16, 19) | measured | 0.9 below depth 18 | 2,482 | 83,063 | 6.7e6 | 524 |
| same | measured | 0.95 | 10,458 | 271,358 | 2.3e7 | 723 |

**Reading.**
- **Depths 0-16, where s(d) is measured, cost about 440 core-hours in the model** (560 branching
  nodes and 2.2M iterations; most of it at depths 4-10, where a node costs 0.5-1 s per iteration and
  the same-cluster children of the large clusters survive). That part does not depend on r, only on
  the survive rule; with the slack x 1.5 it is about 1,700 (2,060 branching nodes, 9.0M iterations).
- **Below depth 16 everything depends on how fast the slack keeps falling.** Measured ratios per
  level: 0.58-0.83 at depths 2-6, 0.89-0.91 at depths 7-9, 0.65 at depths 10-11, then 0.89, 0.81,
  0.97, 0.92, 0.90 at depths 12-16: about 0.9 lately. At r = 0.9 the lower tree adds 220 core-hours;
  at 0.95, 940; if the slack stopped at 0.084M (40 pair swaps are cheaper than that) the fringe would
  run to the bottom: 7,500, and with the larger slack as well, several times that.
- The last two rows came in at the end (depths 17 and 18: ratios 0.93 and 0.78, in line with 0.9) and
  are lower than the rows above them at the same r, because two more rows of the order are the
  search's own: the model moves by a factor of 2 with the order alone.
- **Central reading: 500-2,500 core-hours on this box (5e6-4e7 iterations, 6e4-4e5 nodes), about
  700-3,500 c7a core-hours; range 400 to more than 10,000** [e]. The x 1.5 rows are not a
  pessimistic decoration: the three depth-4 survivors I ran sit 1.0-1.6 times the on-path slack below
  their subtrees' upper bounds.
- For comparison, the tai50b agent's cross-check for tai50b was 50-200 core-hours if its fringe dies
  in one or two levels and 1e4-1e5 if it compounds; here the fringe is seen compounding for two
  levels (depth 4 and 5) and the model says it goes on at a slowly falling rate.

**What it rests on.** One root; the on-path chain to depth 18 (one node per level); 13 + 12 + 14 + 10
first solves at depths 1, 2, 4, 5 and 7 off-path phases; the survive rule checked on those nodes
only at depths 1 and 4; a descent that finds upper bounds on subtree optima (so it can only
under-count survivors where it is poor); one fixed branching order where the search picks its own in
every subtree (the survivor `1:39` branched on location 22, the path on location 6); phases stopped by
`cut.stall = 3` or the 1,000 limit, so s(d) is what these settings leave, not the relaxation's limit
(the root: 5.2M short at 1,000 iterations, 2.9M at 3,000). The cost side assumes off-path survivors
cost what on-path nodes do.

**What would move it.** (a) The slack below depth 18 and off the path: r = 0.9 against 1.0 is a factor
of 12, the slack x 1.5 a factor of 3.5. On tai30b the on-path bound was within 0.01% from depth 3 and
within 1 from depth 15; here it is at 0.01% at depth 18. (b) Longer phases at survivors and a longer root lower
s(d) at a cost per node the model does not trade off. (c) Symmetric-looking arrangements of the
8-point cluster: if the bound cannot tell them apart until they are fixed, that cluster alone is up
to 8! = 40,320 leaves; the model has this through the swap excesses, but only as far as the descent
finds them. (d) c7a against this box: `../tho40_search` measured this box at 0.70 of c7a per matched
node, so core-hours here are about 1.4 c7a core-hours when the box is loaded [e].

## 6. Validity [m]

With the incumbent at best known + 1 (cutoff 608,215,054.05) no node on the path to
`QAPSOLN/tai60b.sln` fathoms and no certified value there exceeds 608,215,054: 18 nodes, depths 1-18,
the largest 608,154,280 at depth 18 (table in section 3). Off the path, the mirror `44:58` (subtree optimum at most 608,217,840)
ends at 605,921,475 and the depth-4 same-cluster survivors at 607,818,084 and 607,983,613, below
their subtrees' upper bounds 608,463k and 608,727k. The root values are below the best known value.
No brute-force check was rerun here; `../regions/README.md` section 3 has those for the patch.

## 7. tai80b root [m `root80`]

`root.sh root80 tai80b 818415043` (`par/root80.par`: the same settings, `tai80b.comb7.txt`, `dnn.depths
= 0:78`, `cut.cp_iters = 5000`). One core; two tai60b processes beside it. The run was stopped by
hand at the end of the session, so its values are the certified running best at each check
(`runs/root80/live.txt`), and the first solve's iteration count was never printed (about 470 [e]
from 3,316 CPU-s at about 7 s).

@@R80TABLE@@

The first solve gave up at 776,697,155 (41.7M short; the published converged plain bound is
786,900,529). @@R80READ@@

## 8. Not measured

- tai60b: a search or a Knuth dive; anything off the path below depth 5; a phase for `1:24`; the survivors' phases beyond 1,000 iterations (all four were still gaining); any
  depth-3 (column-branch) sibling; the root beyond 3,000 phase iterations; the converged plain bound;
  samples warm from the 3,000-iteration root (all samples start from the 1,000-iteration one);
  `dnn.amx`; three Dykstra sweeps / extrapolation (`../cpfast`).
- tai80b: anything but the root; other region sets; the root to a plateau.
- Whether 608,215,054 and 818,415,043 are optimal.

## Files

- `root.sh`, `kid.sh`, `pass1.sh`, `cutset.py`, `chain.sh`, `collect.py`, `mkpar.sh`, `rss.sh`, the
  queues `q*.sh`; `structure.py`, `nearopt.py`, `model.py`, `pathrows.txt`, `onpath.txt`; `live.patch`.
- `par/`: `rootL` (3,000), `rootM` (1,000), `c4M` (comb4), `ab_*` (250-iteration A/B), `kidW`, `kid0`,
  `root80`. They carry absolute paths to `regions/`.
- `regions/`: `tai60b.comb6.txt` (used), `tai60b.comb4.txt`, `tai80b.comb7.txt`.
- `runs/TAG/`: `cmd.txt`, `log.txt`, `nodes.txt`, `dnn.txt` (with the phase history), `live.txt`.
  State files are not committed (277 MB each).
- `logs/`: timing runs (`time60.log`, `time80.log`), queue outputs, `model*.txt`.
