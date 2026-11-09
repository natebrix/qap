# Is tai50b reachable with the certified DNN bound? (Nov 8 2026)

Companion to `../sko49_feasibility/NOTES.md` and `../wil50_feasibility/NOTES.md`,
run on the 2-core cloud Xeon @ 2.1 GHz (AMX; OpenBLAS `scipy-openblas32`, one thread
per process, two processes at a time). Every bound is a certified value printed by
`dnntest`. Scripts and logs: `run.sh`, `root/`, `depth3/`; `../sko49_feasibility/
summarize.py depth3 458821516.05` builds the depth-3 table.

## Summary

- Root DNN bound **430,147,161 certified, 28.7M (6.25%) below the best known
  458,821,517**. QAPLIB's best bound is BBCPOP's 431,090,700 (6.04%); ours is 0.2% under
  it at 1,500 iterations, still rising slowly.
- **2 of 12 random depth-3 nodes fathom** (at 600 and 630 iterations, 332-349 CPU-s); the
  other 10 stop 2.7M-14.4M short (0.58-3.1% of the cutoff, median 6.8M).
- From the root to depth 3 the shortfall falls 5.6x (28.7M to a median of 5.1M over all
  12), about **1.8x per level**. wil50's falls about 1.25x per level; tho40's halves. On
  sko49 none of 12 depth-3 nodes fathomed.
- So the relative root gap overstates how hard the top of the tree is, as it did for tai40b
  (4.5% at the root, solved in 278,465 nodes by Fujii et al. 2021).
- **Tuned (section 5): with the DNN through depth 9 and a 600-iteration give-up, 10 dives
  give 4.1e9 ± 3.0e9 CPU-s, about 47 days on 1,000 cores**: the depth-6-to-8 levels that
  made the first estimate shrink 10-500x, each node costs 8x more, and level 5 now holds
  81% of the CPU. Against the first estimate's 88-day low reading that is about 2x; against
  its 20-dive mean, 45x, most of which is one dive the tuned run did not repeat.
- **The first Knuth estimate (20 dives, section 4) puts the tree at 1.8e11 ± 1.7e11 CPU-s on
  this Xeon (about 5.8 years on 1,000 cores), from one dive carrying 96% of it; without
  that dive, 7.6e9 CPU-s (88 days).** Either way that is far above Fujii et al.'s 5.5e7
  nodes (rule M), and above wil50's roadmap center. The tail comes from depths 7-8, where
  the recipe's DNN gives up at 200 iterations (depth 7) or is not used (depth 8 and
  below): raising the give-up and the DNN depth for tai50b is the first thing to test.

## 0. The instance

- `qaplib/tai50b.dat`: n = 50, Taillard's "real-life-like" type b. The first matrix is
  symmetric (distances, 0-1,843, 98% nonzero); the second is asymmetric (flows, 0-8,071,
  43% nonzero, heavy-tailed: median nonzero flow 58, and the top 10% of nonzero flows carry
  61% of the total). bnbdrv and dnntest symmetrize the second matrix (objective
  values unchanged) and use fathoming tolerance 0.95 (integer, not even), **cutoff
  458,821,516.05**.
- **Best known 458,821,517** (`QAPSOLN/tai50b.sln`, Ro-TS); the permutation evaluates to
  458,821,517 with the first matrix as the distance (757,932,044 the other way round).
  Never proved optimal.
- Objectives are doubles in the solver (`solution.h`), so the 4.6e8 magnitude is no
  overflow risk; the largest possible sum, n² · 1,843 · 8,071 = 3.7e10, also fits.
- No grid and no symmetry file, so neither the symmetry file nor the regional cuts of
  wil50 and sko49 apply. Rule 9 branching, warm starts, pool inheritance and the
  per-node speed work do.

## 1. Root (`root/root.log`)

`dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 458821517 -v`, the sko49 note's
settings: 421,577,777 at 400 iterations, 429,428,065 at 1,000, **430,147,161 at
1,500** (gaining about 1.4M per 1,000 iterations at the end); 833 ms per iteration
with both cores busy, 1,250 CPU-s.

## 2. Depth-3 sample (`depth3/`, 12 random nodes, seed 20261108)

Drawn as the sko49 and wil50 samples were (three random facilities paired with three
random locations, `depth3/draw_nodes.py`), same `dnntest` settings, cutoff
458,821,516.05. "At 400" is the certified bound at the search's give-up point.

| node | spec | lb at 400 | short at 400 | final lb | short | status | iters | CPU s |
|---|---|---|---|---|---|---|---|---|
| 1 | `20:19,28:24,14:35` | 449,539,545 | 9,281,971 | 454,973,581 | 3,847,935 | iter_limit | 1500 | 847 |
| 2 | `0:42,35:27,48:5` | 442,863,825 | 15,957,691 | 449,238,914 | 9,582,602 | iter_limit | 1500 | 902 |
| 3 | `27:0,44:23,6:19` | 449,930,773 | 8,890,743 | 455,582,587 | 3,238,929 | iter_limit | 1500 | 857 |
| 4 | `42:5,25:16,29:23` | 456,140,025 | 2,681,492 | 458,849,521 | fathomed | fathomed | 630 | 349 |
| 5 | `44:27,27:36,0:21` | 441,358,403 | 17,463,113 | 448,190,032 | 10,631,484 | iter_limit | 1500 | 818 |
| 6 | `12:31,27:41,25:11` | 455,943,022 | 2,878,494 | 458,823,193 | fathomed | fathomed | 600 | 332 |
| 7 | `0:44,49:49,25:31` | 451,179,509 | 7,642,007 | 456,142,196 | 2,679,320 | iter_limit | 1500 | 854 |
| 8 | `43:29,6:45,28:31` | 447,419,407 | 11,402,109 | 452,625,028 | 6,196,488 | iter_limit | 1500 | 845 |
| 9 | `4:41,35:32,24:5` | 448,550,160 | 10,271,356 | 454,731,218 | 4,090,298 | iter_limit | 1500 | 821 |
| 10 | `5:40,40:22,44:3` | 440,328,635 | 18,492,881 | 444,384,291 | 14,437,225 | iter_limit | 1500 | 807 |
| 11 | `17:0,31:16,9:33` | 443,884,812 | 14,936,704 | 451,515,658 | 7,305,858 | iter_limit | 1500 | 811 |
| 12 | `25:25,6:40,16:13` | 443,972,299 | 14,849,217 | 451,078,758 | 7,742,758 | iter_limit | 1500 | 872 |

- 538-601 ms per iteration at m = 47 with both cores busy; 810-900 CPU-s per unfathomed
  node at 1,500 iterations.
- The bound gains 4M-8M between iteration 400 and 1,500 at the unfathomed nodes, a larger
  share of the remaining shortfall than on sko49 (10-24 units of 134-337): these nodes are
  further from convergence at the give-up point, so the give-up iteration matters more here.
- Both fathomed nodes were 2.7M-2.9M short at 400 iterations, about the same as the
  closest unfathomed node at the end (node 7). The two classes overlap at the give-up point.

| | tai50b (12 nodes) | sko49 (12 nodes) | wil50 (16 nodes) |
|---|---|---|---|
| root gap, % of cutoff | 6.25% | 3.14% | about 1.4% |
| fathomed at depth 3 | **2 of 12** | 0 of 12 | 2 of 16 |
| unfathomed shortfall, % of cutoff | **0.58-3.1%, median 1.5%** | 0.57-1.44%, median 1.16% | 0.12-0.83% |
| shortfall, root to depth-3 median | **5.6x (about 1.8x per level)** | 2.7x (1.4x) | about 2x (1.25x) |

## 3. What this says, and what is next

- tai50b's bound closes much faster with depth than sko49's or wil50's. That fits its type:
  a minority of large flows carries most of the objective, and fixing the facilities that carry them settles most of
  the objective. Which three facilities a node fixes matters a lot (shortfalls 2.7M to
  14.4M at depth 3), so branching order should matter more here than on the grid instances.
- Per node, a tai50b DNN solve costs what a wil50 one does (same m, 0.55-0.6 s per
  iteration here).
- The tree size is not settled by 12 nodes; section 4 estimates it.

## 4. Knuth estimate (`estimate/`, 20 dives, Nov 8)

`estimate/run_est.sh SEED`: `bnbdrv -e` with the nug30 record recipe
(`../nug30_solve/README.md`: `nug_deep_xbar_all.par`, rule 9 through depth 6, DNN at
depths 1-7 with give-up 200, the Nov 7 environment), the DNN-marginal sampler with a 30%
uniform share (`QAP_EST_XBAR=1 QAP_EST_XBAR_MIX=0.3`), dives from level 1, the per-level
table printed after every dive (`estimate.par`: `1 10000 1 1 1.5 1000`). Seeds 41 and 42,
one per core, two hours each, 10 dives each. Logs in `estimate/logs/`; `dives.csv` (one
row per dive), `dive_levels.csv` (per dive and level) and `levels.csv` (mean and SE per
level) are derived from the 64-bit tables in the logs (dive k's share of a level is k
times the k-dive mean less k-1 times the previous one; level 1's ±0.9 is rounding of the
printed means). The `dives_s*.txt` node column wraps at 2^31 for the largest dives: the
history line printed an `int`, fixed in 3a678aa after these runs started.

| | value (on this Xeon) |
|---|---|
| nodes, mean of 20 dives | **1.8e10 ± 1.7e10** |
| CPU-s, mean of 20 dives | **1.8e11 ± 1.7e11 (about 2,100 days, 5.8 years, on 1,000 cores)** |
| largest dive's share of the CPU sum | 96% (seed 41, dive 6: 3.5e11 nodes, 3.5e12 CPU-s) |
| without that dive (19 dives) | 5.3e8 nodes, 7.6e9 CPU-s (88 days on 1,000 cores) |
| per seed (10 dives each) | seed 41 3.6e11 CPU-s, seed 42 8.0e9 |
| median dive | 9.2e6 nodes, 3.1e8 CPU-s |
| depths where dives ended | 2 (1), 3 (3), 4 (6), 5 (5), 6 (2), 7 (1), 8 (1), 9 (2) |

Per level (mean over the 20 dives ± SE; "reach" is how many dives had nodes there):

| level | nodes | CPU-s | reach |
|---|---|---|---|
| 1 | 50 | 1.2e4 | 20 |
| 2 | 2,806 ± 498 | 2.8e5 ± 5.1e4 | 20 |
| 3 | 1.7e5 ± 5.5e4 | 5.4e7 ± 4.2e7 | 19 |
| 4 | 3.1e6 ± 1.3e6 | 1.6e8 ± 9.3e7 | 16 |
| 5 | 1.4e8 ± 1.2e8 | 4.4e9 ± 2.3e9 | 14 |
| 6 | 4.0e8 ± 3.3e8 | 8.4e9 ± 6.6e9 | 11 |
| 7 | 4.2e8 ± 4.1e8 | 3.9e10 ± 3.8e10 | 7 |
| 8 | 1.7e10 ± 1.7e10 | 1.3e11 ± 1.3e11 | 7 |
| 9 | 2.1e7 ± 1.6e7 | 1.4e8 ± 1.2e8 | 5 |

Reading:

- **The estimate has not converged.** One dive carries 96% of the CPU, the two seeds differ
  45x, and the SE equals the mean. The 88-day figure without the largest dive is a low reading
  (dropping a dive biases the mean down), and the 5.8 years is what the sample says so far;
  a deeper dive than any seen could raise it again.
- **Where the tail is.** In the largest dive, 8.3e9 depth-7 nodes (bounded by the DNN with
  give-up 200) branched into 3.4e11 depth-8 nodes (QPB only), which then fathomed (3 nodes
  at depth 9). Level 8 holds 95% of the estimated nodes and 72% of the CPU over all dives.
  That is the recipe tuned for nug30 (DNN through depth 7, give-up 200) meeting an instance
  whose bound gains a lot after the give-up point: at depth 3 the unfathomed nodes gained
  4M-8M between iteration 400 and 1,500 (section 2).
- **Next, cheapest first:** (1) the same estimate with a later give-up at depths 5-7 and
  the DNN through depth 8 or 9 (`-x 1:9:1:H:1`, H 400-600), to see how much of level 8
  goes away; (2) more dives, ideally started from an exact depth-2 or depth-3 frontier
  (`start_dive_level`), so the top levels add no sampling error; (3) the per-level node
  cost: about 10 CPU-s per node overall here, 0.55-0.6 s per DNN iteration at m = 47.
- Against Fujii et al.'s 5.5e7 nodes (rule M, Newton-bracketing bounds, no QPB): their
  nodes are all DNN-bounded, so the counts don't compare directly; their tai40b took
  1.8e8 core-seconds (278,465 nodes, 1.05e5 s on 1,728 cores).

## 5. Tuned: DNN through depth 9, give-up 600 (`estimate/x9_600/`, 10 dives, Nov 8)

Section 4's run with one change, `-x 1:9:1:600:1` in place of `-x 1:7:1:200:1`
(`X=1:9:1:600:1 TAG=x9_600 estimate/run_est.sh SEED`), the same seeds, sampler and two
hours per core. The nodes cost more, so the two hours gave 10 dives instead of 20. Logs in
`estimate/x9_600/logs/`; `python3 -I estimate/pool.py estimate/x9_600` rebuilds
`dives.csv`, `levels.csv` and `dive_levels.csv` there (it reads `s41/est.out` and
`s42/est.out`; for section 4's run copy `estimate/logs/est_s*.out` and `dives_s*.txt` to
`estimate/base/s*/est.out` and `dives.txt` first).

| | section 4 (`-x 1:7:1:200:1`, 20 dives) | tuned (`-x 1:9:1:600:1`, 10 dives) |
|---|---|---|
| nodes, mean ± SE | 1.8e10 ± 1.7e10 | **5.1e7 ± 2.2e7** |
| CPU-s, mean ± SE | 1.8e11 ± 1.7e11 (2,110 days on 1,000 cores) | **4.1e9 ± 3.0e9 (47 days)** |
| without the largest dive | 7.6e9 (88 days) | 9.0e8 (10 days) |
| largest dive's share of the CPU | 96% | 77% (seed 42 dive 1: 1.1e8 nodes, 3.1e10 CPU-s) |
| median dive | 9.2e6 nodes, 3.1e8 CPU-s | 7.6e6 nodes, 9.5e8 CPU-s |
| CPU-s per node | 10 | 79 |
| CPU by level | 5: 2%, 6: 5%, 7: 21%, **8: 71%** | 4: 11%, **5: 81%**, 6-8: 8% |

Per level (nodes, mean ± SE; reach in parentheses):

| level | section 4 | tuned |
|---|---|---|
| 3 | 1.7e5 ± 5.5e4 (19) | 9.7e4 ± 4.6e4 (10) |
| 4 | 3.1e6 ± 1.3e6 (16) | 1.9e6 ± 8.8e5 (9) |
| 5 | 1.4e8 ± 1.2e8 (14) | 1.6e7 ± 1.1e7 (6) |
| 6 | 4.0e8 ± 3.3e8 (11) | 1.2e5 ± 1.1e5 (2) |
| 7 | 4.2e8 ± 4.1e8 (7) | 7.5e5 ± 5.0e5 (4) |
| 8 | 1.7e10 ± 1.7e10 (7) | 3.2e7 ± 2.2e7 (2) |

Reading:

- **What the tuning does is clear at the deep levels**: levels 5-8 hold 9x-3,000x fewer
  nodes and fewer dives reach them (6 of 10 reach level 5, 2-4 of 10 levels 6-8, against
  14, 11, 7 and 7 of 20). The later give-up fathoms nodes at depths 4-5 that the
  200-iteration give-up left to branch, which is what the depth-3 sample predicted.
- **What it costs**: 79 CPU-s per node against 10, and the median dive is 3x dearer. The
  CPU moves to level 5 (81%), where nodes that still do not fathom now run 600 iterations.
- **The total**: 4.1e9 CPU-s, about 47 days on 1,000 cores, with one dive at 77% and 10
  dives, so not converged either; 10 days without that dive is a low reading. Against
  section 4, the robust comparison is the deep levels, not the totals: section 4's mean is
  one dive the tuned run did not repeat, and a dive like it would still be possible here.
  The tuned figure is the better estimate of what tai50b costs with a give-up suited to it.
- **Next levers, from where the CPU now is**: a depth-dependent give-up (600 where it buys
  fathoms, less at depth 5 where 81% of the CPU now goes, or the gap-ratio give-up
  `QAP_DNN_GAP_RATIO`), and branching lines chosen for this instance's heavy flows. More
  dives, from an exact depth-2 frontier, before relying on the total.
