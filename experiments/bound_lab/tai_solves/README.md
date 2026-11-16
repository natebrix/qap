# tai35b and tai40b solved end to end with the cluster-region cut phase (bound_lab session `tai50b_size`, Nov 16 2026)

Task (from the bound-lab coordinator, the test that `../tai50b_size/README.md` names as missing): the same
code and settings that put tai50b at about 100 core-hours, run to the end on two instances whose optima are
known. **[m]** measured (log named), **[e]** estimated.

## Verdict

| instance | verdict | optimum found and proved | nodes | DNN iterations | CPU-s | wall | first published solve |
|---|---|---|---|---|---|---|---|
| **tai35b** | **solved, valid** | 283,315,445 (= `QAPSOLN/tai35b.sln`), no open node | **922** | 125,232 | 5,275 | 10 min root + 35 min on 3 cores | 2.6 million nodes |
| **tai40b** | **solved, valid** | 637,250,948 (= `QAPSOLN/tai40b.sln`), no open node | **733** | 75,582 | 4,862 | 21 min root + 49 min on 1-2 cores | 278,000 nodes |

- Both from an incumbent one above the optimum, with tai50b's settings plus `dnn.beta_sched = w=0.36,wc=2.88`.
- **No sign of an invalid cut.** On the optima's paths 32 and 37 nodes were bounded and none fathomed; the
  largest certified bounds there are 283,315,444.9994 and 637,250,947.9995. No search ended above the optimum.
- **Each instance was solved a second time by a different route** (another region family, another root, so
  different cuts and pools at every node): tai35b with `combh` from a 1,690-iteration root in 890 nodes, tai40b
  with `comb` from a 1,970-iteration root in 726 nodes. Same optima, no open node, no bound on the paths above
  the optimum (section 4).
- **A solver fault was found and patched on the way** (section 5): with `wc=2.88`, a child's first solve from
  its parent's state and pool sometimes fails outright (a certified value of -3e11, or none). The bound stays
  valid, but the node gets no phase and its children inherit the failure; tai35b's tree blew up at depth 8
  until the failed solve was redone without the pool, and cold if that fails too (`pool_retry.patch`). It did
  not occur on tai40b or in the tai50b session.
  A tai50b run with `wc` needs the patch.
- What it means for tai50b (section 6): the code and the mechanism hold on two instances first solved at
  great cost, and both trees are one path with a small fringe, as the tai50b report assumed. They are easier
  instances than tai50b (roots 0.89% and 0.18% under the optimum against 1.15%), so they do not test its number.

## 1. Setup

- Solver: `latest/` + `../regions/regions.patch` + `../tai50b_size/est_resume.patch`, the build of the tai50b
  session; from the point named in each table, that build + `pool_retry.patch` (section 5), built outside the
  repository. No patch is committed to `latest/`.
- Settings: `../tai50b_size/par/A_wc.par` with the instance's depth range and region file (`mkpar.sh`; `par/s35.par`,
  `par/s40.par`), i.e. tai50b's block with `dnn.beta_sched = w=0.36,wc=2.88`: `dnn.depths = 0:n-2`,
  `dnn.giveup = 100`, `pool.on = 1`, `cut.side = 0`, `cut.zone = 1`, `cut.cp = 100:1000`, `cut.warmy = 1`,
  `cut.rounds = 12`, `cut.iters = 200`, `cut.freeze = 0`, `cut.stall = 3`. The root: the same with
  `cut.stall = 0`, `cut.cp_iters = 6000` (`par/root35.par`, `par/root40.par`).
- Incumbent: the optimum + 1 (`-u 283315446`, `-u 637250949`), so each search has to find the optimum and
  prove it; cutoff = optimum + 0.05.
- The optima, checked: `QAPSOLN/tai35b.sln` and `tai40b.sln` evaluate to 283,315,445 and 637,250,948 with
  `sum A[i,k] B[p(i),p(k)]` (first cell of the setup log; the permutations 0-based in `logs/*.perm.txt`).
- The whole search, `bfs.sh`: level by level. The open nodes at depth D are dealt into state files, each
  bounded by one `bnbdrv -r FILE -d D -o OUT` process (2 to 4 at a time on the 4 cores), and the children left
  open are merged for depth D+1 (`splitstate.py -merge`); at the last level named the files are searched to the
  end. Every node starts warm from its parent with its pool, through the state files and the parent store, as
  under the coordinator. `collect.py` sums the node traces (one line per node) and DNN traces.
- Machine: as `../tai50b_size/logs/machine.txt` (4 cores, Xeon 2.3 GHz). CPU-s are the solver's own; the
  searches shared the cores with each other, so wall times are indicative.

## 2. Regions and roots [m]

Structure (`gen_regions.py -stats`, `logs/structure.txt`): tai35b 3 clusters (15, 2, 18), 84% of the MDS variance
on the first axis; tai40b 3 clusters (12, 9, 19), 4 at the next level, 64% on the first axis: both are 2-D,
unlike tai50b (97%).

Families tried at the root, 1,500 phase iterations after the first solve (`runs/r35_*`, `runs/r40_*`; region
files in `regions/`, built by `../regions/gen_regions.py`):

| family | regions | tai35b: certified at 1,690 iterations, CPU-s | regions | tai40b: certified at 1,970 iterations, CPU-s |
|---|---|---|---|---|
| `comb` = `line:1,top:K,linegap:4:0.03` (tai50b's) | 42 | 279,891,163, 155 | 48 | 633,078,560, 358 |
| `combh` = `comb` + `hier` (every dendrogram cluster) | 74 | 279,907,632, 157 | 85 | **633,363,894**, 369 |
| `line:4` | 126 | 279,741,502, 196 | 146 | 633,071,916, 432 |
| `line:2` | 67 | 279,334,871, 164 | | |

- **tai35b: `comb`**, tai50b's family (`combh` is 16k better, within noise of nothing). **tai40b: `combh`**:
  285k better than `comb` at the same cost, on a 4.2M gap. The other family is each instance's second route
  (section 4).
- Long roots (6,000 phase iterations): **tai35b 280,789,962** (2,525,483 = 0.89% under the optimum; 6,190
  iterations, 583 CPU-s, `runs/root35`); **tai40b 636,110,765** (1,140,183 = 0.18% under; 6,470 iterations,
  1,227 CPU-s, `runs/root40`). With `wc=2.88` the roots are well above the regions session's (279.33M and
  627.37M at 2,300 iterations without it).

## 3. The searches [m]

### tai35b (`runs/bfs35`, `logs/solve_bfs35.txt`)

Root + 921 nodes = **922 nodes, 125,232 DNN iterations, 5,275 CPU-s** (root 583 + 4,692); the optimum
283,315,445 found at depth 32 (`New objective value: 283315445`); every process ended with exit 0 and no
open node. Wall: root 10 minutes, search 35 minutes on 3 cores (05:44-05:54; 06:00-06:28 and 06:38-06:45: depths
7 to 9 were run twice, section 5).

| depth | nodes | branched | fathomed: QPB / DNN first solve / cut phase | enumerated | discarded at pop | DNN bounds | phases | DNN iterations | CPU-s |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 35 | 3 | 0 / 31 / 1 | 0 | 0 | 35 | 4 | 6,509 | 507 |
| 2 | 102 | 3 | 0 / 98 / 1 | 0 | 0 | 102 | 4 | 8,429 | 632 |
| 3 | 99 | 5 | 0 / 90 / 4 | 0 | 0 | 99 | 8 | 15,517 | 927 |
| 4 | 160 | 5 | 0 / 154 / 1 | 0 | 0 | 160 | 6 | 16,816 | 837 |
| 5 | 155 | 2 | 0 / 149 / 4 | 0 | 0 | 155 | 6 | 19,676 | 823 |
| 6 | 60 | 2 | 0 / 58 / 0 | 0 | 0 | 60 | 2 | 6,650 | 218 |
| 7 | 58 | 1 | 0 / 56 / 1 | 0 | 0 | 58 | 2 | 10,080 | 340 |
| 8 | 28 | 1 | 0 / 27 / 0 | 0 | 0 | 28 | 1 | 5,265 | 122 |
| 9 | 27 | 1 | 0 / 26 / 0 | 0 | 0 | 27 | 1 | 2,709 | 47 |
| 10 | 26 | 1 | 0 / 25 / 0 | 0 | 0 | 26 | 1 | 1,923 | 27 |
| 11 | 25 | 1 | 0 / 24 / 0 | 0 | 0 | 25 | 1 | 1,813 | 19 |
| 12 | 24 | 1 | 0 / 23 / 0 | 0 | 0 | 24 | 1 | 1,744 | 16 |
| 13 | 23 | 1 | 0 / 22 / 0 | 0 | 0 | 23 | 1 | 3,519 | 73 |
| 14 | 22 | 1 | 0 / 21 / 0 | 0 | 0 | 22 | 1 | 2,441 | 37 |
| 15 | 19 | 1 | 8 / 10 / 0 | 0 | 2 | 11 | 1 | 1,055 | 18 |
| 16 | 11 | 1 | 6 / 4 / 0 | 0 | 9 | 5 | 1 | 1,293 | 12 |
| 17 | 6 | 1 | 0 / 5 / 0 | 0 | 13 | 6 | 1 | 912 | 9 |
| 18 | 3 | 1 | 0 / 2 / 0 | 0 | 15 | 3 | 1 | 998 | 8 |
| 19 | 2 | 1 | 0 / 1 / 0 | 0 | 15 | 2 | 1 | 850 | 5 |
| 20 | 5 | 1 | 2 / 2 / 0 | 0 | 11 | 3 | 1 | 855 | 4 |
| 21 | 13 | 1 | 7 / 5 / 0 | 0 | 2 | 6 | 1 | 952 | 2 |
| 22 | 4 | 1 | 1 / 2 / 0 | 0 | 10 | 3 | 1 | 1,387 | 4 |
| 23 | 4 | 1 | 1 / 2 / 0 | 0 | 9 | 3 | 1 | 879 | 2 |
| 24 | 1 | 1 | 0 / 0 / 0 | 0 | 10 | 1 | 1 | 1,060 | 1 |
| 25 | 1 | 1 | 0 / 0 / 0 | 0 | 3 | 1 | 1 | 820 | 1 |
| 26 | 1 | 1 | 0 / 0 / 0 | 0 | 0 | 1 | 1 | 830 | 0 |
| 27 | 2 | 1 | 1 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 28 | 1 | 1 | 0 / 0 / 0 | 0 | 0 | 1 | 1 | 820 | 0 |
| 29 | 1 | 1 | 0 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 30 | 1 | 1 | 0 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 31 | 1 | 1 | 0 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 32 | 1 | 0 | 0 / 0 / 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| all | 921 | 45 | 26 / 837 / 12 | 1 | 99 | 894 | 56 | 119,042 | 4,692 |

The optimum's path, depth: how far the node's certified bound is under the optimum:
1: 2,545,608, 2: 1,567,542, 3: 781,191, 4: 582,504, 5: 188,570, 6: 78,249, 7: 45,917, 8: 13,334, 9: 11,070, 10: 11,943, 11: 2,431, 12: 1,700, 13: 0.0, 14: 0.0, 15: 2,278, 16: 0.0, 17: 4,787, 18: 0.0, 19: 2, 20: 52, 21: 9, 22: 0.0, 23: 6, 24: 0.0, 25: 0.2, 26: 0.0, 27: 0.1, 28: 0.0, 29: 0.0, 30: 0.0, 31: 0.0.

The tree is the path plus a short-lived fringe: 3 to 5 nodes branch per level at depths 1 to 4 (the on-path
node, the twin `19:12,18:26` that swaps the first two fixes, and near misses of the phase limit), 2 at depths
5 and 6, and one from depth 7 on.

### tai40b (`runs/bfs40`, `logs/solve_bfs40.txt`)

Root + 732 nodes = **733 nodes, 75,582 DNN iterations, 4,862 CPU-s** (root 1,227 + 3,635); the optimum
637,250,948 found at depth 37; every process ended with exit 0 and no open node. Wall: root 21 minutes,
search 49 minutes on 1 to 2 cores (06:26-07:15).

| depth | nodes | branched | fathomed: QPB / DNN first solve / cut phase | enumerated | discarded at pop | DNN bounds | phases | DNN iterations | CPU-s |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 40 | 1 | 0 / 37 / 2 | 0 | 0 | 40 | 3 | 5,489 | 730 |
| 2 | 39 | 1 | 0 / 38 / 0 | 0 | 0 | 39 | 1 | 3,709 | 405 |
| 3 | 38 | 1 | 0 / 37 / 0 | 0 | 0 | 38 | 1 | 3,242 | 321 |
| 4 | 37 | 1 | 0 / 36 / 0 | 0 | 0 | 37 | 1 | 2,545 | 245 |
| 5 | 36 | 1 | 0 / 35 / 0 | 0 | 0 | 36 | 1 | 2,122 | 200 |
| 6 | 35 | 1 | 0 / 34 / 0 | 0 | 0 | 35 | 1 | 2,115 | 175 |
| 7 | 34 | 1 | 0 / 33 / 0 | 0 | 0 | 34 | 1 | 1,867 | 141 |
| 8 | 33 | 2 | 0 / 30 / 1 | 0 | 0 | 33 | 3 | 6,196 | 429 |
| 9 | 64 | 1 | 0 / 61 / 2 | 0 | 0 | 64 | 3 | 5,864 | 341 |
| 10 | 31 | 1 | 0 / 30 / 0 | 0 | 0 | 31 | 1 | 1,761 | 98 |
| 11 | 30 | 1 | 0 / 29 / 0 | 0 | 0 | 30 | 1 | 1,865 | 99 |
| 12 | 29 | 1 | 0 / 28 / 0 | 0 | 0 | 29 | 1 | 1,923 | 80 |
| 13 | 28 | 1 | 4 / 23 / 0 | 0 | 0 | 24 | 1 | 2,032 | 67 |
| 14 | 27 | 1 | 8 / 18 / 0 | 0 | 0 | 19 | 1 | 1,690 | 53 |
| 15 | 26 | 1 | 0 / 25 / 0 | 0 | 0 | 26 | 1 | 1,956 | 40 |
| 16 | 25 | 1 | 0 / 24 / 0 | 0 | 0 | 25 | 1 | 1,546 | 32 |
| 17 | 24 | 1 | 0 / 23 / 0 | 0 | 0 | 24 | 1 | 1,704 | 31 |
| 18 | 23 | 1 | 6 / 16 / 0 | 0 | 0 | 17 | 1 | 1,467 | 27 |
| 19 | 22 | 1 | 14 / 7 / 0 | 0 | 0 | 8 | 1 | 1,857 | 32 |
| 20 | 13 | 1 | 6 / 6 / 0 | 0 | 0 | 7 | 1 | 2,096 | 32 |
| 21 | 14 | 1 | 0 / 13 / 0 | 0 | 0 | 14 | 1 | 1,140 | 12 |
| 22 | 17 | 1 | 4 / 12 / 0 | 0 | 0 | 13 | 1 | 1,108 | 11 |
| 23 | 18 | 1 | 10 / 7 / 0 | 0 | 0 | 8 | 1 | 959 | 8 |
| 24 | 11 | 1 | 5 / 5 / 0 | 0 | 0 | 6 | 1 | 956 | 4 |
| 25 | 11 | 1 | 0 / 10 / 0 | 0 | 0 | 11 | 1 | 1,463 | 8 |
| 26 | 4 | 1 | 0 / 3 / 0 | 0 | 6 | 4 | 1 | 1,330 | 6 |
| 27 | 1 | 1 | 0 / 0 / 0 | 0 | 8 | 1 | 1 | 960 | 3 |
| 28 | 5 | 1 | 0 / 4 / 0 | 0 | 3 | 5 | 1 | 920 | 2 |
| 29 | 3 | 1 | 0 / 2 / 0 | 0 | 0 | 3 | 1 | 840 | 1 |
| 30 | 2 | 1 | 0 / 1 / 0 | 0 | 3 | 2 | 1 | 1,070 | 1 |
| 31 | 2 | 1 | 0 / 1 / 0 | 0 | 0 | 2 | 1 | 1,270 | 1 |
| 32 | 3 | 1 | 2 / 0 / 0 | 0 | 1 | 1 | 1 | 810 | 0 |
| 33 | 1 | 1 | 0 / 0 / 0 | 0 | 2 | 1 | 1 | 810 | 0 |
| 34 | 1 | 1 | 0 / 0 / 0 | 0 | 1 | 1 | 1 | 810 | 0 |
| 35 | 2 | 1 | 1 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 36 | 2 | 1 | 1 / 0 / 0 | 0 | 0 | 1 | 1 | 810 | 0 |
| 37 | 1 | 0 | 0 / 0 / 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| all | 732 | 37 | 61 / 628 / 5 | 1 | 24 | 670 | 42 | 69,112 | 3,635 |

The optimum's path, depth: shortfall: 1: 368,475, 2: 261,654, 3: 209,374, 4: 315,406, 5: 216,367, 6: 165,298, 7: 137,943, 8: 108,875, 9: 78,701, 10: 57,563, 11: 36,379, 12: 50,892, 13: 35,137, 14: 378,360, 15: 148,329, 16: 110,722, 17: 42,491, 18: 19,381, 19: 5,903, 20: 6, 21: 805, 22: 4,025, 23: 21, 24: 5, 25: 0.0, 26: 0.1, 27: 0.0, 28: 0.0, 29: 0.1, 30: 0.0, 31: 0.0, 32: 0.0, 33: 0.0, 34: 0.0, 35: 0.0, 36: 0.0.

One path: one node branches per level (two at depth 8), and its siblings fathom in the pooled first solve in
40 to 110 iterations.

### What the sizing dives predicted

| | dives | predicted below the root: nodes, CPU-s | the search took | 
|---|---|---|---|
| tai35b | 36 from the root's 35 children (`runs/d35`, `logs/dives35.txt`) | 36, 272 | 921, 4,692 |
| tai40b | 25 from the root's 40 children (`runs/d40`, `logs/dives40.txt`) | 39, 328 | 732, 3,635 |

Both said "fits easily", which was all they were for, and both were 20 to 25 times low in nodes: 60 of the 61
dives drew a child that fathoms at once, and the one that drew a branching child then drew a fathoming
grandchild. A tree that is one path is invisible to a few dozen dives from the top. This is the effect the
tai50b report stratified its dives against (its sets by spine level); its unstratified estimate E1 should be
read with that in mind.

## 4. Validity [m]

- **The answers.** Both searches start at the optimum + 1, report `New objective value` = the optimum
  (283315445 in `runs/bfs35/d14_c*/log.txt`, 637250948 in `runs/bfs40/d16_c*/log.txt`), end with every process at
  exit 0 and no `node` line in the last level's state files. The values of the `QAPSOLN` permutations were
  recomputed from the `.dat` files and equal them.
- **The path check** (`collect.py`, in `logs/solve_bfs35.txt` and `logs/solve_bfs40.txt` with every node): a node
  all of whose fixes agree with the optimal permutation must not be fathomed, and no DNN or cut bound there may
  exceed the optimum. tai35b: 32 nodes on the path (depths 1 to 32), all branched or enumerated, largest
  certified bound 283,315,444.9994. tai40b: 37 nodes, largest 637,250,947.9995. **0 violations.** The roots:
  280,789,962 and 636,110,765.
- The check is sharp where it matters: from depth 13 (tai35b) and 20 (tai40b) the bounds on the path are within
  10 of the optimum, most within 0.1, and never above it, through 30 more solves each with an inherited pool.
- **A second route for each** (`runs/bfs35h`, `runs/bfs40c`; `logs/solve_bfs35h.txt`, `logs/solve_bfs40c.txt`),
  the other region family from its own short root (the 1,500-phase-iteration roots of section 2), so the root's
  pool, every inherited pool and every phase differ:

  | search | regions | root | nodes | DNN iterations | CPU-s | optimum found, no open node | path: nodes, violations, largest bound |
  |---|---|---|---|---|---|---|---|
  | tai35b main | `comb` (42) | 280,789,962 | 922 | 125,232 | 5,275 | yes | 32, 0, 283,315,444.9994 |
  | tai35b second | `combh` (74) | 279,907,632 | 890 | 112,396 | 4,882 | yes | 32, 0, 283,315,444.9999 |
  | tai40b main | `combh` (85) | 636,110,765 | 733 | 75,582 | 4,862 | yes | 37, 0, 637,250,947.9995 |
  | tai40b second | `comb` (48) | 633,078,560 | 726 | 69,745 | 4,551 | yes | 37, 0, 637,250,947.9998 |

  The second routes are no bigger than the main ones although their roots are 0.9M and 3.0M lower: the
  6,000-iteration roots did not pay for themselves on trees this small (tai40b: 1,227 against 359 CPU-s at the root).
- Not a check of every cut: a cut that removed a solution worse than the optimum only would change node
  counts, not answers, and would not show here. What these runs exclude is a cut that removes the optimum on
  these two instances with these regions.

## 5. A solver fault found on the way: pooled first solves that fail [m]

The first tai35b search blew up at depth 8: 84 nodes left 1,539 open children (`runs/bfs35_unpatched`). Cause:
at depth 7 a node's **pooled first solve** (the parent's cuts frozen in the objective) ended with the certified
value -3.1e11 against an inherited bound of 2.8e8. That is a valid bound and a useless one; being more than
`cut.zone` x |incumbent| below the cutoff it got no phase, the node branched, and its children inherited the
failed state and pool: 29 of the 84 depth-8 first solves and 86 of 216 at depth 9 ended with no certificate
at all (`first lb -inf`). Nothing was fathomed wrongly; the bound was simply lost below that node.

- `pool_retry.patch` (30 lines in `latest/bnb.C`), two steps. When a pooled first solve ends unfathomed more than
  1% of |incumbent| below the bound the node inherited, the node is solved again without the pool; when that
  solve, or any warm first solve, fails the same way, the node is solved cold. The phase then starts from the
  solve that worked. tai35b's main search ran with the first step from depth 7: 5 retries, all at depth 7, and
  in 4 of them the plain solve fathoms the node outright (`runs/bfs35/d7_c*/err.txt`). Its second route showed
  the first step is not enough: at depth 7 two retried solves failed again (-inf, -2.2e11) and the tree blew
  up the same way (84 nodes at depth 8 left 624 open, `runs/bfs35h_retry1`); with the cold step both nodes
  fathom (283,538,000 and 283,317,000 against the cutoff 283,315,445.05) and the search ends in 890 nodes
  (5 pool retries at depths 3, 4, 7 and 2 cold ones).
- It had happened once before, at depth 3 (`15:5,18:12,19:26`: -2.6e8), where the node branched without a phase
  and its children recovered: an avoidable subtree inside the 922 nodes.
- **It is new with `wc=2.88`**: none of the 359 node bounds of the tai50b session (no `wc`) shows it. tai40b had
  none in either search (0 retries in about 1,300 warm solves). tai35b: 6 in 894 in the main search, 5 in 862 in
  the second. I did not find the mechanism (the failures are in children of nodes whose phase ran long at the
  phase's penalty; the child's solve starts from that state at its own penalty, and conditioning the parent's
  Y on an assignment it gives almost no weight is the other suspect).
  **A tai50b run with `wc=2.88` needs this patch or its equivalent**; without `wc` it costs phase iterations
  instead (`../tai50b_size`, section 9b).

## 6. What this says about the tai50b estimate

- The code and settings prove both instances, by two region families each, with every bound on the optima's
  paths at or below the optimum. That removes the report's caveat 5 as far as two instances can.
- **Node counts against the first published solves** (Fujii et al., as the coordinator quotes them): tai35b
  922 against 2.6 million, tai40b 733 against 278,000.
- The trees have the shape the tai50b report assumed: the optimum's path, on which the bound's gap falls by
  about half per level at first (tai35b: 2.5M, 1.6M, 0.78M, 0.58M, 0.19M at depths 1 to 5) and is exact from
  about depth 13 to 25; a fringe of nodes holding near-optimal solutions that is a few nodes wide for a few
  levels and gone by depth 7 (tai35b) or nearly absent (tai40b); and siblings that fathom in the first solve.
  Neither tree compounds.
- They are easier than tai50b in the way that matters: the roots end 0.89% and 0.18% under the optimum,
  tai50b's 1.15% under its best known value with a flatter landscape (150 2-swaps within 1M). So these runs
  support the mechanism and the code, not the number: tai50b's 60 to 250 core-hours still rests on its own
  thinly sampled levels 4 and 5.

## 7. Not done

- tai50b itself.
- A search with the cut phase switched off below some depth (the other "different route" the task offered); the
  second routes here change the region family, the root and every pool.
- The mechanism of the failed pooled solves, and `pool_retry.patch` on anything but these two instances.
- A rerun of the main searches from the root with one binary. tai35b's main search ran depths 1 to 6 before the
  patch (one failed pooled solve at depth 3, which cost a subtree) and depths 7 on with its first step; tai40b's
  ran depth 1 before it, and nothing failed there or later. tai35b's second route ran depths 1 to 6 with the
  first step and 7 on with both; tai40b's second route with the first step throughout (0 retries).

## Files

| | |
|---|---|
| `verdict.md`, `verdict_second.md`, `validity.md`, `validity_second.md`, `README.tmpl.md` | the hand-written parts of this file |
| `par/` | `s35`, `s40` (the searches), `s35h`, `s40c` (second routes), `root35`, `root40`, `r35_*`, `r40_*` (roots by family) |
| `regions/` | the region files used (`comb`, `combh`, `line2`, `line4` per instance) |
| `pool_retry.patch` | section 5; applies after `../regions/regions.patch` and `../tai50b_size/est_resume.patch` |
| `mkpar.sh`, `run_root.sh`, `run_dive.sh`, `bfs.sh`, `run_piece.sh` (unused), `splitstate.py`, `collect.py`, `mkreadme.py` | par files, the root with state, dives, the level-by-level search, state files, the summary, this file |
| `runs/TAG/dD_cK/` | one process of a search: `cmd.txt`, `log.txt`, `nodes.txt` (node trace), `dnn.txt`, `err.txt` (retries); `runs/TAG/progress.txt` the levels with their times; `runs/bfs35_unpatched`, `runs/bfs35h_retry1`: the levels that blew up (section 5) |
| `logs/solve_TAG.txt` | `collect.py`'s output: the tree by depth, every node on the optimum's path with its bound |
| `logs/dives35.txt`, `dives40.txt`, `structure.txt`, `*.perm.txt` | the sizing dives, the cluster statistics, the optimal permutations |

State files and parents' DNN states are not kept (`run_root.sh` and `bfs.sh` rewrite them).
