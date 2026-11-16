# tai35b and tai40b solved end to end with the cluster-region cut phase (bound_lab session `tai50b_size`, Nov 16 2026)

Task (from the bound-lab coordinator, the test that `../tai50b_size/README.md` names as missing): the same
code and settings that put tai50b at about 100 core-hours, run to the end on two instances whose optima are
known. **[m]** measured (log named), **[e]** estimated.

## Verdict

@FILE:verdict.md@

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

@TREE:bfs35@

The optimum's path, depth: how far the node's certified bound is under the optimum:
@PATH:bfs35@.

The tree is the path plus a short-lived fringe: 3 to 5 nodes branch per level at depths 1 to 4 (the on-path
node, the twin `19:12,18:26` that swaps the first two fixes, and near misses of the phase limit), 2 at depths
5 and 6, and one from depth 7 on.

### tai40b (`runs/bfs40`, `logs/solve_bfs40.txt`)

Root + 732 nodes = **733 nodes, 75,582 DNN iterations, 4,862 CPU-s** (root 1,227 + 3,635); the optimum
637,250,948 found at depth 37; every process ended with exit 0 and no open node. Wall: root 21 minutes,
search 49 minutes on 1 to 2 cores (06:26-07:15).

@TREE:bfs40@

The optimum's path, depth: shortfall: @PATH:bfs40@.

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

@FILE:validity.md@

## 5. A solver fault found on the way: pooled first solves that fail [m]

The first tai35b search blew up at depth 8: 84 nodes left 1,539 open children (`runs/bfs35_unpatched`). Cause:
at depth 7 a node's **pooled first solve** (the parent's cuts frozen in the objective) ended with the certified
value -3.1e11 against an inherited bound of 2.8e8. That is a valid bound and a useless one; being more than
`cut.zone` x |incumbent| below the cutoff it got no phase, the node branched, and its children inherited the
failed state and pool: 29 of the 84 depth-8 first solves and 86 of 216 at depth 9 ended with no certificate
at all (`first lb -inf`). Nothing was fathomed wrongly; the bound was simply lost below that node.

- `pool_retry.patch` (15 lines in `latest/bnb.C`): when a pooled first solve ends unfathomed more than 1% of
  |incumbent| below the bound the node inherited, the node is solved again without the pool, and its phase starts
  from that solve. From depth 7 the tai35b search ran with it: 5 retries, all at depth 7, and in 4 of them the plain solve fathoms
  the node outright (`runs/bfs35/d7_c*/err.txt`).
- It had happened once before, at depth 3 (`15:5,18:12,19:26`: -2.6e8), where the node branched without a phase
  and its children recovered: an avoidable subtree inside the 922 nodes.
- **It is new with `wc=2.88`**: none of the 359 node bounds of the tai50b session (no `wc`) shows it, tai40b
  had none (0 retries in 630 pooled solves from depth 2), tai35b's main search 6 in 894 (and its second route more: section 4). I did not find the mechanism
  (the pool's multipliers were set under the phase's penalty and are reused in a solve at the child's).
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
- A rerun of tai35b's main search from the root with one binary: its depths 1 to 6 ran before the patch
  (with one failed pooled solve at depth 3), depths 7 on after it. The second route ran with the patch throughout.
- tai40b's main search ran depth 1 before the patch; no pooled solve failed there or later.

## Files

| | |
|---|---|
| `par/` | `s35`, `s40` (the searches), `s35h`, `s40c` (second routes), `root35`, `root40`, `r35_*`, `r40_*` (roots by family) |
| `regions/` | the region files used (`comb`, `combh`, `line2`, `line4` per instance) |
| `pool_retry.patch` | section 5; applies after `../regions/regions.patch` and `../tai50b_size/est_resume.patch` |
| `mkpar.sh`, `run_root.sh`, `run_dive.sh`, `bfs.sh`, `run_piece.sh` (unused), `splitstate.py`, `collect.py`, `mkreadme.py` | par files, the root with state, dives, the level-by-level search, state files, the summary, this file |
| `runs/TAG/dD_cK/` | one process of a search: `cmd.txt`, `log.txt`, `nodes.txt` (node trace), `dnn.txt`, `err.txt` (retries); `runs/TAG/progress.txt` the levels with their times |
| `logs/solve_TAG.txt` | `collect.py`'s output: the tree by depth, every node on the optimum's path with its bound |
| `logs/dives35.txt`, `dives40.txt`, `structure.txt`, `*.perm.txt` | the sizing dives, the cluster statistics, the optimal permutations |

State files and parents' DNN states are not kept (`run_root.sh` and `bfs.sh` rewrite them).
