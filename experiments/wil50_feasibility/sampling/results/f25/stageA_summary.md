# Stage A: depths 0-2 of wil50t.dat (IN PROGRESS; 180 tasks)

- binary: commit `8dd01bc3bdd2`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 180 of 180 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,623 + 0 + 12 = 1,634 | 108.9 | 0 | 0 | 635 | 0 |
| 2 | 164 | 34 | 164 | 34 | 50,160 + 116,710 = 166,870 | 17,183 + 43,768 + 50 = 61,001 | 372.0 | 164 | 87 | 6,080 | 0 |
| all | 180 | 34 | 179 | 34 | 54,660 + 116,710 = 171,370 | 18,806 + 43,768 + 62 = 62,636 | | | 87 | | |

Task wall time summed over slots: 63,309 s; CPU/wall 0.989.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 504.4 / 549.5 / 609.4 | 108.8 |
| 2 | hopeless | 164 | 300 / 310 / 330 | 65.0 / 238.6 / 432.5 | 104.7 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 34 | 20 | 18 | 15 | 0 | 0 | 37 | 29 | 8 | 3 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 34 | 65.0 / 122.3 / 156.0 | -0.3 / -0.1 / -0.0 | 722 | 13,363 |
| 2 | iter_limit | 4 | 153.1 / 157.1 / 160.9 | 0.7 / 2.9 / 4.8 | 3400 | 5,074 |
| 2 | stalled | 49 | 156.7 / 208.6 / 242.0 | 8.0 / 84.8 / 135.3 | 1120 | 25,331 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 25 (15)
- depth 1: 15 nodes, 15 branched: 15 on a row
  `25:0` row 45 (49); `25:1` row 45 (49); `25:2` row 45 (49); `25:3` row 45 (49); `25:4` row 45 (49); `25:10` row 45 (49); `25:11` row 45 (49); `25:12` row 45 (49); `25:13` row 45 (49); `25:14` row 45 (49); `25:20` row 21 (29); `25:21` row 21 (29); `25:22` row 21 (29); `25:23` row 21 (29); `25:24` row 21 (29)
- depth 2: 164 nodes, 130 branched: 61 on a col, 69 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `25:0` | 609.4 | row 45 | 49 | 11 | 0 | col 9 x5, col 40 x3, col 49 x2, col 10 x1 | 528 |
| `25:1` | 583.3 | row 45 | 49 | 13 | 0 | row 21 x5, col 9 x3, col 40 x3, col 0 x2 | 624 |
| `25:2` | 553.9 | row 45 | 49 | 18 | 4 | row 21 x8, col 9 x2, col 49 x2, col 0 x1 ... | 672 |
| `25:3` | 546.3 | row 45 | 49 | 10 | 1 | row 21 x5, col 40 x2, col 0 x1, col 49 x1 | 432 |
| `25:4` | 549.5 | row 45 | 49 | 15 | 5 | row 21 x9, col 9 x1 | 480 |
| `25:10` | 602.8 | row 45 | 49 | 12 | 0 | row 21 x5, col 0 x4, col 9 x3 | 576 |
| `25:11` | 564.2 | row 45 | 49 | 12 | 2 | row 21 x4, col 9 x3, col 0 x2, col 49 x1 | 480 |
| `25:12` | 522.7 | row 45 | 49 | 14 | 3 | row 21 x5, col 9 x3, col 49 x2, col 40 x1 | 528 |
| `25:13` | 512.0 | row 45 | 49 | 9 | 3 | col 9 x4, col 40 x1, row 21 x1 | 288 |
| `25:14` | 514.8 | row 45 | 49 | 13 | 6 | row 21 x3, col 49 x2, col 0 x1, col 9 x1 | 336 |
| `25:20` | 603.3 | row 21 | 29 | 9 | 0 | row 45 x6, col 0 x2, row 31 x1 | 372 |
| `25:21` | 557.8 | row 21 | 29 | 6 | 2 | row 45 x3, col 0 x1 | 152 |
| `25:22` | 515.7 | row 21 | 29 | 8 | 4 | row 45 x4 | 192 |
| `25:23` | 504.4 | row 21 | 29 | 8 | 2 | row 45 x6 | 248 |
| `25:24` | 506.0 | row 21 | 29 | 6 | 2 | row 45 x4 | 172 |

Children per depth-2 node: 0: 34 nodes, 28: 8 nodes, 48: 122 nodes
**Depth-2 sample: 164 of 635 nodes, 6,080 children; expected frontier at depth 3: 23,541 nodes (1 SE 831)** [e: 635 x the sample's mean children]
Root line forced: facility 25 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 940 | 864 | 700 | 0 | 0 | 1,696 | 1,372 | 364 | 144 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 0.7, 10% 31.5, median 257.6, 90% 343.7, max 432.5. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1586.4, median 1895.1, max 2191.3.
Parents kept with their DNN state for stage B: 0 of 164 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.06703831182045436; parents of the 400 that stage B must rebuild: 126

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
