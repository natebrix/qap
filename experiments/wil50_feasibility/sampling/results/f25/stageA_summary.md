# Stage A: depths 0-2 of wil50t.dat (IN PROGRESS; 53 tasks)

- binary: commit `8dd01bc3bdd2`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 53 of 53 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,623 + 0 + 12 = 1,634 | 108.9 | 0 | 0 | 635 | 0 |
| 2 | 37 | 7 | 37 | 7 | 11,290 + 20,624 = 31,914 | 3,942 + 7,962 + 9 = 11,914 | 322.0 | 37 | 17 | 1,420 | 0 |
| all | 53 | 7 | 52 | 7 | 15,790 + 20,624 = 36,414 | 5,565 + 7,962 + 22 = 13,549 | | | 17 | | |

Task wall time summed over slots: 13,699 s; CPU/wall 0.989.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 504.4 / 549.5 / 609.4 | 108.8 |
| 2 | hopeless | 37 | 300 / 310 / 310 | 65.0 / 246.8 / 395.2 | 106.2 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 4 | 2 | 4 | 0 | 0 | 10 | 7 | 3 | 0 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 7 | 65.0 / 113.5 / 133.8 | -0.3 / -0.1 / -0.0 | 513 | 1,736 |
| 2 | iter_limit | 1 | 159.1 / 159.1 / 159.1 | 4.8 / 4.8 / 4.8 | 3400 | 1,284 |
| 2 | stalled | 9 | 170.3 / 213.5 / 241.8 | 9.1 / 84.8 / 135.3 | 1110 | 4,942 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 25 (15)
- depth 1: 15 nodes, 15 branched: 15 on a row
  `25:0` row 45 (49); `25:1` row 45 (49); `25:2` row 45 (49); `25:3` row 45 (49); `25:4` row 45 (49); `25:10` row 45 (49); `25:11` row 45 (49); `25:12` row 45 (49); `25:13` row 45 (49); `25:14` row 45 (49); `25:20` row 21 (29); `25:21` row 21 (29); `25:22` row 21 (29); `25:23` row 21 (29); `25:24` row 21 (29)
- depth 2: 37 nodes, 30 branched: 12 on a col, 18 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `25:0` | 609.4 | row 45 | 49 | 1 | 0 | col 9 x1 | 48 |
| `25:1` | 583.3 | row 45 | 49 | 3 | 0 | row 21 x1, col 0 x1, col 9 x1 | 144 |
| `25:2` | 553.9 | row 45 | 49 | 3 | 0 | row 21 x2, col 0 x1 | 144 |
| `25:3` | 546.3 | row 45 | 49 | 2 | 0 | row 21 x1, col 40 x1 | 96 |
| `25:4` | 549.5 | row 45 | 49 | 4 | 2 | row 21 x2 | 96 |
| `25:10` | 602.8 | row 45 | 49 | 2 | 0 | row 21 x2 | 96 |
| `25:11` | 564.2 | row 45 | 49 | 4 | 0 | row 21 x2, col 9 x1, col 0 x1 | 192 |
| `25:12` | 522.7 | row 45 | 49 | 6 | 1 | col 9 x3, row 21 x2 | 240 |
| `25:13` | 512.0 | row 45 | 49 | 0 | 0 |  | 0 |
| `25:14` | 514.8 | row 45 | 49 | 3 | 3 |  | 0 |
| `25:20` | 603.3 | row 21 | 29 | 4 | 0 | col 0 x2, row 45 x2 | 192 |
| `25:21` | 557.8 | row 21 | 29 | 0 | 0 |  | 0 |
| `25:22` | 515.7 | row 21 | 29 | 3 | 1 | row 45 x2 | 96 |
| `25:23` | 504.4 | row 21 | 29 | 1 | 0 | row 45 x1 | 28 |
| `25:24` | 506.0 | row 21 | 29 | 1 | 0 | row 45 x1 | 48 |

Children per depth-2 node: 0: 7 nodes, 28: 1 nodes, 48: 29 nodes
**Depth-2 sample: 37 of 635 nodes, 1,420 children; expected frontier at depth 3: 24,370 nodes (1 SE 1,908)** [e: 635 x the sample's mean children]
Root line forced: facility 25 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 192 | 96 | 192 | 0 | 0 | 460 | 336 | 144 | 0 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 4.8, 10% 36.0, median 266.8, 90% 375.4, max 395.2. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1588.9, median 1908.1, max 2171.1.
Parents kept with their DNN state for stage B: 0 of 37 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.2945995638394044; parents of the 400 that stage B must rebuild: 30

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
