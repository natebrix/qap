# Stage A: depths 0-2 of wil50t.dat (finished; 651 tasks)

- binary: commit `a667727601db`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.04
- symmetry file loaded in 651 of 651 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,762 + 0 + 13 = 1,775 | 118.3 | 0 | 0 | 635 | 0 |
| 2 | 635 | 331 | 635 | 331 | 187,936 + 455,344 = 643,280 | 69,816 + 183,881 + 139 = 253,837 | 399.7 | 635 | 463 | 14,192 | 0 |
| all | 651 | 331 | 650 | 331 | 192,436 + 455,344 = 647,780 | 71,578 + 183,881 + 153 = 255,612 | | | 463 | | |

Task wall time summed over slots: 258,128 s; CPU/wall 0.990.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 333.1 / 460.7 / 570.1 | 117.6 |
| 2 | fathomed | 57 | 23 / 100 / 984 | -15.9 / -0.3 / -0.0 | 39.0 |
| 2 | hopeless | 578 | 300 / 300 / 910 | 4.0 / 146.1 / 457.8 | 112.1 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 331 | 75 | 72 | 42 | 0 | 0 | 54 | 32 | 15 | 13 | 1 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 274 | 4.0 / 83.5 / 170.5 | -3.5 / -0.1 / -0.0 | 393 | 70,645 |
| 2 | iter_limit | 15 | 131.1 / 141.7 / 169.0 | 0.1 / 2.6 / 7.5 | 3400 | 20,509 |
| 2 | stalled | 174 | 134.2 / 190.7 / 243.9 | 6.5 / 69.2 / 144.9 | 1160 | 92,728 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 21 (15)
- depth 1: 15 nodes, 15 branched: 8 on a col, 7 on a row
  `21:0` row 45 (49); `21:1` col 0 (49); `21:2` col 0 (49); `21:3` col 0 (49); `21:4` col 0 (49); `21:10` col 0 (49); `21:11` row 45 (49); `21:12` col 0 (49); `21:13` col 0 (49); `21:14` col 0 (49); `21:20` row 45 (29); `21:21` row 45 (29); `21:22` row 31 (29); `21:23` row 31 (29); `21:24` row 31 (29)
- depth 2: 635 nodes, 304 branched: 183 on a col, 121 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `21:0` | 565.2 | row 45 | 49 | 49 | 17 | col 9 x14, col 40 x12, col 49 x5, col 1 x1 | 1,536 |
| `21:1` | 541.0 | col 0 | 49 | 49 | 13 | row 45 x19, col 9 x8, col 40 x5, col 49 x3 ... | 1,728 |
| `21:2` | 494.6 | col 0 | 49 | 49 | 19 | col 9 x12, row 45 x11, col 49 x7 | 1,440 |
| `21:3` | 440.7 | col 0 | 49 | 49 | 29 | col 9 x12, row 45 x6, col 49 x2 | 960 |
| `21:4` | 402.7 | col 0 | 49 | 49 | 39 | col 9 x9, col 49 x1 | 480 |
| `21:10` | 570.1 | col 0 | 49 | 49 | 11 | col 40 x22, row 45 x15, col 49 x1 | 1,824 |
| `21:11` | 524.8 | row 45 | 49 | 49 | 26 | col 9 x9, col 0 x6, col 49 x5, col 40 x3 | 1,104 |
| `21:12` | 460.7 | col 0 | 49 | 49 | 21 | row 45 x14, col 40 x10, col 49 x3, col 9 x1 | 1,344 |
| `21:13` | 402.6 | col 0 | 49 | 49 | 33 | row 45 x9, col 9 x4, col 40 x2, col 49 x1 | 768 |
| `21:14` | 360.5 | col 0 | 49 | 49 | 44 | col 49 x3, col 9 x2 | 240 |
| `21:20` | 559.2 | row 45 | 29 | 29 | 15 | col 40 x5, col 49 x4, row 31 x3, col 0 x2 | 612 |
| `21:21` | 509.3 | row 45 | 29 | 29 | 19 | col 49 x8, row 34 x1, row 40 x1 | 440 |
| `21:22` | 441.9 | row 31 | 29 | 29 | 5 | row 45 x17, row 40 x6, row 0 x1 | 992 |
| `21:23` | 379.0 | row 31 | 29 | 29 | 13 | row 45 x11, row 40 x4, row 0 x1 | 648 |
| `21:24` | 333.1 | row 31 | 29 | 29 | 27 | row 45 x2 | 76 |

Children per depth-2 node: 0: 331 nodes, 28: 20 nodes, 48: 284 nodes
**Frontier at depth 3: 14,192 nodes** (Burnside orbit count 29,580 for wil50 at depth 3, roadmap/tree_size.md section 1, before any fathoming)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3,500 | 3,336 | 1,976 | 0 | 0 | 2,492 | 1,496 | 720 | 624 | 48 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 0.1, 10% 18.5, median 106.0, 90% 349.8, max 457.8. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1431.4, median 1871.2, max 2270.4.
Parents kept with their DNN state for stage B: 251 of 635 (tau history [[0.01454, 0], [0.0175, 34], [0.04, 77]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.028712535323134263; parents of the 400 that stage B must rebuild: 9

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
