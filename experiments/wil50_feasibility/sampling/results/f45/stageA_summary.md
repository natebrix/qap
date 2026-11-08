# Stage A: depths 0-2 of wil50t.dat (IN PROGRESS; 158 tasks)

- binary: commit `c4f4dc4760f6`, latest/ modified: False, sources newer than the binary: none (built 2026-11-07 18:05:31)
- `bnbdrv -p nug_deep_xbar9_all.par -x 1:4:1:300:1 -s latest/SYM/wil50t.sym -u 48816 -d <depth> -r <one node> -o ...`; environment: fullrun/settings.sh unchanged
- stage-B sample seed 20261107, parents kept below hash 0.0
- symmetry file loaded in 158 of 158 solver logs ("read symmetry info from ...")

| depth | nodes bounded | fathomed | DNN bounded | DNN fathomed | DNN iterations: first solves + cut phases | CPU s: first solves + cut phases + QPB and branching = all | CPU s per node | warm-started | cut phases | children kept | children eliminated in branching |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0 | 0 | 0 + 0 = 0 | 0 + 0 + 1 = 1 | 1.1 | 0 | 0 | 15 | 0 |
| 1 | 15 | 0 | 15 | 0 | 4,500 + 0 = 4,500 | 1,711 + 0 + 12 = 1,723 | 114.9 | 0 | 0 | 635 | 0 |
| 2 | 142 | 85 | 142 | 85 | 45,855 + 94,243 = 140,098 | 16,453 + 37,007 + 25 = 53,485 | 376.7 | 142 | 94 | 2,656 | 0 |
| all | 158 | 85 | 157 | 85 | 50,355 + 94,243 = 144,598 | 18,165 + 37,007 + 38 = 55,209 | | | 94 | | |

Task wall time summed over slots: 55,725 s; CPU/wall 0.991.

Cutoff 48814.05. First DNN solve at each depth (status counts; certified bound at the end of the solve):

| depth | status | nodes | iterations (min / median / max) | shortfall (min / median / max) | CPU s per solve (median) |
|---|---|---|---|---|---|
| 1 | hopeless | 15 | 300 / 300 / 300 | 343.7 / 456.8 / 574.4 | 114.5 |
| 2 | fathomed | 25 | 24 / 300 / 1500 | -1.8 / -0.1 / -0.0 | 122.3 |
| 2 | hopeless | 117 | 300 / 310 / 380 | 22.4 / 128.5 / 444.9 | 108.7 |

Depth-2 nodes by shortfall after their DNN bound (first solve, or the cut phase where one ran; the cut zone is 244):

| fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 85 | 20 | 10 | 4 | 0 | 0 | 9 | 8 | 3 | 3 | 0 | 0 | 0 |

Cut phases (CP form), by outcome (every phase in `phases.tsv`):

| depth | outcome | phases | first-solve shortfall (min / median / max) | after the phase (min / median / max) | iterations (median) | CPU s (sum) |
|---|---|---|---|---|---|---|
| 2 | fathomed | 60 | 22.4 / 74.9 / 139.2 | -1.6 / -0.1 / -0.0 | 349 | 13,910 |
| 2 | iter_limit | 4 | 135.1 / 144.9 / 176.0 | 2.7 / 3.0 / 3.3 | 3400 | 5,443 |
| 2 | stalled | 30 | 127.8 / 170.1 / 228.8 | 6.3 / 43.9 / 126.2 | 1410 | 17,654 |

Branching lines (from the children's assignments: `row i` = facility i of the instance as given placed at each child's location, `col j` = location j):

- depth 0: 1 nodes, 1 branched: 1 on a row
  `root` row 45 (15)
- depth 1: 15 nodes, 15 branched: 3 on a col, 12 on a row
  `45:0` col 49 (49); `45:1` row 21 (49); `45:2` row 21 (49); `45:3` row 21 (49); `45:4` row 21 (49); `45:10` row 21 (49); `45:11` row 21 (49); `45:12` row 21 (49); `45:13` col 0 (49); `45:14` col 0 (49); `45:20` row 21 (29); `45:21` row 21 (29); `45:22` row 21 (29); `45:23` row 31 (29); `45:24` row 31 (29)
- depth 2: 142 nodes, 57 branched: 45 on a col, 12 on a row

Per depth-1 node: its first solve, its line, and its depth-2 children (bounded so far, fathomed, their lines, the frontier nodes they left):

| depth-1 node | shortfall | line | children | bounded | fathomed | depth-2 lines | frontier nodes |
|---|---|---|---|---|---|---|---|
| `45:0` | 574.4 | col 49 | 49 | 8 | 0 | row 21 x6, col 40 x2 | 384 |
| `45:1` | 534.0 | row 21 | 49 | 12 | 6 | col 49 x4, col 0 x1, col 9 x1 | 288 |
| `45:2` | 494.1 | row 21 | 49 | 16 | 12 | col 49 x3, col 1 x1 | 192 |
| `45:3` | 456.8 | row 21 | 49 | 8 | 6 | col 0 x1, col 49 x1 | 96 |
| `45:4` | 432.7 | row 21 | 49 | 13 | 8 | col 40 x2, col 49 x2, col 9 x1 | 240 |
| `45:10` | 560.3 | row 21 | 49 | 11 | 7 | col 0 x4 | 192 |
| `45:11` | 505.6 | row 21 | 49 | 11 | 7 | col 49 x4 | 192 |
| `45:12` | 449.4 | row 21 | 49 | 13 | 8 | col 49 x5 | 240 |
| `45:13` | 406.0 | col 0 | 49 | 7 | 3 | col 49 x3, row 21 x1 | 192 |
| `45:14` | 379.6 | col 0 | 49 | 13 | 10 | col 49 x2, col 40 x1 | 144 |
| `45:20` | 549.6 | row 21 | 29 | 7 | 5 | col 0 x1, col 9 x1 | 96 |
| `45:21` | 491.9 | row 21 | 29 | 6 | 4 | col 9 x1, row 40 x1 | 76 |
| `45:22` | 427.7 | row 21 | 29 | 7 | 5 | col 9 x1, col 49 x1 | 96 |
| `45:23` | 376.2 | row 31 | 29 | 6 | 1 | col 0 x2, row 25 x1, row 40 x1, row 21 x1 | 200 |
| `45:24` | 343.7 | row 31 | 29 | 4 | 3 | row 40 x1 | 28 |

Children per depth-2 node: 0: 85 nodes, 28: 4 nodes, 48: 53 nodes
**Depth-2 sample: 142 of 635 nodes, 2,656 children; expected frontier at depth 3: 11,877 nodes (1 SE 1,084)** [e: 635 x the sample's mean children]
Root line forced: facility 45 in the 15 J1 cells (the rule's line was row 21; out.rule.state)

Frontier nodes by inherited shortfall (cutoff minus the parent's best bound, an upper bound on the node's own shortfall before it is bounded):

| 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 920 | 460 | 192 | 0 | 0 | 432 | 364 | 144 | 144 | 0 | 0 | 0 |

Frontier nodes' inherited bound (the parent's best bound), shortfall: min 2.7, 10% 7.8, median 90.8, 90% 354.9, max 444.9. QPB predicted bound (parent's QPB + U_ij), shortfall: min 1558.9, median 1902.9, max 2265.6.
Parents kept with their DNN state for stage B: 0 of 142 (tau history [[0.0, 0]]: [tau, tasks done when it was set]); the 400th smallest hash is 0.14072210393983334; parents of the 400 that stage B must rebuild: 57

Sanity checks:
- the root: 15 children (expected 15 orbits)
- depth 2 before bounding: 635 candidate children of the 15 depth-1 nodes (635 kept, 0 eliminated in branching; expected 10 x 49 + 5 x 29 = 635 when no depth-1 node fathoms)
