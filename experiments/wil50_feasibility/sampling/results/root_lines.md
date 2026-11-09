# Root lines: facility 21 (the rule's line, full stage A) against forced lines

## Depth-1 first solves against the laptop's children

Ours: the production first solve (QPB warm start from the root, cold DNN, give-up at 300 iterations when hopeless; Xeon, OpenBLAS). Laptop: the same children cold with `-i 300` (`tuning/lines/root/children_i300/`, rev 165780c) and `-i 1,500` (`children/`, rev 2d26817). Bound (shortfall to the cutoff 48,814.05).

| cell | f21 ours, 300 it | f21 laptop, 300 it | f21 laptop, 1,500 it | f21 line, children | f25 ours, 300 it | f25 laptop, 300 it | f25 laptop, 1,500 it | f25 line, children | f45 ours, 300 it | f45 laptop, 300 it | f45 laptop, 1,500 it | f45 line, children |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 48248.9 (565.2) | 48249.1 (565.0) | 48305.8 (508.3) | row 45, 49 | 48204.7 (609.4) | 48204.8 (609.2) | 48255.1 (558.9) | row 45, 49 | 48239.6 (574.4) | 48239.9 (574.2) | 48298.9 (515.1) | col 49, 49 |
| 1 | 48273.1 (541.0) | 48273.4 (540.7) | 48328.1 (486.0) | col 0, 49 | 48230.8 (583.3) | 48230.9 (583.1) | 48281.5 (532.6) | row 45, 49 | 48280.0 (534.0) | 48280.2 (533.8) | 48335.2 (478.8) | row 21, 49 |
| 2 | 48319.4 (494.6) | 48319.9 (494.1) | 48379.6 (434.4) | col 0, 49 | 48260.1 (553.9) | 48260.3 (553.7) | 48310.5 (503.5) | row 45, 49 | 48320.0 (494.1) | 48320.4 (493.6) | 48383.7 (430.3) | row 21, 49 |
| 3 | 48373.3 (440.7) | 48373.7 (440.4) | 48435.0 (379.1) | col 0, 49 | 48267.8 (546.3) | 48267.8 (546.3) | 48320.1 (493.9) | row 45, 49 | 48357.3 (456.8) | 48357.5 (456.5) | 48425.6 (388.4) | row 21, 49 |
| 4 | 48411.4 (402.7) | 48411.5 (402.5) | 48470.4 (343.7) | col 0, 49 | 48264.5 (549.5) | 48264.5 (549.5) | 48319.0 (495.1) | row 45, 49 | 48381.3 (432.7) | 48381.5 (432.5) | 48449.5 (364.5) | row 21, 49 |
| 10 | 48244.0 (570.1) | 48244.2 (569.9) | 48300.6 (513.5) | col 0, 49 | 48211.3 (602.8) | 48211.5 (602.6) | 48262.8 (551.2) | row 45, 49 | 48253.8 (560.3) | 48254.0 (560.1) | 48310.3 (503.7) | row 21, 49 |
| 11 | 48289.2 (524.8) | 48289.5 (524.5) | 48347.2 (466.9) | row 45, 49 | 48249.9 (564.2) | 48250.0 (564.0) | 48303.2 (510.9) | row 45, 49 | 48308.4 (505.6) | 48308.7 (505.4) | 48364.7 (449.3) | row 21, 49 |
| 12 | 48353.3 (460.7) | 48353.6 (460.5) | 48411.1 (403.0) | col 0, 49 | 48291.4 (522.7) | 48291.5 (522.6) | 48342.6 (471.5) | row 45, 49 | 48364.7 (449.4) | 48365.0 (449.0) | 48423.5 (390.5) | row 21, 49 |
| 13 | 48411.5 (402.6) | 48411.9 (402.2) | 48472.7 (341.3) | col 0, 49 | 48302.0 (512.0) | 48302.1 (511.9) | 48355.1 (459.0) | row 45, 49 | 48408.0 (406.0) | 48408.4 (405.7) | 48471.9 (342.2) | col 0, 49 |
| 14 | 48453.6 (360.5) | 48453.8 (360.2) | 48512.4 (301.6) | col 0, 49 | 48299.2 (514.8) | 48299.4 (514.7) | 48354.6 (459.4) | row 45, 49 | 48434.4 (379.6) | 48434.8 (379.3) | 48500.1 (314.0) | col 0, 49 |
| 20 | 48254.9 (559.2) | 48254.8 (559.3) | 48309.1 (505.0) | row 45, 29 | 48210.7 (603.3) | 48210.9 (603.2) | 48261.8 (552.2) | row 21, 29 | 48264.5 (549.6) | 48264.4 (549.6) | 48319.4 (494.7) | row 21, 29 |
| 21 | 48304.7 (509.3) | 48304.6 (509.4) | 48361.6 (452.5) | row 45, 29 | 48256.2 (557.8) | 48256.3 (557.7) | 48307.3 (506.8) | row 21, 29 | 48322.2 (491.9) | 48322.0 (492.0) | 48377.9 (436.2) | row 21, 29 |
| 22 | 48372.1 (441.9) | 48372.1 (441.9) | 48429.5 (384.6) | row 31, 29 | 48298.3 (515.7) | 48298.5 (515.6) | 48348.7 (465.4) | row 21, 29 | 48386.3 (427.7) | 48386.5 (427.6) | 48441.0 (373.1) | row 21, 29 |
| 23 | 48435.1 (379.0) | 48435.2 (378.8) | 48494.5 (319.5) | row 31, 29 | 48309.7 (504.4) | 48309.8 (504.2) | 48362.3 (451.8) | row 21, 29 | 48437.8 (376.2) | 48438.0 (376.1) | 48493.9 (320.1) | row 31, 29 |
| 24 | 48480.9 (333.1) | 48481.0 (333.0) | 48537.9 (276.2) | row 31, 29 | 48308.1 (506.0) | 48308.2 (505.9) | 48361.9 (452.2) | row 21, 29 | 48470.4 (343.7) | 48470.5 (343.6) | 48525.9 (288.2) | row 31, 29 |

Lines cost over the 15 children, sum 47^(gap/130), and the ratio to facility 45's:

| facility | ours, 300 it | laptop, 300 it | laptop, 1,500 it | largest difference ours - laptop at 300 it |
|---|---|---|---|---|
| 21 | 7.85e+07, 1.13x | 7.81e+07, 1.13x | 1.48e+07, 1.17x | 0.48 |
| 25 | 3.07e+08, 4.42x | 3.06e+08, 4.42x | 6.72e+07, 5.28x | 0.21 |
| 45 | 6.96e+07, 1.00x | 6.93e+07, 1.00x | 1.27e+07, 1.00x | 0.41 |

## Depth 2

The sample: the 200 depth-2 task ids with the smallest seeded order hash (`stage_a.py --sample-last`), the same ids in every run. Lines cost: 635 x mean over the nodes of 47^(max(gap, 0)/130), gap = the node's shortfall after its first solve and cut phase (a fathomed node counts 1), with a 95% bootstrap interval (4,000 resamples); expected frontier: 635 x the mean children (± 1 SE).

| line | depth-2 nodes | fathomed (fraction) | cut phases | expected frontier | lines cost (95%) | lines cost, a miss at least 1 + children (95%) | CPU s per node: first solve + phase = all |
|---|---|---|---|---|---|---|---|
| f21, all (population) | 635 of 635 | 331 (0.521; 1-sigma 0.501-0.541) | 463 | 14,192 ± 0 | 5.67e+06 (3.49e+06-8.29e+06) | 5.67e+06 (3.5e+06-8.29e+06) | 110 + 290 = 400 |
| f21, matched sample | 200 of 635 | 102 (0.510; 1-sigma 0.475-0.545) | 146 | 14,491 ± 876 | 6.42e+06 (3.25e+06-1.02e+07) | 6.43e+06 (3.26e+06-1.02e+07) | 112 + 281 = 393 |
| f25, sample | 200 of 635 | 37 (0.185; 1-sigma 0.159-0.214) | 106 | 24,206 ± 694 | 6.79e+06 (3.92e+06-1.05e+07) | 6.8e+06 (3.93e+06-1.05e+07) | 105 + 271 = 376 |
| f45, sample | 200 of 635 | 113 (0.565; 1-sigma 0.530-0.600) | 137 | 12,814 ± 866 | 6.32e+06 (1.7e+06-1.25e+07) | 6.33e+06 (1.7e+06-1.25e+07) | 114 + 265 = 379 |

Shortfall histogram of the depth-2 nodes (after the first solve and cut phase; the zone is 244):

| nodes | fathomed | 0-50 | 50-100 | 100-150 | 150-200 | 200-244 | 244-300 | 300-350 | 350-400 | 400-450 | 450-500 | 500-600 | >= 600 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f21 all | 331 | 75 | 72 | 42 | 0 | 0 | 54 | 32 | 15 | 13 | 1 | 0 | 0 |
| f21 matched sample | 102 | 26 | 18 | 14 | 0 | 0 | 16 | 10 | 9 | 5 | 0 | 0 | 0 |
| f25 sample | 37 | 25 | 27 | 17 | 0 | 0 | 50 | 32 | 8 | 4 | 0 | 0 | 0 |
| f45 sample | 113 | 27 | 19 | 8 | 0 | 0 | 17 | 9 | 3 | 4 | 0 | 0 | 0 |

## The line at each depth-1 node and its children

Rows-only nodes (the cell in the grid's middle row, 20-24: the J2 symmetry) branch by QPB row lookahead (`select_row_lookahead`, no DNN marginals); the others by rule 9's lookahead between the DNN-marginal row and column (`select_ind_lookahead` with `xbar_lines`). Children: the depth-2 nodes bounded (all for the full run, the sample's otherwise).

| line | depth-1 node | how | line | children | bounded | fathomed | gap min / median / max |
|---|---|---|---|---|---|---|---|
| f21 | `21:0` | rule 9 xbar lines | row 45 | 49 | 49 | 17 | 5 / 199 / 426 |
| f21 | `21:1` | rule 9 xbar lines | col 0 | 49 | 49 | 13 | 4 / 250 / 431 |
| f21 | `21:2` | rule 9 xbar lines | col 0 | 49 | 49 | 19 | 11 / 105 / 381 |
| f21 | `21:3` | rule 9 xbar lines | col 0 | 49 | 49 | 29 | 3 / 63 / 319 |
| f21 | `21:4` | rule 9 xbar lines | col 0 | 49 | 49 | 39 | 32 / 75 / 272 |
| f21 | `21:10` | rule 9 xbar lines | col 0 | 49 | 49 | 11 | 19 / 269 / 458 |
| f21 | `21:11` | rule 9 xbar lines | row 45 | 49 | 49 | 26 | 0 / 136 / 391 |
| f21 | `21:12` | rule 9 xbar lines | col 0 | 49 | 49 | 21 | 6 / 92 / 348 |
| f21 | `21:13` | rule 9 xbar lines | col 0 | 49 | 49 | 33 | 0 / 26 / 275 |
| f21 | `21:14` | rule 9 xbar lines | col 0 | 49 | 49 | 44 | 19 / 47 / 80 |
| f21 | `21:20` | rows-only QPB lookahead | row 45 | 29 | 29 | 15 | 3 / 257 / 422 |
| f21 | `21:21` | rows-only QPB lookahead | row 45 | 29 | 29 | 19 | 17 / 262 / 365 |
| f21 | `21:22` | rows-only QPB lookahead | row 31 | 29 | 29 | 5 | 14 / 115 / 299 |
| f21 | `21:23` | rows-only QPB lookahead | row 31 | 29 | 29 | 13 | 2 / 56 / 95 |
| f21 | `21:24` | rows-only QPB lookahead | row 31 | 29 | 29 | 27 | 0 / 1 / 1 |
| f25 | `25:0` | rule 9 xbar lines | row 45 | 49 | 15 | 0 | 81 / 290 / 416 |
| f25 | `25:1` | rule 9 xbar lines | row 45 | 49 | 17 | 0 | 5 / 258 / 384 |
| f25 | `25:2` | rule 9 xbar lines | row 45 | 49 | 20 | 4 | 1 / 256 / 348 |
| f25 | `25:3` | rule 9 xbar lines | row 45 | 49 | 13 | 1 | 9 / 255 / 324 |
| f25 | `25:4` | rule 9 xbar lines | row 45 | 49 | 16 | 5 | 54 / 267 / 350 |
| f25 | `25:10` | rule 9 xbar lines | row 45 | 49 | 15 | 0 | 8 / 245 / 433 |
| f25 | `25:11` | rule 9 xbar lines | row 45 | 49 | 17 | 2 | 7 / 133 / 375 |
| f25 | `25:12` | rule 9 xbar lines | row 45 | 49 | 18 | 4 | 34 / 191 / 319 |
| f25 | `25:13` | rule 9 xbar lines | row 45 | 49 | 10 | 3 | 24 / 256 / 313 |
| f25 | `25:14` | rule 9 xbar lines | row 45 | 49 | 17 | 7 | 2 / 98 / 308 |
| f25 | `25:20` | rows-only QPB lookahead | row 21 | 29 | 10 | 0 | 7 / 286 / 415 |
| f25 | `25:21` | rows-only QPB lookahead | row 21 | 29 | 6 | 2 | 14 / 289 / 326 |
| f25 | `25:22` | rows-only QPB lookahead | row 21 | 29 | 9 | 4 | 109 / 276 / 306 |
| f25 | `25:23` | rows-only QPB lookahead | row 21 | 29 | 9 | 3 | 30 / 171 / 258 |
| f25 | `25:24` | rows-only QPB lookahead | row 21 | 29 | 8 | 2 | 39 / 168 / 294 |
| f45 | `45:0` | rule 9 xbar lines | col 49 | 49 | 15 | 1 | 29 / 309 / 449 |
| f45 | `45:1` | rule 9 xbar lines | row 21 | 49 | 17 | 8 | 22 / 70 / 358 |
| f45 | `45:2` | rule 9 xbar lines | row 21 | 49 | 20 | 13 | 41 / 244 / 348 |
| f45 | `45:3` | rule 9 xbar lines | row 21 | 49 | 13 | 7 | 54 / 82 / 271 |
| f45 | `45:4` | rule 9 xbar lines | row 21 | 49 | 16 | 9 | 3 / 32 / 78 |
| f45 | `45:10` | rule 9 xbar lines | row 21 | 49 | 15 | 10 | 24 / 77 / 292 |
| f45 | `45:11` | rule 9 xbar lines | row 21 | 49 | 17 | 10 | 8 / 279 / 337 |
| f45 | `45:12` | rule 9 xbar lines | row 21 | 49 | 18 | 11 | 14 / 137 / 305 |
| f45 | `45:13` | rule 9 xbar lines | col 0 | 49 | 10 | 6 | 10 / 62 / 91 |
| f45 | `45:14` | rule 9 xbar lines | col 0 | 49 | 17 | 14 | 6 / 29 / 53 |
| f45 | `45:20` | rows-only QPB lookahead | row 21 | 29 | 10 | 7 | 126 / 257 / 402 |
| f45 | `45:21` | rows-only QPB lookahead | row 21 | 29 | 6 | 4 | 3 / 162 / 321 |
| f45 | `45:22` | rows-only QPB lookahead | row 21 | 29 | 9 | 6 | 64 / 135 / 259 |
| f45 | `45:23` | rows-only QPB lookahead | row 31 | 29 | 9 | 2 | 3 / 43 / 89 |
| f45 | `45:24` | rows-only QPB lookahead | row 31 | 29 | 8 | 5 | 3 / 27 / 39 |
