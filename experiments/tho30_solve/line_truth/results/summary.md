# Line score against subtree truth, tho30 depth-1 nodes (2026-11-09)

Recipe `../cutd2_cp.par`; (a) cold first solve of every child, give-up 300, no cut phase, shortfall = cutoff - lb summed over the children where positive; (b) every child's subtree to completion with `-F`, summed plus 1. Ranks: 1 = best (smallest).

## Node 5:2 (rule 9 chose col0)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| col9 | 0.1098 | 6 | 10836.0 | 19/29 | 1 | 864 | 1 | 369 | 2 |
| col29 | 0.1207 | 4 | 13437.1 | 19/29 | 3 | 894 | 2 | 367 | 1 |
| row18 | 0.1076 | 7 | 21027.3 | 8/29 | 5 | 1004 | 3 | 501 | 4 |
| row28 | 0.1463 | 1 | 23599.1 | 6/29 | 6 | 1223 | 4 | 575 | 6 |
| col0 (rule 9) | 0.1225 | 3 | 13223.9 | 18/29 | 2 | 1267 | 5 | 408 | 3 |
| row17 | 0.1051 | 8 | 17157.8 | 14/29 | 4 | 1352 | 6 | 505 | 5 |
| row19 | 0.1105 | 5 | 44545.1 | 0/29 | 8 | 1993 | 7 | 722 | 7 |
| row25 | 0.1299 | 2 | 27643.7 | 14/29 | 7 | 3091 | 8 | 944 | 8 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.71, marginal score -0.21; with (b) CPU: (a) 0.88, score -0.29.
Rule 9's line col0: rank 5 of 8 by (b) nodes, 3 by (b) CPU; the best line col9 has 0.68x its nodes and 0.90x its CPU.

## Node 5:1 (rule 9 chose col0)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| col9 | 0.1229 | 6 | 14589.5 | 18/29 | 1 | 1850 | 1 | 641 | 2 |
| col29 | 0.1343 | 4 | 16680.4 | 18/29 | 3 | 1911 | 2 | 612 | 1 |
| col19 | 0.1170 | 8 | 20549.7 | 15/29 | 5 | 2206 | 3 | 782 | 4 |
| row18 | 0.1215 | 7 | 26578.3 | 6/29 | 6 | 2210 | 4 | 881 | 6 |
| col0 (rule 9) | 0.1548 | 1 | 15537.8 | 18/29 | 2 | 2286 | 5 | 701 | 3 |
| row17 | 0.1293 | 5 | 18839.7 | 14/29 | 4 | 2391 | 6 | 794 | 5 |
| row28 | 0.1537 | 2 | 31568.9 | 0/29 | 8 | 2536 | 7 | 943 | 7 |
| row25 | 0.1374 | 3 | 31256.7 | 11/29 | 7 | 4789 | 8 | 1421 | 8 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.71, marginal score -0.55; with (b) CPU: (a) 0.88, score -0.14.
Rule 9's line col0: rank 5 of 8 by (b) nodes, 3 by (b) CPU; the best line col9 has 0.81x its nodes and 0.91x its CPU.

## Node 5:3 (rule 9 chose col29)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| col29 (rule 9) | 0.1166 | 3 | 8988.1 | 20/29 | 3 | 385 | 1 | 230 | 2 |
| col9 | 0.1063 | 7 | 8137.9 | 21/29 | 1 | 385 | 2 | 225 | 1 |
| col0 | 0.1081 | 5 | 8770.2 | 19/29 | 2 | 493 | 3 | 246 | 3 |
| row18 | 0.1058 | 8 | 14252.1 | 11/29 | 4 | 496 | 4 | 359 | 4 |
| row28 | 0.1469 | 1 | 16709.1 | 12/29 | 5 | 767 | 5 | 466 | 5 |
| row19 | 0.1137 | 4 | 28342.4 | 0/29 | 7 | 966 | 6 | 513 | 6 |
| row25 | 0.1438 | 2 | 22655.4 | 16/29 | 6 | 1812 | 7 | 751 | 7 |
| col14 | 0.1070 | 6 | 36879.5 | 7/29 | 8 | 3084 | 8 | 1219 | 8 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.93, marginal score -0.22; with (b) CPU: (a) 0.95, score -0.26.
Rule 9's line col29: rank 1 of 8 by (b) nodes, 2 by (b) CPU; the best line col29 has 1.00x its nodes and 1.00x its CPU.
