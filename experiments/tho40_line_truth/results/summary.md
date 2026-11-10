# Line score against subtree truth, tho40t depth-3 nodes (2026-11-09)

Recipe `infra/phase1/recipes/nug30.par`; (a) cold first solve of every child, give-up 300, no cut phase, shortfall = cutoff - lb summed over the children where positive; (b) every child's subtree to completion with `-F`, summed plus 1. Ranks: 1 = best (smallest).

## Node 26:3,31:7,37:0 (rule 9 chose col39)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| row2 | 0.1649 | 1 | 1158.6 | 35/37 | 1 | 110 | 1 | 431 | 2 |
| col32 | 0.1237 | 3 | 1604.2 | 35/37 | 2 | 110 | 2 | 415 | 1 |
| col39 (rule 9) | 0.1351 | 2 | 2385.1 | 33/37 | 3 | 182 | 3 | 446 | 3 |
| col1 | 0.1034 | 6 | 2738.0 | 33/37 | 4 | 182 | 4 | 501 | 4 |
| col2 | 0.0952 | 8 | 3456.0 | 32/37 | 5 | 253 | 5 | 572 | 5 |
| row30 | 0.1162 | 4 | 5695.5 | 31/37 | 7 | 254 | 6 | 723 | 6 |
| row14 | 0.0989 | 7 | 7248.2 | 30/37 | 8 | 290 | 7 | 736 | 7 |
| row25 | 0.1089 | 5 | 4566.5 | 29/37 | 6 | 361 | 8 | 776 | 8 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.92, marginal score 0.60; with (b) CPU: (a) 0.90, score 0.62.
Rule 9's line col39: rank 3 of 8 by (b) nodes, 3 by (b) CPU; the best line row2 has 0.60x its nodes and 0.97x its CPU.

## Node 25:32,26:8,39:0 (rule 9 chose col7)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| row35 | 0.0884 | 8 | 493.0 | 35/37 | 1 | 110 | 1 | 650 | 6 |
| row2 | 0.1536 | 1 | 1331.9 | 34/37 | 2 | 146 | 2 | 662 | 7 |
| row30 | 0.0994 | 6 | 1525.6 | 35/37 | 3 | 180 | 3 | 588 | 2 |
| col7 (rule 9) | 0.1311 | 2 | 2721.2 | 33/37 | 6 | 182 | 4 | 603 | 3 |
| col39 | 0.0931 | 7 | 2798.1 | 33/37 | 7 | 182 | 5 | 615 | 4 |
| col24 | 0.1238 | 4 | 2402.8 | 34/37 | 4 | 216 | 6 | 537 | 1 |
| row37 | 0.1276 | 3 | 3082.5 | 33/37 | 8 | 252 | 7 | 770 | 8 |
| col16 | 0.1087 | 5 | 2434.6 | 33/37 | 5 | 252 | 8 | 629 | 5 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.77, marginal score -0.23; with (b) CPU: (a) 0.07, score -0.19.
Rule 9's line col7: rank 4 of 8 by (b) nodes, 3 by (b) CPU; the best line row35 has 0.60x its nodes and 1.08x its CPU.

## Node 24:7,26:1,39:0 (rule 9 chose col39)

| line | marginal score | rank | (a) shortfall sum | fathom cold | rank | (b) nodes | rank | (b) CPU s | rank |
|---|---|---|---|---|---|---|---|---|---|
| col32 | 0.1078 | 5 | 3857.6 | 32/37 | 1 | 253 | 1 | 607 | 1 |
| col39 (rule 9) | 0.1247 | 3 | 6398.0 | 31/37 | 3 | 359 | 2 | 704 | 2 |
| row35 | 0.1582 | 1 | 9884.2 | 29/37 | 6 | 431 | 3 | 834 | 3 |
| col15 | 0.1115 | 4 | 7158.5 | 28/37 | 4 | 432 | 4 | 844 | 4 |
| row2 | 0.1419 | 2 | 6369.5 | 29/37 | 2 | 466 | 5 | 919 | 5 |
| row18 | 0.0861 | 8 | 8787.1 | 22/37 | 5 | 578 | 6 | 1337 | 6 |
| row14 | 0.0940 | 7 | 15922.8 | 24/37 | 7 | 856 | 7 | 1433 | 7 |
| row7 | 0.1025 | 6 | 19300.6 | 20/37 | 8 | 1175 | 8 | 1839 | 8 |

Spearman rank correlation with (b) nodes: (a) shortfall sum 0.76, marginal score 0.55; with (b) CPU: (a) 0.76, score 0.55.
Rule 9's line col39: rank 2 of 8 by (b) nodes, 2 by (b) CPU; the best line col32 has 0.70x its nodes and 0.86x its CPU.
