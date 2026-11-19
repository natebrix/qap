# The landscape signal against the line truth sets

## tho30 5:1: slack s = 3,875, the node's own excess 736, 232 children in 8 lines

Child by child, branches iff ex300 <= s: tp 91, fp 52, fn 2 (precision 0.64, recall 0.98); children that branch 93 of 232; median x of the branched 0.45, of the rest 1.07; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 92, fp 77, fn 1

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| col9 | 1850 | 8 | 16 | 29 | 8 | 20 | 29 | 12.3 | 11.8 | 1.46e+04 |
| col29 | 1911 | 9 | 13 | 28 | 9 | 16 | 29 | 12 | 11.6 | 1.67e+04 |
| col19 | 2206 | 11 | 16 | 29 | 9 | 19 | 29 | 12.2 | 11.8 | 2.05e+04 |
| row18 | 2210 | 14 | 19 | 29 | 9 | 26 | 29 | 13.4 | 13 | 2.66e+04 |
| col0 | 2286 | 9 | 14 | 27 | 10 | 15 | 27 | 11.1 | 10.6 | 1.55e+04 |
| row17 | 2391 | 10 | 21 | 29 | 10 | 25 | 29 | 13.9 | 13.7 | 1.88e+04 |
| row28 | 2536 | 18 | 26 | 29 | 11 | 28 | 29 | 14.5 | 14.4 | 3.16e+04 |
| row25 | 4789 | 14 | 18 | 27 | 16 | 20 | 29 | 13.3 | 13.3 | 3.13e+04 |

Spearman with the truth / top line over best: cnt1 +0.62 / 1.03x, cnt1_d0 +0.97 / 1.00x, cnt2 -0.29 / 1.91x, cntg +0.40 / 1.24x, cntg2 -0.08 / 1.24x, soft +0.52 / 1.24x, model +0.62 / 1.24x, a +0.71 / 1.00x

## tho30 5:2: slack s = 3,719, the node's own excess 864, 232 children in 8 lines

Child by child, branches iff ex300 <= s: tp 92, fp 22, fn 4 (precision 0.81, recall 0.96); children that branch 96 of 232; median x of the branched 0.66, of the rest 1.29; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 96, fp 59, fn 0

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| col9 | 864 | 5 | 9 | 28 | 5 | 17 | 28 | 9.95 | 9.17 | 1.08e+04 |
| col29 | 894 | 8 | 9 | 28 | 5 | 12 | 28 | 9.64 | 8.93 | 1.34e+04 |
| row18 | 1004 | 10 | 11 | 29 | 6 | 16 | 29 | 10 | 9.21 | 2.1e+04 |
| row28 | 1223 | 13 | 14 | 29 | 6 | 25 | 29 | 11.8 | 11.1 | 2.36e+04 |
| col0 | 1267 | 8 | 13 | 24 | 8 | 19 | 27 | 9.53 | 8.8 | 1.32e+04 |
| row17 | 1352 | 10 | 15 | 29 | 7 | 20 | 29 | 11.3 | 10.6 | 1.72e+04 |
| row19 | 1993 | 29 | 29 | 29 | 11 | 29 | 29 | 14.6 | 14.5 | 4.45e+04 |
| row25 | 3091 | 13 | 14 | 24 | 6 | 17 | 28 | 10.7 | 10.2 | 2.76e+04 |

Spearman with the truth / top line over best: cnt1 +0.86 / 1.02x, cnt1_d0 +0.72 / 1.02x, cnt2 -0.05 / 2.52x, cntg +0.54 / 1.03x, cntg2 +0.16 / 1.47x, soft +0.52 / 1.47x, model +0.52 / 1.47x, a +0.71 / 1.00x

## tho30 5:3: slack s = 3,537, the node's own excess 1,708, 232 children in 8 lines

Child by child, branches iff ex300 <= s: tp 47, fp 8, fn 34 (precision 0.85, recall 0.58); children that branch 81 of 232; median x of the branched 0.89, of the rest 1.53; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 80, fp 66, fn 1

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| col29 | 385 | 5 | 8 | 26 | 0 | 13 | 28 | 7.77 | 6.89 | 8.99e+03 |
| col9 | 385 | 5 | 4 | 27 | 0 | 15 | 28 | 7.5 | 6.5 | 8.14e+03 |
| col0 | 493 | 5 | 5 | 26 | 0 | 15 | 27 | 7.44 | 6.48 | 8.77e+03 |
| row18 | 496 | 8 | 6 | 29 | 0 | 17 | 29 | 7.93 | 6.9 | 1.43e+04 |
| row28 | 767 | 9 | 7 | 29 | 0 | 17 | 29 | 8.31 | 7.3 | 1.67e+04 |
| row19 | 966 | 19 | 6 | 29 | 0 | 28 | 29 | 9.87 | 8.76 | 2.83e+04 |
| row25 | 1812 | 11 | 6 | 23 | 0 | 17 | 25 | 7.45 | 6.65 | 2.27e+04 |
| col14 | 3084 | 19 | 13 | 29 | 0 | 24 | 29 | 10.4 | 9.62 | 3.69e+04 |

Spearman with the truth / top line over best: cnt1 +0.44 / 1.00x, cnt1_d0 +nan / 2.72x, cnt2 +0.30 / 4.71x, cntg +0.86 / 1.00x, cntg2 +0.28 / 4.71x, soft +0.51 / 1.28x, model +0.65 / 1.28x, a +0.93 / 1.00x

## tho40 24:7,26:1,39:0: slack s = 3,088, the node's own excess 5,954, 296 children in 8 lines

Child by child, branches iff ex300 <= s: tp 0, fp 0, fn 81 (precision 0.00, recall 0.00); children that branch 81 of 296; median x of the branched 2.85, of the rest 3.55; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 43, fp 35, fn 38

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| col32 | 253 | 5 | 0 | 1 | 0 | 11 | 30 | 1.73 | 1.47 | 3.86e+03 |
| col39 | 359 | 6 | 0 | 0 | 0 | 7 | 27 | 1.42 | 1.25 | 6.4e+03 |
| row35 | 431 | 8 | 0 | 0 | 0 | 5 | 31 | 1.43 | 1.26 | 9.88e+03 |
| col15 | 432 | 9 | 0 | 2 | 0 | 13 | 26 | 1.85 | 1.56 | 7.16e+03 |
| row2 | 466 | 8 | 0 | 0 | 0 | 8 | 26 | 1.4 | 1.24 | 6.37e+03 |
| row18 | 578 | 15 | 0 | 0 | 0 | 8 | 36 | 1.51 | 1.32 | 8.79e+03 |
| row14 | 856 | 13 | 0 | 0 | 0 | 10 | 32 | 1.58 | 1.36 | 1.59e+04 |
| row7 | 1175 | 17 | 0 | 0 | 0 | 16 | 35 | 2.01 | 1.67 | 1.93e+04 |

Spearman with the truth / top line over best: cnt1 +nan / 2.25x, cnt1_d0 +nan / 2.25x, cnt2 -0.45 / 2.55x, cntg +0.40 / 1.70x, cntg2 +0.54 / 1.77x, soft +0.31 / 1.84x, model +0.31 / 1.84x, a +0.76 / 1.00x

## tho40 25:32,26:8,39:0: slack s = 2,899, the node's own excess 6,396, 296 children in 8 lines

Child by child, branches iff ex300 <= s: tp 0, fp 0, fn 26 (precision 0.00, recall 0.00); children that branch 26 of 296; median x of the branched 2.81, of the rest 3.73; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 17, fp 62, fn 9

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| row35 | 110 | 2 | 0 | 0 | 0 | 10 | 30 | 1.13 | 1.04 | 493 |
| row2 | 146 | 3 | 0 | 0 | 0 | 7 | 27 | 0.978 | 0.925 | 1.33e+03 |
| row30 | 180 | 2 | 0 | 0 | 0 | 15 | 32 | 1.51 | 1.31 | 1.53e+03 |
| col7 | 182 | 4 | 0 | 0 | 0 | 7 | 24 | 0.969 | 0.916 | 2.72e+03 |
| col39 | 182 | 4 | 0 | 0 | 0 | 8 | 26 | 1.07 | 0.994 | 2.8e+03 |
| col24 | 216 | 3 | 0 | 0 | 0 | 9 | 26 | 1.03 | 0.962 | 2.4e+03 |
| row37 | 252 | 4 | 0 | 0 | 0 | 10 | 30 | 1.28 | 1.15 | 3.08e+03 |
| col16 | 252 | 4 | 0 | 0 | 0 | 13 | 27 | 1.26 | 1.13 | 2.43e+03 |

Spearman with the truth / top line over best: cnt1 +nan / 1.73x, cnt1_d0 +nan / 1.73x, cnt2 +nan / 1.73x, cntg +0.21 / 1.49x, cntg2 -0.23 / 1.65x, soft +0.22 / 1.65x, model +0.22 / 1.65x, a +0.77 / 1.00x

## tho40 26:3,31:7,37:0: slack s = 2,565, the node's own excess 5,020, 296 children in 8 lines

Child by child, branches iff ex300 <= s: tp 0, fp 0, fn 38 (precision 0.00, recall 0.00); children that branch 38 of 296; median x of the branched 2.46, of the rest 3.71; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 34, fp 35, fn 4

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| row2 | 110 | 2 | 0 | 0 | 0 | 5 | 25 | 1.18 | 1.07 | 1.16e+03 |
| col32 | 110 | 2 | 0 | 0 | 0 | 6 | 22 | 1.05 | 0.972 | 1.6e+03 |
| col39 | 182 | 4 | 0 | 0 | 0 | 3 | 17 | 1 | 0.937 | 2.39e+03 |
| col1 | 182 | 4 | 0 | 0 | 0 | 7 | 24 | 1.32 | 1.17 | 2.74e+03 |
| col2 | 253 | 5 | 0 | 0 | 0 | 8 | 29 | 1.31 | 1.17 | 3.46e+03 |
| row30 | 254 | 6 | 0 | 0 | 0 | 10 | 34 | 1.64 | 1.41 | 5.7e+03 |
| row14 | 290 | 7 | 0 | 1 | 0 | 15 | 26 | 1.72 | 1.46 | 7.25e+03 |
| row25 | 361 | 8 | 0 | 0 | 0 | 15 | 31 | 1.75 | 1.48 | 4.57e+03 |

Spearman with the truth / top line over best: cnt1 +nan / 1.98x, cnt1_d0 +nan / 1.98x, cnt2 +0.42 / 1.89x, cntg +0.88 / 1.65x, cntg2 +0.72 / 1.65x, soft +0.86 / 1.65x, model +0.86 / 1.65x, a +0.92 / 1.00x

## tai60b 44:45,43:58: slack s = 1,760,000, the node's own excess 0, 232 children in 4 lines

Child by child, branches iff ex300 <= s: tp 48, fp 15, fn 0 (precision 0.76, recall 1.00); children that branch 48 of 232; median x of the branched 0.27, of the rest 5.22; with the node's excess subtracted (xg = (ex - ex_node) / s <= 1): tp 48, fp 15, fn 0

| line | truth | branched | cnt1 | cnt2 | cnt1_d0 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|---|---|
| r1 | 19600 | 6 | 12 | 15 | 12 | 12 | 15 | 8.27 | 8.32 | nan |
| r45 | 22900 | 10 | 17 | 36 | 15 | 17 | 36 | 15.9 | 15.5 | nan |
| c26 | 38800 | 20 | 20 | 32 | 20 | 20 | 32 | 19.2 | 19.6 | nan |
| r8 | 67100 | 12 | 14 | 14 | 14 | 14 | 14 | 10.6 | 11.1 | nan |

Spearman with the truth / top line over best: cnt1 +0.40 / 1.00x, cnt1_d0 +0.40 / 1.00x, cnt2 -0.40 / 3.42x, cntg +0.40 / 1.00x, cntg2 -0.40 / 3.42x, soft +0.40 / 1.00x, model +0.40 / 1.00x

## Summary

| node | cnt1 | cnt1_d0 | cnt2 | cntg | cntg2 | soft | model | a |
|---|---|---|---|---|---|---|---|---|
| tho30 5:1 | +0.62 / 1.03x | +0.97 / 1.00x | -0.29 / 1.91x | +0.40 / 1.24x | -0.08 / 1.24x | +0.52 / 1.24x | +0.62 / 1.24x | +0.71 / 1.00x |
| tho30 5:2 | +0.86 / 1.02x | +0.72 / 1.02x | -0.05 / 2.52x | +0.54 / 1.03x | +0.16 / 1.47x | +0.52 / 1.47x | +0.52 / 1.47x | +0.71 / 1.00x |
| tho30 5:3 | +0.44 / 1.00x | nan / 2.72x | +0.30 / 4.71x | +0.86 / 1.00x | +0.28 / 4.71x | +0.51 / 1.28x | +0.65 / 1.28x | +0.93 / 1.00x |
| tho40 24:7,26:1,39:0 | nan / 2.25x | nan / 2.25x | -0.45 / 2.55x | +0.40 / 1.70x | +0.54 / 1.77x | +0.31 / 1.84x | +0.31 / 1.84x | +0.76 / 1.00x |
| tho40 25:32,26:8,39:0 | nan / 1.73x | nan / 1.73x | nan / 1.73x | +0.21 / 1.49x | -0.23 / 1.65x | +0.22 / 1.65x | +0.22 / 1.65x | +0.77 / 1.00x |
| tho40 26:3,31:7,37:0 | nan / 1.98x | nan / 1.98x | +0.42 / 1.89x | +0.88 / 1.65x | +0.72 / 1.65x | +0.86 / 1.65x | +0.86 / 1.65x | +0.92 / 1.00x |
| tai60b 44:45,43:58 | +0.40 / 1.00x | +0.40 / 1.00x | -0.40 / 3.42x | +0.40 / 1.00x | -0.40 / 3.42x | +0.40 / 1.00x | +0.40 / 1.00x | - |
| mean, tho | +0.64 / 1.50x (worst 2.25x) | +0.85 / 1.78x (worst 2.72x) | -0.01 / 2.55x (worst 4.71x) | +0.55 / 1.35x (worst 1.70x) | +0.23 / 2.08x (worst 4.71x) | +0.49 / 1.52x (worst 1.84x) | +0.53 / 1.52x (worst 1.84x) | +0.80 / 1.00x (worst 1.00x) |
| mean, tai | +0.40 / 1.00x (worst 1.00x) | +0.40 / 1.00x (worst 1.00x) | -0.40 / 3.42x (worst 3.42x) | +0.40 / 1.00x (worst 1.00x) | -0.40 / 3.42x (worst 3.42x) | +0.40 / 1.00x (worst 1.00x) | +0.40 / 1.00x (worst 1.00x) | - |
