- tho30 5:1: best line col9 1,850 nodes, worst 2.59x, rule 9's col0 1.24x; QPB children sum (P = 1): col0 1.24x (+0.48); exp, L = 32: col9 1.00x (+0.74)
- tho30 5:2: best line col9 864 nodes, worst 3.58x, rule 9's col0 1.47x; QPB children sum (P = 1): col0 1.47x (+0.19); exp, L = 32: col9 1.00x (+0.52)
- tho30 5:3: best line col29 385 nodes, worst 8.01x, rule 9's col29 1.00x; QPB children sum (P = 1): col0 1.28x (+0.71); exp, L = 32: col9 1.00x (+0.67)
- tho40 24:7,26:1,39:0: best line col32 253 nodes, worst 4.64x, rule 9's col39 1.42x; QPB children sum (P = 1): col32 1.00x (+0.83); exp, L = 32: col32 1.00x (+0.83)
- tho40 25:32,26:8,39:0: best line row35 110 nodes, worst 2.29x, rule 9's col7 1.65x; QPB children sum (P = 1): col16 2.29x (-0.30); exp, L = 32: col16 2.29x (-0.05)
- tho40 26:3,31:7,37:0: best line row2 110 nodes, worst 3.28x, rule 9's col39 1.65x; QPB children sum (P = 1): col32 1.00x (+0.82); exp, L = 32: col32 1.00x (+0.94)

| children bounded by | score | mean Spearman with the truth | best line first | its first line / the best, mean (worst) |
|---|---|---|---|---|
| QPB | sum (P = 1) | +0.45 | 2/6 | 1.38x (2.29x) |
| QPB | power 8 | +0.57 | 4/6 | 1.25x (2.29x) |
| QPB | power 32 | +0.58 | 5/6 | 1.22x (2.29x) |
| QPB | exp, L = 16 | +0.60 | 5/6 | 1.22x (2.29x) |
| QPB | exp, L = 32 | +0.61 | 5/6 | 1.22x (2.29x) |
| QPB | largest shortfall | +0.47 | 2/6 | 1.82x (2.51x) |
| DNN | sum (P = 1) | +0.80 | 6/6 | 1.00x (1.00x) |
| DNN | power 8 | +0.62 | 5/6 | 1.14x (1.84x) |
| DNN | power 32 | +0.55 | 4/6 | 1.36x (2.31x) |
| DNN | exp, L = 16 | +0.80 | 6/6 | 1.00x (1.00x) |
| DNN | exp, L = 32 | +0.68 | 5/6 | 1.14x (1.84x) |
| DNN | largest shortfall | +0.52 | 4/6 | 1.36x (2.31x) |

Child by child (1584 children): Spearman of the subtree's nodes with the child's QPB shortfall +0.06, with its DNN first solve's +0.90; QPB's with the DNN's -0.07.
