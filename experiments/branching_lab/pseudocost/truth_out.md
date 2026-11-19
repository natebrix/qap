# Step 3: ranking the candidate lines at the line-truth nodes (truth.py)

## tho30 5:1

training: tho30-dnn, 12,687 children; node first solve 145112.0, cutoff 149,934.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col9 | 1,850 | 1 | 1 | 1 | 2 | 1 | 3,424 |
| col29 | 1,911 | 2 | 2 | 2 | 3 | 3 | 3,951 |
| col19 | 2,206 | 3 | 3 | 3 | 4 | 5 | 4,231 |
| row18 | 2,210 | 4 | 4 | 6 | 8 | 6 | 7,676 |
| col0 | 2,286 | 5 | 5 | 4 | 1 | 2 | 8,503 |
| row17 | 2,391 | 6 | 6 | 7 | 6 | 4 | 9,050 |
| row28 | 2,536 | 7 | 7 | 5 | 7 | 8 | 9,116 |
| row25 | 4,789 | 8 | 8 | 8 | 5 | 7 | 15,314 |

Spearman with the truth: gbt +1.00, pc +0.88, QPB lookahead (rule 4) +0.48, (a) the children's first DNN solve +0.71

Child level, the node's 232 candidate children against their (a) solves (132 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.78 / 0.92; pc +0.42 / 0.74; -(bound + U_ij) +0.72 / 0.79; x_ij (marginal; larger = more likely open) +0.67 / 0.73

## tho30 5:2

training: tho30-dnn, 13,674 children; node first solve 145062.2, cutoff 149,934.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col9 | 864 | 1 | 3 | 1 | 2 | 1 | 5,938 |
| col29 | 894 | 2 | 1 | 2 | 3 | 3 | 4,857 |
| row18 | 1,004 | 3 | 2 | 6 | 8 | 5 | 5,247 |
| row28 | 1,223 | 4 | 6 | 5 | 6 | 6 | 11,264 |
| col0 | 1,267 | 5 | 5 | 3 | 1 | 2 | 7,502 |
| row17 | 1,352 | 6 | 4 | 7 | 7 | 4 | 7,152 |
| row19 | 1,993 | 7 | 7 | 4 | 5 | 8 | 12,765 |
| row25 | 3,091 | 8 | 8 | 8 | 4 | 7 | 17,330 |

Spearman with the truth: gbt +0.83, pc +0.71, QPB lookahead (rule 4) +0.19, (a) the children's first DNN solve +0.71

Child level, the node's 232 candidate children against their (a) solves (134 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.84 / 0.94; pc +0.44 / 0.71; -(bound + U_ij) +0.78 / 0.81; x_ij (marginal; larger = more likely open) +0.71 / 0.78

## tho30 5:3

training: tho30-dnn, 14,422 children; node first solve 145345.4, cutoff 149,934.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col29 | 385 | 1.5 | 4 | 2 | 3 | 3 | 11,670 |
| col9 | 385 | 1.5 | 6 | 1 | 2 | 1 | 12,180 |
| col0 | 493 | 3 | 5 | 4 | 1 | 2 | 11,829 |
| row18 | 496 | 4 | 1 | 6 | 7 | 4 | 3,775 |
| row28 | 767 | 5 | 3 | 5 | 5 | 5 | 9,567 |
| row19 | 966 | 6 | 7 | 3 | 6 | 7 | 16,715 |
| row25 | 1,812 | 7 | 8 | 8 | 4 | 6 | 25,659 |
| col14 | 3,084 | 8 | 2 | 7 | 8 | 8 | 8,204 |

Spearman with the truth: gbt +0.02, pc +0.80, QPB lookahead (rule 4) +0.71, (a) the children's first DNN solve +0.93

Child level, the node's 232 candidate children against their (a) solves (126 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.70 / 0.84; pc +0.44 / 0.74; -(bound + U_ij) +0.78 / 0.84; x_ij (marginal; larger = more likely open) +0.76 / 0.82

## tho40 24:7,26:1,39:0

training: tho40ab-32, 22,649 children; node first solve 237670.1, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col32 | 253 | 1 | 3 | 8 | 1 | 1 | 1,813 |
| col39 | 359 | 2 | 5 | 2 | 2 | 3 | 3,593 |
| row35 | 431 | 3 | 7 | 5 | 5 | 6 | 12,640 |
| col15 | 432 | 4 | 4 | 1 | 3 | 4 | 2,762 |
| row2 | 466 | 5 | 8 | 7 | 4 | 2 | 40,163 |
| row18 | 578 | 6 | 6 | 6 | 8 | 5 | 5,670 |
| row14 | 856 | 7 | 1 | 3 | 6 | 7 | 0 |
| row7 | 1,175 | 8 | 2 | 4 | 7 | 8 | 1,448 |

Spearman with the truth: gbt -0.31, pc -0.17, QPB lookahead (rule 4) +0.86, (a) the children's first DNN solve +0.76

Child level, the node's 296 candidate children against their (a) solves (81 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.35 / 0.65; pc +0.17 / 0.53; -(bound + U_ij) +0.64 / 0.92; x_ij (marginal; larger = more likely open) +0.71 / 0.93

training: tho40-dnn less the truth subtrees, 28,886 children; node first solve 237670.1, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col32 | 253 | 1 | 3 | 7 | 1 | 1 | 2,814 |
| col39 | 359 | 2 | 4 | 3 | 2 | 3 | 4,002 |
| row35 | 431 | 3 | 8 | 6 | 5 | 6 | 15,277 |
| col15 | 432 | 4 | 6 | 5 | 3 | 4 | 5,684 |
| row2 | 466 | 5 | 7 | 8 | 4 | 2 | 13,057 |
| row18 | 578 | 6 | 5 | 2 | 8 | 5 | 5,040 |
| row14 | 856 | 7 | 1 | 1 | 6 | 7 | 67 |
| row7 | 1,175 | 8 | 2 | 4 | 7 | 8 | 496 |

Spearman with the truth: gbt -0.36, pc -0.48, QPB lookahead (rule 4) +0.86, (a) the children's first DNN solve +0.76

Child level, the node's 296 candidate children against their (a) solves (81 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.46 / 0.70; pc +0.35 / 0.67; -(bound + U_ij) +0.64 / 0.92; x_ij (marginal; larger = more likely open) +0.71 / 0.93

training: both, 51,535 children; node first solve 237670.1, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| col32 | 253 | 1 | 3 | 8 | 1 | 1 | 2,688 |
| col39 | 359 | 2 | 5 | 5 | 2 | 3 | 4,249 |
| row35 | 431 | 3 | 8 | 6 | 5 | 6 | 14,569 |
| col15 | 432 | 4 | 4 | 2 | 3 | 4 | 4,128 |
| row2 | 466 | 5 | 7 | 7 | 4 | 2 | 8,940 |
| row18 | 578 | 6 | 6 | 3 | 8 | 5 | 4,971 |
| row14 | 856 | 7 | 1 | 1 | 6 | 7 | 122 |
| row7 | 1,175 | 8 | 2 | 4 | 7 | 8 | 713 |

Spearman with the truth: gbt -0.36, pc -0.62, QPB lookahead (rule 4) +0.86, (a) the children's first DNN solve +0.76

Child level, the node's 296 candidate children against their (a) solves (81 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.49 / 0.72; pc +0.36 / 0.69; -(bound + U_ij) +0.64 / 0.92; x_ij (marginal; larger = more likely open) +0.71 / 0.93

## tho40 25:32,26:8,39:0

training: tho40ab-32, 22,649 children; node first solve 238042.0, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row35 | 110 | 1 | 7 | 4 | 6 | 1 | 13,187 |
| row2 | 146 | 2 | 8 | 7 | 5 | 2 | 20,739 |
| row30 | 180 | 3 | 2 | 3 | 7 | 3 | 774 |
| col7 | 182 | 4.5 | 3 | 8 | 2 | 6 | 971 |
| col39 | 182 | 4.5 | 4 | 6 | 4 | 7 | 2,079 |
| col24 | 216 | 6 | 5 | 5 | 3 | 4 | 2,158 |
| row37 | 252 | 7.5 | 1 | 1 | 8 | 8 | 146 |
| col16 | 252 | 7.5 | 6 | 2 | 1 | 5 | 2,261 |

Spearman with the truth: gbt -0.46, pc -0.47, QPB lookahead (rule 4) -0.30, (a) the children's first DNN solve +0.77

Child level, the node's 296 candidate children against their (a) solves (26 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.19 / 0.64; pc +0.11 / 0.60; -(bound + U_ij) +0.40 / 0.90; x_ij (marginal; larger = more likely open) +0.44 / 0.95

training: tho40-dnn less the truth subtrees, 29,063 children; node first solve 238042.0, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row35 | 110 | 1 | 6 | 6 | 6 | 1 | 4,905 |
| row2 | 146 | 2 | 7 | 7 | 5 | 2 | 5,346 |
| row30 | 180 | 3 | 8 | 5 | 7 | 3 | 5,703 |
| col7 | 182 | 4.5 | 1 | 8 | 2 | 6 | 1,104 |
| col39 | 182 | 4.5 | 2 | 2 | 4 | 7 | 1,545 |
| col24 | 216 | 6 | 3 | 3 | 3 | 4 | 2,630 |
| row37 | 252 | 7.5 | 5 | 4 | 8 | 8 | 4,515 |
| col16 | 252 | 7.5 | 4 | 1 | 1 | 5 | 3,491 |

Spearman with the truth: gbt -0.46, pc -0.64, QPB lookahead (rule 4) -0.30, (a) the children's first DNN solve +0.77

Child level, the node's 296 candidate children against their (a) solves (26 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.41 / 0.93; pc +0.40 / 0.87; -(bound + U_ij) +0.40 / 0.90; x_ij (marginal; larger = more likely open) +0.44 / 0.95

training: both, 51,712 children; node first solve 238042.0, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row35 | 110 | 1 | 7 | 6 | 6 | 1 | 5,761 |
| row2 | 146 | 2 | 8 | 7 | 5 | 2 | 6,556 |
| row30 | 180 | 3 | 6 | 5 | 7 | 3 | 4,350 |
| col7 | 182 | 4.5 | 1 | 8 | 2 | 6 | 878 |
| col39 | 182 | 4.5 | 4 | 4 | 4 | 7 | 2,057 |
| col24 | 216 | 6 | 3 | 3 | 3 | 4 | 1,814 |
| row37 | 252 | 7.5 | 5 | 2 | 8 | 8 | 3,095 |
| col16 | 252 | 7.5 | 2 | 1 | 1 | 5 | 1,681 |

Spearman with the truth: gbt -0.67, pc -0.78, QPB lookahead (rule 4) -0.30, (a) the children's first DNN solve +0.77

Child level, the node's 296 candidate children against their (a) solves (26 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.39 / 0.88; pc +0.32 / 0.81; -(bound + U_ij) +0.40 / 0.90; x_ij (marginal; larger = more likely open) +0.44 / 0.95

## tho40 26:3,31:7,37:0

training: tho40ab-32, 22,649 children; node first solve 238250.3, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row2 | 110 | 1.5 | 8 | 7 | 4 | 1 | 21,797 |
| col32 | 110 | 1.5 | 4 | 8 | 1 | 2 | 553 |
| col39 | 182 | 3.5 | 6 | 3 | 2 | 3 | 1,259 |
| col1 | 182 | 3.5 | 5 | 1 | 3 | 4 | 1,135 |
| col2 | 253 | 5 | 7 | 2 | 7 | 5 | 1,360 |
| row30 | 254 | 6 | 3 | 6 | 5 | 7 | 550 |
| row14 | 290 | 7 | 1 | 4 | 6 | 8 | 0 |
| row25 | 361 | 8 | 2 | 5 | 8 | 6 | 20 |

Spearman with the truth: gbt -0.71, pc -0.28, QPB lookahead (rule 4) +0.82, (a) the children's first DNN solve +0.92

Child level, the node's 296 candidate children against their (a) solves (38 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.29 / 0.61; pc +0.16 / 0.46; -(bound + U_ij) +0.50 / 0.94; x_ij (marginal; larger = more likely open) +0.52 / 0.95

training: tho40-dnn less the truth subtrees, 29,063 children; node first solve 238250.3, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row2 | 110 | 1.5 | 8 | 8 | 4 | 1 | 5,229 |
| col32 | 110 | 1.5 | 5 | 7 | 1 | 2 | 1,115 |
| col39 | 182 | 3.5 | 4 | 5 | 2 | 3 | 1,022 |
| col1 | 182 | 3.5 | 6 | 4 | 3 | 4 | 1,349 |
| col2 | 253 | 5 | 3 | 2 | 7 | 5 | 748 |
| row30 | 254 | 6 | 7 | 6 | 5 | 7 | 1,928 |
| row14 | 290 | 7 | 2 | 1 | 6 | 8 | 120 |
| row25 | 361 | 8 | 1 | 3 | 8 | 6 | 66 |

Spearman with the truth: gbt -0.69, pc -0.75, QPB lookahead (rule 4) +0.82, (a) the children's first DNN solve +0.92

Child level, the node's 296 candidate children against their (a) solves (38 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.44 / 0.80; pc +0.25 / 0.65; -(bound + U_ij) +0.50 / 0.94; x_ij (marginal; larger = more likely open) +0.52 / 0.95

training: both, 51,712 children; node first solve 238250.3, cutoff 240,514.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | (a) the children's first DNN solve rank | gbt shortfall sum |
|---|---|---|---|---|---|---|---|
| row2 | 110 | 1.5 | 8 | 7 | 4 | 1 | 5,005 |
| col32 | 110 | 1.5 | 3 | 8 | 1 | 2 | 1,126 |
| col39 | 182 | 3.5 | 4 | 5 | 2 | 3 | 1,146 |
| col1 | 182 | 3.5 | 6 | 4 | 3 | 4 | 1,682 |
| col2 | 253 | 5 | 5 | 3 | 7 | 5 | 1,278 |
| row30 | 254 | 6 | 7 | 6 | 5 | 7 | 2,202 |
| row14 | 290 | 7 | 1 | 1 | 6 | 8 | 59 |
| row25 | 361 | 8 | 2 | 2 | 8 | 6 | 222 |

Spearman with the truth: gbt -0.49, pc -0.82, QPB lookahead (rule 4) +0.82, (a) the children's first DNN solve +0.92

Child level, the node's 296 candidate children against their (a) solves (38 open): Spearman within a line / AUC open vs fathomed over the node: gbt +0.40 / 0.75; pc +0.25 / 0.63; -(bound + U_ij) +0.50 / 0.94; x_ij (marginal; larger = more likely open) +0.52 / 0.95

## tai60b 44:45,43:58

training: tai60b-261117, 175,721 children; node first solve 606420321.9, cutoff 608,215,054.05

| line | truth nodes | rank | gbt rank | pc rank | QPB lookahead (rule 4) rank | gbt shortfall sum |
|---|---|---|---|---|---|---|
| row1 | 19,600 | 1 | 4 | 4 | 2 | 4,214,802 |
| row45 | 22,900 | 2 | 2 | 2 | 4 | 1,352,240 |
| col26 | 38,800 | 3 | 1 | 1 | 3 | 837,118 |
| row8 | 67,100 | 4 | 3 | 3 | 1 | 3,618,530 |

Spearman with the truth: gbt -0.40, pc -0.40, QPB lookahead (rule 4) -0.40

## Over the nodes

| training | score | nodes | mean Spearman | best line 1st | its top line / best, mean (worst) |
|---|---|---|---|---|---|
| tho30 | gbt: predicted shortfall sum | 3 | +0.62 | 1/3 | 1.11x (1.29x) |
| tho30 | pc: predicted shortfall sum | 3 | +0.80 | 2/3 | 1.00x (1.00x) |
| tho30 | QPB lookahead (rule 4) | 3 | +0.46 | 0/3 | 1.33x (1.47x) |
| tho30 | (a) the children's first DNN solve | 3 | +0.79 | 2/3 | 1.00x (1.00x) |
| tho40 [tho40ab-32] | gbt: predicted shortfall sum | 3 | -0.49 | 0/3 | 2.77x (3.38x) |
| tho40 [tho40ab-32] | pc: predicted shortfall sum | 3 | -0.30 | 0/3 | 1.88x (2.29x) |
| tho40 [tho40ab-32] | QPB lookahead (rule 4) | 3 | +0.46 | 1/3 | 1.43x (2.29x) |
| tho40 [tho40ab-32] | (a) the children's first DNN solve | 3 | +0.82 | 3/3 | 1.00x (1.00x) |
| tho40 [tho40-dnn less the truth subtrees] | gbt: predicted shortfall sum | 3 | -0.50 | 0/3 | 2.77x (3.38x) |
| tho40 [tho40-dnn less the truth subtrees] | pc: predicted shortfall sum | 3 | -0.62 | 0/3 | 2.77x (3.38x) |
| tho40 [tho40-dnn less the truth subtrees] | QPB lookahead (rule 4) | 3 | +0.46 | 1/3 | 1.43x (2.29x) |
| tho40 [tho40-dnn less the truth subtrees] | (a) the children's first DNN solve | 3 | +0.82 | 3/3 | 1.00x (1.00x) |
| tho40 [both] | gbt: predicted shortfall sum | 3 | -0.51 | 0/3 | 2.56x (3.38x) |
| tho40 [both] | pc: predicted shortfall sum | 3 | -0.74 | 0/3 | 2.77x (3.38x) |
| tho40 [both] | QPB lookahead (rule 4) | 3 | +0.46 | 1/3 | 1.43x (2.29x) |
| tho40 [both] | (a) the children's first DNN solve | 3 | +0.82 | 3/3 | 1.00x (1.00x) |
| tai60b | gbt: predicted shortfall sum | 1 | -0.40 | 0/1 | 1.98x (1.98x) |
| tai60b | pc: predicted shortfall sum | 1 | -0.40 | 0/1 | 1.98x (1.98x) |
| tai60b | QPB lookahead (rule 4) | 1 | -0.40 | 0/1 | 3.42x (3.42x) |

## Child level over the nodes (Spearman within a line / AUC over the node)

| training | score | nodes | Spearman within a line | AUC open |
|---|---|---|---|---|
| tho30 | gbt | 3 | +0.77 | 0.90 |
| tho30 | pc | 3 | +0.43 | 0.73 |
| tho30 | -(bound + U_ij) | 3 | +0.76 | 0.82 |
| tho30 | x_ij (marginal; larger = more likely open) | 3 | +0.71 | 0.77 |
| tho40 [tho40ab-32] | gbt | 3 | +0.28 | 0.64 |
| tho40 [tho40ab-32] | pc | 3 | +0.15 | 0.53 |
| tho40 [tho40ab-32] | -(bound + U_ij) | 3 | +0.51 | 0.92 |
| tho40 [tho40ab-32] | x_ij (marginal; larger = more likely open) | 3 | +0.56 | 0.95 |
| tho40 [tho40-dnn less the truth subtrees] | gbt | 3 | +0.44 | 0.81 |
| tho40 [tho40-dnn less the truth subtrees] | pc | 3 | +0.33 | 0.73 |
| tho40 [tho40-dnn less the truth subtrees] | -(bound + U_ij) | 3 | +0.51 | 0.92 |
| tho40 [tho40-dnn less the truth subtrees] | x_ij (marginal; larger = more likely open) | 3 | +0.56 | 0.95 |
| tho40 [both] | gbt | 3 | +0.42 | 0.78 |
| tho40 [both] | pc | 3 | +0.31 | 0.71 |
| tho40 [both] | -(bound + U_ij) | 3 | +0.51 | 0.92 |
| tho40 [both] | x_ij (marginal; larger = more likely open) | 3 | +0.56 | 0.95 |
