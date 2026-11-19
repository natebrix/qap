# Step 2: out-of-sample prediction of a child's first solve (predict.py)

## tai60b-261117

### time: train on the first 25% of the trace

train 82,787 children; test 228,279 children of 5,472 parents (4,916 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.50 | -0.00 | 0.07 | 0.26 | 0.07 |
| qpb | 0.97 | 0.39 | 0.66 | 0.99 | 0.25 |
| pc | 0.90 | 0.32 | 0.65 | 0.83 | 0.21 |
| lin | 0.94 | 0.36 | 0.57 | 0.87 | 0.22 |
| lin+q | 0.97 | 0.39 | 0.69 | 0.96 | 0.24 |
| gbt | 0.92 | 0.33 | 0.72 | 0.87 | 0.22 |
| gbt+q | 0.92 | 0.33 | 0.74 | 0.87 | 0.22 |

### time: train on the first 50% of the trace

train 165,573 children; test 150,512 children of 3,845 parents (3,402 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.50 | 0.00 | 0.08 | 0.26 | 0.08 |
| qpb | 0.97 | 0.41 | 0.71 | 0.98 | 0.28 |
| pc | 0.86 | 0.31 | 0.65 | 0.82 | 0.24 |
| lin | 0.93 | 0.37 | 0.59 | 0.86 | 0.25 |
| lin+q | 0.98 | 0.42 | 0.76 | 0.99 | 0.29 |
| gbt | 0.93 | 0.38 | 0.73 | 0.90 | 0.26 |
| gbt+q | 0.94 | 0.38 | 0.73 | 0.90 | 0.26 |

### depth: train on depths <= 11, test below

train 108,002 children; test 223,145 children of 5,379 parents (4,820 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.50 | 0.00 | 0.07 | 0.26 | 0.07 |
| qpb | 0.97 | 0.39 | 0.68 | 0.99 | 0.25 |
| pc | 0.88 | 0.31 | 0.53 | 0.83 | 0.21 |
| lin | 0.94 | 0.36 | 0.60 | 0.90 | 0.22 |
| lin+q | 0.97 | 0.39 | 0.69 | 0.96 | 0.24 |
| gbt | 0.88 | 0.31 | 0.61 | 0.85 | 0.21 |
| gbt+q | 0.92 | 0.34 | 0.64 | 0.89 | 0.22 |

## tai50b-261116
