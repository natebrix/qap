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
| gbt-w | 0.93 | 0.35 | 0.75 | 0.90 | 0.23 |
| gbt-w+q | 0.93 | 0.35 | 0.74 | 0.91 | 0.23 |

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
| gbt-w | 0.95 | 0.39 | 0.73 | 0.94 | 0.27 |
| gbt-w+q | 0.93 | 0.38 | 0.69 | 0.90 | 0.26 |

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
| gbt-w | 0.85 | 0.28 | 0.59 | 0.80 | 0.20 |
| gbt-w+q | 0.86 | 0.29 | 0.62 | 0.81 | 0.20 |

## tai50b-261116

### time: train on the first 25% of the trace

train 2,054 children; test 6,080 children of 160 parents (154 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.00 | 0.03 | 0.24 | 0.05 |
| qpb | 0.91 | 0.31 | 0.40 | 0.82 | 0.17 |
| pc | 0.55 | 0.05 | 0.11 | 0.38 | 0.08 |
| lin | 0.92 | 0.31 | 0.69 | 0.88 | 0.19 |
| lin+q | 0.92 | 0.31 | 0.69 | 0.88 | 0.19 |
| gbt | 0.84 | 0.25 | 0.60 | 0.59 | 0.12 |
| gbt+q | 0.80 | 0.21 | 0.59 | 0.59 | 0.12 |
| gbt-w | 0.80 | 0.21 | 0.60 | 0.58 | 0.12 |
| gbt-w+q | 0.89 | 0.29 | 0.65 | 0.76 | 0.16 |

### time: train on the first 50% of the trace

train 4,107 children; test 4,070 children of 114 parents (109 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.48 | -0.02 | 0.03 | 0.25 | 0.05 |
| qpb | 0.96 | 0.34 | 0.48 | 0.97 | 0.19 |
| pc | 0.63 | 0.12 | 0.29 | 0.61 | 0.12 |
| lin | 0.93 | 0.31 | 0.74 | 0.92 | 0.18 |
| lin+q | 0.93 | 0.31 | 0.74 | 0.92 | 0.18 |
| gbt | 0.96 | 0.34 | 0.80 | 0.94 | 0.18 |
| gbt+q | 0.93 | 0.31 | 0.75 | 0.83 | 0.16 |
| gbt-w | 0.95 | 0.32 | 0.77 | 0.91 | 0.18 |
| gbt-w+q | 0.95 | 0.33 | 0.76 | 0.90 | 0.18 |

### depth: train on depths <= 6, test below

train 2,810 children; test 5,404 children of 145 parents (140 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.55 | 0.04 | 0.06 | 0.29 | 0.06 |
| qpb | 0.93 | 0.32 | 0.42 | 0.87 | 0.18 |
| pc | 0.55 | 0.05 | 0.07 | 0.38 | 0.08 |
| lin | 0.94 | 0.33 | 0.72 | 0.94 | 0.19 |
| lin+q | 0.94 | 0.33 | 0.73 | 0.93 | 0.19 |
| gbt | 0.81 | 0.21 | 0.64 | 0.61 | 0.13 |
| gbt+q | 0.80 | 0.20 | 0.61 | 0.58 | 0.12 |
| gbt-w | 0.83 | 0.23 | 0.62 | 0.63 | 0.13 |
| gbt-w+q | 0.81 | 0.22 | 0.64 | 0.63 | 0.13 |

## tai50b-r9

### time: train on the first 25% of the trace

train 1,500 children; test 4,234 children of 106 parents (104 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.48 | -0.01 | 0.04 | 0.26 | 0.05 |
| qpb | 0.86 | 0.25 | 0.29 | 0.68 | 0.13 |
| pc | 0.60 | 0.06 | 0.12 | 0.35 | 0.06 |
| lin | 0.85 | 0.23 | 0.56 | 0.72 | 0.13 |
| lin+q | 0.84 | 0.23 | 0.55 | 0.69 | 0.13 |
| gbt | 0.65 | 0.10 | 0.13 | 0.47 | 0.09 |
| gbt+q | 0.54 | 0.03 | 0.06 | 0.34 | 0.06 |
| gbt-w | 0.57 | 0.05 | 0.13 | 0.42 | 0.08 |
| gbt-w+q | 0.66 | 0.10 | 0.10 | 0.41 | 0.08 |

### time: train on the first 50% of the trace

train 3,000 children; test 2,927 children of 77 parents (76 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.52 | 0.01 | 0.06 | 0.26 | 0.04 |
| qpb | 0.90 | 0.27 | 0.32 | 0.80 | 0.14 |
| pc | 0.71 | 0.13 | 0.26 | 0.45 | 0.08 |
| lin | 0.85 | 0.23 | 0.60 | 0.73 | 0.12 |
| lin+q | 0.86 | 0.24 | 0.60 | 0.75 | 0.13 |
| gbt | 0.81 | 0.20 | 0.62 | 0.71 | 0.12 |
| gbt+q | 0.82 | 0.21 | 0.63 | 0.67 | 0.11 |
| gbt-w | 0.79 | 0.19 | 0.59 | 0.67 | 0.11 |
| gbt-w+q | 0.85 | 0.23 | 0.58 | 0.72 | 0.12 |

### depth: train on depths <= 5, test below

train 2,180 children; test 3,820 children of 97 parents (95 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.53 | 0.01 | 0.05 | 0.30 | 0.06 |
| qpb | 0.86 | 0.25 | 0.30 | 0.67 | 0.13 |
| pc | 0.70 | 0.15 | 0.29 | 0.49 | 0.09 |
| lin | 0.83 | 0.23 | 0.56 | 0.67 | 0.13 |
| lin+q | 0.83 | 0.23 | 0.56 | 0.67 | 0.13 |
| gbt | 0.80 | 0.20 | 0.55 | 0.62 | 0.12 |
| gbt+q | 0.78 | 0.19 | 0.54 | 0.61 | 0.12 |
| gbt-w | 0.82 | 0.22 | 0.53 | 0.68 | 0.13 |
| gbt-w+q | 0.88 | 0.26 | 0.47 | 0.73 | 0.14 |

## tho40ab-32

### time: train on the first 25% of the trace

train 5,662 children; test 16,300 children of 467 parents (142 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.51 | 0.01 | 0.05 | 0.26 | 0.07 |
| qpb | 0.98 | 0.39 | 0.80 | 0.97 | 0.25 |
| pc | 0.88 | 0.31 | 0.44 | 0.81 | 0.21 |
| lin | 0.91 | 0.32 | 0.35 | 0.87 | 0.22 |
| lin+q | 0.90 | 0.32 | 0.32 | 0.87 | 0.22 |
| gbt | 0.81 | 0.26 | 0.45 | 0.81 | 0.21 |
| gbt+q | 0.84 | 0.28 | 0.56 | 0.85 | 0.22 |
| gbt-w | 0.76 | 0.23 | 0.47 | 0.75 | 0.19 |
| gbt-w+q | 0.90 | 0.33 | 0.67 | 0.92 | 0.24 |

### time: train on the first 50% of the trace

train 11,324 children; test 10,789 children of 311 parents (93 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.51 | 0.00 | 0.08 | 0.26 | 0.06 |
| qpb | 0.99 | 0.38 | 0.82 | 0.99 | 0.24 |
| pc | 0.97 | 0.37 | 0.65 | 0.95 | 0.24 |
| lin | 0.89 | 0.31 | 0.31 | 0.85 | 0.21 |
| lin+q | 0.92 | 0.33 | 0.47 | 0.85 | 0.21 |
| gbt | 0.96 | 0.37 | 0.67 | 0.96 | 0.24 |
| gbt+q | 0.96 | 0.36 | 0.76 | 0.94 | 0.23 |
| gbt-w | 0.96 | 0.36 | 0.69 | 0.94 | 0.23 |
| gbt-w+q | 0.99 | 0.39 | 0.86 | 0.99 | 0.24 |

### depth: train on depths <= 5, test below

train 8,997 children; test 13,652 children of 395 parents (78 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.44 | -0.04 | 0.01 | 0.17 | 0.03 |
| qpb | 1.00 | 0.36 | 0.91 | 1.00 | 0.20 |
| pc | 0.95 | 0.33 | 0.53 | 0.91 | 0.18 |
| lin | 0.85 | 0.25 | 0.17 | 0.81 | 0.16 |
| lin+q | 0.89 | 0.28 | 0.38 | 0.81 | 0.16 |
| gbt | 0.83 | 0.25 | 0.57 | 0.86 | 0.17 |
| gbt+q | 0.88 | 0.28 | 0.64 | 0.88 | 0.18 |
| gbt-w | 0.89 | 0.29 | 0.59 | 0.91 | 0.18 |
| gbt-w+q | 0.74 | 0.20 | 0.61 | 0.83 | 0.17 |

## tho40ab-0

### time: train on the first 25% of the trace

train 5,610 children; test 15,530 children of 445 parents (125 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.01 | 0.07 | 0.24 | 0.06 |
| qpb | 0.99 | 0.39 | 0.83 | 0.96 | 0.25 |
| pc | 0.96 | 0.37 | 0.54 | 0.92 | 0.24 |
| lin | 0.89 | 0.31 | 0.30 | 0.86 | 0.23 |
| lin+q | 0.93 | 0.34 | 0.42 | 0.86 | 0.23 |
| gbt | 0.92 | 0.34 | 0.66 | 0.91 | 0.24 |
| gbt+q | 0.89 | 0.33 | 0.69 | 0.90 | 0.24 |
| gbt-w | 0.92 | 0.34 | 0.56 | 0.92 | 0.24 |
| gbt-w+q | 0.94 | 0.36 | 0.69 | 0.93 | 0.24 |

### time: train on the first 50% of the trace

train 11,219 children; test 9,750 children of 282 parents (82 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.51 | 0.01 | 0.07 | 0.29 | 0.07 |
| qpb | 0.99 | 0.38 | 0.87 | 0.98 | 0.24 |
| pc | 0.97 | 0.37 | 0.59 | 0.94 | 0.23 |
| lin | 0.87 | 0.29 | 0.23 | 0.83 | 0.20 |
| lin+q | 0.91 | 0.32 | 0.43 | 0.83 | 0.20 |
| gbt | 0.98 | 0.37 | 0.70 | 0.98 | 0.24 |
| gbt+q | 0.96 | 0.36 | 0.79 | 0.96 | 0.23 |
| gbt-w | 0.98 | 0.37 | 0.72 | 0.97 | 0.23 |
| gbt-w+q | 0.99 | 0.38 | 0.88 | 0.97 | 0.23 |

### depth: train on depths <= 5, test below

train 8,925 children; test 13,513 children of 391 parents (79 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.48 | -0.01 | 0.04 | 0.26 | 0.05 |
| qpb | 1.00 | 0.36 | 0.92 | 1.00 | 0.19 |
| pc | 0.95 | 0.33 | 0.50 | 0.98 | 0.19 |
| lin | 0.86 | 0.25 | 0.16 | 0.82 | 0.16 |
| lin+q | 0.91 | 0.29 | 0.42 | 0.82 | 0.16 |
| gbt | 0.83 | 0.25 | 0.63 | 0.83 | 0.16 |
| gbt+q | 0.74 | 0.18 | 0.46 | 0.72 | 0.14 |
| gbt-w | 0.81 | 0.24 | 0.60 | 0.86 | 0.17 |
| gbt-w+q | 0.92 | 0.31 | 0.64 | 0.93 | 0.18 |

## tho40-dnn

### time: train on the first 25% of the trace

train 7,311 children; test 21,019 children of 672 parents (234 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.46 | -0.03 | 0.14 | 0.25 | 0.08 |
| pc | 0.87 | 0.35 | 0.45 | 0.75 | 0.23 |
| lin | 0.92 | 0.42 | 0.51 | 0.81 | 0.25 |
| gbt | 0.93 | 0.42 | 0.70 | 0.82 | 0.25 |
| gbt-w | 0.92 | 0.43 | 0.66 | 0.82 | 0.25 |

### time: train on the first 50% of the trace

train 14,622 children; test 13,552 children of 438 parents (153 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.00 | 0.16 | 0.27 | 0.08 |
| pc | 0.88 | 0.36 | 0.50 | 0.76 | 0.22 |
| lin | 0.94 | 0.45 | 0.54 | 0.85 | 0.24 |
| gbt | 0.94 | 0.44 | 0.73 | 0.83 | 0.24 |
| gbt-w | 0.93 | 0.44 | 0.71 | 0.83 | 0.24 |

### depth: train on depths <= 5, test below

train 13,786 children; test 15,459 children of 472 parents (106 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.46 | -0.01 | 0.08 | 0.21 | 0.04 |
| pc | 0.91 | 0.29 | 0.42 | 0.86 | 0.15 |
| lin | 0.95 | 0.33 | 0.49 | 0.96 | 0.17 |
| gbt | 0.92 | 0.30 | 0.61 | 0.88 | 0.15 |
| gbt-w | 0.93 | 0.32 | 0.61 | 0.92 | 0.16 |

## tho30-dnn

### time: train on the first 25% of the trace

train 3,696 children; test 10,749 children of 436 parents (130 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.51 | 0.01 | 0.13 | 0.27 | 0.11 |
| pc | 0.85 | 0.34 | 0.62 | 0.72 | 0.29 |
| lin | 0.90 | 0.39 | 0.52 | 0.73 | 0.29 |
| gbt | 0.86 | 0.36 | 0.64 | 0.78 | 0.31 |
| gbt-w | 0.90 | 0.40 | 0.66 | 0.79 | 0.32 |

### time: train on the first 50% of the trace

train 7,391 children; test 7,049 children of 282 parents (81 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.02 | 0.10 | 0.20 | 0.08 |
| pc | 0.94 | 0.42 | 0.63 | 0.81 | 0.33 |
| lin | 0.92 | 0.40 | 0.58 | 0.73 | 0.30 |
| gbt | 0.94 | 0.43 | 0.67 | 0.82 | 0.33 |
| gbt-w | 0.96 | 0.44 | 0.73 | 0.87 | 0.35 |

### depth: train on depths <= 4, test below

train 9,973 children; test 4,810 children of 203 parents (34 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.47 | -0.02 | 0.04 | 0.25 | 0.05 |
| pc | 0.92 | 0.32 | 0.62 | 0.90 | 0.17 |
| lin | 0.84 | 0.25 | 0.38 | 0.68 | 0.13 |
| gbt | 0.84 | 0.26 | 0.53 | 0.80 | 0.15 |
| gbt-w | 0.91 | 0.32 | 0.56 | 0.90 | 0.17 |

## nug30-trace-a

### time: train on the first 25% of the trace

train 922 children; test 2,605 children of 98 parents (30 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.45 | -0.03 | 0.11 | 0.24 | 0.09 |
| qpb | 0.99 | 0.46 | 0.84 | 0.86 | 0.34 |
| pc | 0.80 | 0.27 | 0.49 | 0.59 | 0.23 |
| lin | 0.92 | 0.40 | 0.39 | 0.78 | 0.30 |
| lin+q | 0.96 | 0.44 | 0.52 | 0.88 | 0.34 |
| gbt | 0.88 | 0.33 | 0.51 | 0.66 | 0.26 |
| gbt+q | 0.99 | 0.46 | 0.86 | 0.89 | 0.35 |
| gbt-w | 0.89 | 0.33 | 0.57 | 0.70 | 0.27 |
| gbt-w+q | 0.99 | 0.46 | 0.94 | 0.89 | 0.35 |

### time: train on the first 50% of the trace

train 1,843 children; test 1,581 children of 58 parents (17 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.47 | -0.03 | 0.04 | 0.22 | 0.05 |
| qpb | 1.00 | 0.38 | 0.90 | 0.96 | 0.22 |
| pc | 0.88 | 0.28 | 0.66 | 0.70 | 0.16 |
| lin | 0.93 | 0.33 | 0.37 | 0.93 | 0.21 |
| lin+q | 0.97 | 0.36 | 0.50 | 1.00 | 0.23 |
| gbt | 0.86 | 0.29 | 0.54 | 0.78 | 0.18 |
| gbt+q | 0.99 | 0.38 | 0.90 | 1.00 | 0.23 |
| gbt-w | 0.91 | 0.31 | 0.42 | 0.81 | 0.18 |
| gbt-w+q | 1.00 | 0.38 | 0.93 | 1.00 | 0.23 |

### depth: train on depths <= 3, test below

train 1,858 children; test 1,829 children of 70 parents (16 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.54 | 0.04 | 0.09 | 0.19 | 0.04 |
| qpb | 0.99 | 0.37 | 0.72 | 1.00 | 0.19 |
| pc | 0.78 | 0.19 | 0.34 | 0.57 | 0.11 |
| lin | 0.92 | 0.32 | 0.38 | 1.00 | 0.19 |
| lin+q | 0.94 | 0.33 | 0.44 | 0.95 | 0.18 |
| gbt | 0.82 | 0.23 | 0.41 | 0.71 | 0.14 |
| gbt+q | 0.87 | 0.28 | 0.59 | 0.81 | 0.16 |
| gbt-w | 0.85 | 0.25 | 0.59 | 0.76 | 0.15 |
| gbt-w+q | 0.92 | 0.31 | 0.72 | 0.90 | 0.17 |

## Across runs of one instance

### train tai50b-261116, test tai50b-r9

train 8,214 children; test 6,000 children of 143 parents (138 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.48 | -0.02 | 0.09 | 0.27 | 0.08 |
| qpb | 0.84 | 0.23 | 0.31 | 0.52 | 0.15 |
| pc | 0.90 | 0.32 | 0.66 | 0.78 | 0.22 |
| lin | 0.79 | 0.20 | 0.54 | 0.52 | 0.15 |
| lin+q | 0.85 | 0.25 | 0.57 | 0.60 | 0.17 |
| gbt | 0.95 | 0.35 | 0.80 | 0.85 | 0.24 |
| gbt+q | 0.94 | 0.34 | 0.79 | 0.84 | 0.24 |
| gbt-w | 0.95 | 0.35 | 0.82 | 0.84 | 0.24 |
| gbt-w+q | 0.95 | 0.35 | 0.82 | 0.85 | 0.24 |

### train tai50b-r9, test tai50b-261116

train 6,000 children; test 8,214 children of 205 parents (196 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.02 | 0.07 | 0.24 | 0.07 |
| qpb | 0.89 | 0.28 | 0.39 | 0.64 | 0.18 |
| pc | 0.89 | 0.32 | 0.71 | 0.82 | 0.23 |
| lin | 0.90 | 0.30 | 0.65 | 0.74 | 0.21 |
| lin+q | 0.90 | 0.31 | 0.65 | 0.75 | 0.21 |
| gbt | 0.97 | 0.37 | 0.88 | 0.88 | 0.24 |
| gbt+q | 0.98 | 0.37 | 0.84 | 0.88 | 0.24 |
| gbt-w | 0.97 | 0.37 | 0.86 | 0.89 | 0.25 |
| gbt-w+q | 0.98 | 0.38 | 0.88 | 0.90 | 0.25 |

### train tho40ab-0, test tho40ab-32

train 22,438 children; test 22,649 children of 644 parents (192 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.01 | 0.09 | 0.25 | 0.09 |
| qpb | 0.97 | 0.42 | 0.77 | 0.85 | 0.31 |
| pc | 0.99 | 0.43 | 0.76 | 0.87 | 0.32 |
| lin | 0.94 | 0.38 | 0.56 | 0.81 | 0.30 |
| lin+q | 0.96 | 0.41 | 0.63 | 0.81 | 0.30 |
| gbt | 0.97 | 0.41 | 0.86 | 0.84 | 0.31 |
| gbt+q | 0.99 | 0.43 | 0.95 | 0.86 | 0.32 |
| gbt-w | 0.98 | 0.42 | 0.86 | 0.84 | 0.31 |
| gbt-w+q | 0.99 | 0.43 | 0.95 | 0.86 | 0.31 |

### train tho40-dnn, test tho40ab-32

train 29,245 children; test 22,649 children of 644 parents (192 with an open child)

| score | AUC | Spearman | R-precision | recall, top quarter | precision, top quarter |
|---|---|---|---|---|---|
| random | 0.49 | -0.01 | 0.09 | 0.25 | 0.09 |
| pc | 0.92 | 0.37 | 0.45 | 0.76 | 0.28 |
| lin | 0.94 | 0.38 | 0.52 | 0.81 | 0.30 |
| gbt | 0.96 | 0.40 | 0.65 | 0.82 | 0.30 |
| gbt-w | 0.93 | 0.39 | 0.64 | 0.81 | 0.30 |

## Summary: AUC within parent (open vs fathomed), each score

| test | split | random | qpb | pc | lin | lin+q | gbt | gbt+q | gbt-w | gbt-w+q |
|---|---|---|---|---|---|---|---|---|---|---|
| tai60b-261117 | time 25% | 0.50 | 0.97 | 0.90 | 0.94 | 0.97 | 0.92 | 0.92 | 0.93 | 0.93 |
| tai60b-261117 | time 50% | 0.50 | 0.97 | 0.86 | 0.93 | 0.98 | 0.93 | 0.94 | 0.95 | 0.93 |
| tai60b-261117 | depth <= 11 | 0.50 | 0.97 | 0.88 | 0.94 | 0.97 | 0.88 | 0.92 | 0.85 | 0.86 |
| tai50b-261116 | time 25% | 0.49 | 0.91 | 0.55 | 0.92 | 0.92 | 0.84 | 0.80 | 0.80 | 0.89 |
| tai50b-261116 | time 50% | 0.48 | 0.96 | 0.63 | 0.93 | 0.93 | 0.96 | 0.93 | 0.95 | 0.95 |
| tai50b-261116 | depth <= 6 | 0.55 | 0.93 | 0.55 | 0.94 | 0.94 | 0.81 | 0.80 | 0.83 | 0.81 |
| tai50b-r9 | time 25% | 0.48 | 0.86 | 0.60 | 0.85 | 0.84 | 0.65 | 0.54 | 0.57 | 0.66 |
| tai50b-r9 | time 50% | 0.52 | 0.90 | 0.71 | 0.85 | 0.86 | 0.81 | 0.82 | 0.79 | 0.85 |
| tai50b-r9 | depth <= 5 | 0.53 | 0.86 | 0.70 | 0.83 | 0.83 | 0.80 | 0.78 | 0.82 | 0.88 |
| tho40ab-32 | time 25% | 0.51 | 0.98 | 0.88 | 0.91 | 0.90 | 0.81 | 0.84 | 0.76 | 0.90 |
| tho40ab-32 | time 50% | 0.51 | 0.99 | 0.97 | 0.89 | 0.92 | 0.96 | 0.96 | 0.96 | 0.99 |
| tho40ab-32 | depth <= 5 | 0.44 | 1.00 | 0.95 | 0.85 | 0.89 | 0.83 | 0.88 | 0.89 | 0.74 |
| tho40ab-0 | time 25% | 0.49 | 0.99 | 0.96 | 0.89 | 0.93 | 0.92 | 0.89 | 0.92 | 0.94 |
| tho40ab-0 | time 50% | 0.51 | 0.99 | 0.97 | 0.87 | 0.91 | 0.98 | 0.96 | 0.98 | 0.99 |
| tho40ab-0 | depth <= 5 | 0.48 | 1.00 | 0.95 | 0.86 | 0.91 | 0.83 | 0.74 | 0.81 | 0.92 |
| tho40-dnn | time 25% | 0.46 | - | 0.87 | 0.92 | - | 0.93 | - | 0.92 | - |
| tho40-dnn | time 50% | 0.49 | - | 0.88 | 0.94 | - | 0.94 | - | 0.93 | - |
| tho40-dnn | depth <= 5 | 0.46 | - | 0.91 | 0.95 | - | 0.92 | - | 0.93 | - |
| tho30-dnn | time 25% | 0.51 | - | 0.85 | 0.90 | - | 0.86 | - | 0.90 | - |
| tho30-dnn | time 50% | 0.49 | - | 0.94 | 0.92 | - | 0.94 | - | 0.96 | - |
| tho30-dnn | depth <= 4 | 0.47 | - | 0.92 | 0.84 | - | 0.84 | - | 0.91 | - |
| nug30-trace-a | time 25% | 0.45 | 0.99 | 0.80 | 0.92 | 0.96 | 0.88 | 0.99 | 0.89 | 0.99 |
| nug30-trace-a | time 50% | 0.47 | 1.00 | 0.88 | 0.93 | 0.97 | 0.86 | 0.99 | 0.91 | 1.00 |
| nug30-trace-a | depth <= 3 | 0.54 | 0.99 | 0.78 | 0.92 | 0.94 | 0.82 | 0.87 | 0.85 | 0.92 |
| tai50b-261116 -> tai50b-r9 | run | 0.48 | 0.84 | 0.90 | 0.79 | 0.85 | 0.95 | 0.94 | 0.95 | 0.95 |
| tai50b-r9 -> tai50b-261116 | run | 0.49 | 0.89 | 0.89 | 0.90 | 0.90 | 0.97 | 0.98 | 0.97 | 0.98 |
| tho40ab-0 -> tho40ab-32 | run | 0.49 | 0.97 | 0.99 | 0.94 | 0.96 | 0.97 | 0.99 | 0.98 | 0.99 |
| tho40-dnn -> tho40ab-32 | run | 0.49 | - | 0.92 | 0.94 | - | 0.96 | - | 0.93 | - |

## Summary: Spearman within parent with the shortfall

| test | split | random | qpb | pc | lin | lin+q | gbt | gbt+q | gbt-w | gbt-w+q |
|---|---|---|---|---|---|---|---|---|---|---|
| tai60b-261117 | time 25% | -0.00 | 0.39 | 0.32 | 0.36 | 0.39 | 0.33 | 0.33 | 0.35 | 0.35 |
| tai60b-261117 | time 50% | 0.00 | 0.41 | 0.31 | 0.37 | 0.42 | 0.38 | 0.38 | 0.39 | 0.38 |
| tai60b-261117 | depth <= 11 | 0.00 | 0.39 | 0.31 | 0.36 | 0.39 | 0.31 | 0.34 | 0.28 | 0.29 |
| tai50b-261116 | time 25% | -0.00 | 0.31 | 0.05 | 0.31 | 0.31 | 0.25 | 0.21 | 0.21 | 0.29 |
| tai50b-261116 | time 50% | -0.02 | 0.34 | 0.12 | 0.31 | 0.31 | 0.34 | 0.31 | 0.32 | 0.33 |
| tai50b-261116 | depth <= 6 | 0.04 | 0.32 | 0.05 | 0.33 | 0.33 | 0.21 | 0.20 | 0.23 | 0.22 |
| tai50b-r9 | time 25% | -0.01 | 0.25 | 0.06 | 0.23 | 0.23 | 0.10 | 0.03 | 0.05 | 0.10 |
| tai50b-r9 | time 50% | 0.01 | 0.27 | 0.13 | 0.23 | 0.24 | 0.20 | 0.21 | 0.19 | 0.23 |
| tai50b-r9 | depth <= 5 | 0.01 | 0.25 | 0.15 | 0.23 | 0.23 | 0.20 | 0.19 | 0.22 | 0.26 |
| tho40ab-32 | time 25% | 0.01 | 0.39 | 0.31 | 0.32 | 0.32 | 0.26 | 0.28 | 0.23 | 0.33 |
| tho40ab-32 | time 50% | 0.00 | 0.38 | 0.37 | 0.31 | 0.33 | 0.37 | 0.36 | 0.36 | 0.39 |
| tho40ab-32 | depth <= 5 | -0.04 | 0.36 | 0.33 | 0.25 | 0.28 | 0.25 | 0.28 | 0.29 | 0.20 |
| tho40ab-0 | time 25% | -0.01 | 0.39 | 0.37 | 0.31 | 0.34 | 0.34 | 0.33 | 0.34 | 0.36 |
| tho40ab-0 | time 50% | 0.01 | 0.38 | 0.37 | 0.29 | 0.32 | 0.37 | 0.36 | 0.37 | 0.38 |
| tho40ab-0 | depth <= 5 | -0.01 | 0.36 | 0.33 | 0.25 | 0.29 | 0.25 | 0.18 | 0.24 | 0.31 |
| tho40-dnn | time 25% | -0.03 | - | 0.35 | 0.42 | - | 0.42 | - | 0.43 | - |
| tho40-dnn | time 50% | -0.00 | - | 0.36 | 0.45 | - | 0.44 | - | 0.44 | - |
| tho40-dnn | depth <= 5 | -0.01 | - | 0.29 | 0.33 | - | 0.30 | - | 0.32 | - |
| tho30-dnn | time 25% | 0.01 | - | 0.34 | 0.39 | - | 0.36 | - | 0.40 | - |
| tho30-dnn | time 50% | -0.02 | - | 0.42 | 0.40 | - | 0.43 | - | 0.44 | - |
| tho30-dnn | depth <= 4 | -0.02 | - | 0.32 | 0.25 | - | 0.26 | - | 0.32 | - |
| nug30-trace-a | time 25% | -0.03 | 0.46 | 0.27 | 0.40 | 0.44 | 0.33 | 0.46 | 0.33 | 0.46 |
| nug30-trace-a | time 50% | -0.03 | 0.38 | 0.28 | 0.33 | 0.36 | 0.29 | 0.38 | 0.31 | 0.38 |
| nug30-trace-a | depth <= 3 | 0.04 | 0.37 | 0.19 | 0.32 | 0.33 | 0.23 | 0.28 | 0.25 | 0.31 |
| tai50b-261116 -> tai50b-r9 | run | -0.02 | 0.23 | 0.32 | 0.20 | 0.25 | 0.35 | 0.34 | 0.35 | 0.35 |
| tai50b-r9 -> tai50b-261116 | run | -0.02 | 0.28 | 0.32 | 0.30 | 0.31 | 0.37 | 0.37 | 0.37 | 0.38 |
| tho40ab-0 -> tho40ab-32 | run | -0.01 | 0.42 | 0.43 | 0.38 | 0.41 | 0.41 | 0.43 | 0.42 | 0.43 |
| tho40-dnn -> tho40ab-32 | run | -0.01 | - | 0.37 | 0.38 | - | 0.40 | - | 0.39 | - |
