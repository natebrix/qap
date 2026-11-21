# The root lookahead for the grids

Nov 21 2026. In progress. [m] = measured here, [e] = estimate or inference.

## What it is

At a grid's root every row (facility) is a candidate, and under symmetry a row has only the children the
symmetry leaves (the columns of J1: 9 at nug30, 10 at tho30, 12 at tho40). `root.patch` makes the DNN lookahead
bound only those children (`QAP_LA=1:150:0:0 QAP_LA_SHORTLIST=a`, the root given a DNN bound: `grid_root.par`).
All rows cost about half a CPU hour at n = 30 and 3.8 core-hours at tho40 [m].

## The score: predicted tree size, not the sum of shortfalls

- **The sum of the children's shortfalls is the wrong score** [m]: at nug30 its best row (r6) gives 4,875 nodes
  against 4,365 for production's r0, which it ranks seventh (`runs/nug30_r6`, `runs/nug30_prod`).
- **A child's shortfall predicts the log of its subtree** [m]: over the 12 depth-1 children of tho40's r26,
  the shortfall after 150 iterations against ln(estimated nodes) from the estimate tho40-est2-261119:
  Pearson +0.94, R^2 0.88 (`est_by_child.py`). Over the 9 children of nug30's r0 against their subtrees in the
  whole tree: R^2 0.90; the 10 of tho30's r5: R^2 0.49 (`expscore.py`).
- **So a row is scored by the tree it predicts**: the sum over its children of exp(a + b x shortfall), a and b
  fitted where a tree or an estimate gives the children's subtree sizes.

## Predictions, written before the trees ran (17:55 UTC)

Fits from one row's children each (`expscore.py`), so the predictions outside that row's range of shortfalls
are extrapolations.

| grid | row | by the sum of shortfalls | by predicted size | predicted nodes | actual |
|---|---|---|---|---|---|
| nug30 | r0 (production) | rank 7 | rank 4 | 4,208 (in sample) | 4,365 and 4,285 [m] |
| nug30 | r23 | 4 | 1 | 3,211 | running |
| nug30 | r6 | 1 | 17 | 9,168 | 4,875 [m, known before the fit] |
| nug30 | r2 | 30 | 29 | 16,404 | running |
| tho30 | r5 (production) | 2 | 1 | 13,264 (in sample) | 14,815 [m] |
| tho30 | r25 | 1 | 30 | 333,945 | to run |
| tho30 | r4 | 5 | 2 | 16,613 | to run |

r6 already says the exponential over-predicts beyond the fitted range (its largest child is 150 short, the
fit's largest 105): right direction, wrong size. The trees are run with the trace on (`tree_tr.sh`) so that
every tree adds its children to the fit.

Decision rule, fixed with the first trees (`tree.sh`): over the rows run, the rank correlation between the
score and the tree's nodes; the root lookahead pays on a grid if its best row's tree has at most 95% of
production's nodes; a score carries no information if the correlation is 0.3 or less.
