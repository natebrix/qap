# The root lookahead for the grids

Nov 21 2026. [m] = measured here, [e] = estimate or inference.

**Scoring every row at a grid's root is cheap, and scoring a row by the tree its children predict sorts good
rows from bad ones (nug30: 3,623 to 10,165 nodes over six rows; tho30: 13,812 to more than 21,000 over
three). It does not order the best few: its error on a row's tree is a factor of about 1.3. At tho40 it puts production's row 15th of 40 and predicts about
30% fewer nodes for row 2: a bet at the edge of that error, which a second estimate can settle.**

## What it is

At a grid's root every row (facility) is a candidate, and under symmetry a row has only the children the
symmetry leaves (the columns of J1: 9 at nug30, 10 at tho30 and sko49, 12 at tho40, 15 at wil50).
`root.patch` makes the DNN lookahead bound only those children; `QAP_LA=1:150:0:0 QAP_LA_SHORTLIST=a` scores
every row; the root needs a DNN bound for the children's warm start (`grid_root.par`: `dnn.depths` from 0).
The patch goes on cognac 0e6c2ecf + `experiments/line_signals/experiment_build.patch` + `../lsla/lsla.patch` +
`../tail/tail.patch`.

| instance | rows x children | CPU for all rows | children fathomed within 150 iterations |
|---|---|---|---|
| nug30 | 30 x 9 | 1,726 s [m] | 4 |
| tho30 | 30 x 10 | 1,813 s [m] | 1 |
| tho40 | 40 x 12 | 3.8 core-hours [m `../lsla/grids`] | none |
| sko49 | 49 x 10 | 8.9 core-hours [m `../lsla/grids`] | none |
| wil50 | 50 x 15 | 14.8 core-hours [m `../lsla/grids`] | none |

## The score

- **The sum of the children's shortfalls is the wrong score** [m]. Its first row at nug30 (r6) gives 4,875
  nodes against production's 4,365, and at tho30 (r25) at least 21,095 trace lines against 14,900.
- **A child's shortfall predicts the log of its subtree** [m]: over the 12 depth-1 children of tho40's r26
  against the estimate tho40-est2-261119's subtree under each (`est_by_child.py`): Pearson +0.94, R^2 0.88.
  Over 42 children of five nug30 trees: R^2 0.90 (exponential), residual a factor 1.56 a child (`models.py`).
- **So a row is scored by the tree it predicts**: 1 + the sum over its children of the size the fit gives
  each (`rootfit.py`, `models.py`). One child far from fathoming outweighs several near it, which the plain
  sum cannot see.
- **The form**: exponential (ln n = a + b s) and power law (ln n = a + p ln s) agree inside the range of
  shortfalls the fit saw and part beyond it; the exponential overshoots there (r6, whose worst child is 150
  short against the first fit's 105: predicted 2.6 times the tree). Neither is better at the top of the
  ranking. A term for symmetry left at the child helps the exponential a little (1.44 to 1.34).
- **The scatter is not the lookahead's**: the log of a child's subtree correlates 0.96 with the lookahead's
  shortfall, 0.95 with the child's own depth-1 bound in the tree and 0.96 with its branching children at
  depth 2 [m, 33 children]. One level of bound information carries a factor of about 1.5 a child, and a
  row's tree is mostly its one or two worst children.

## nug30: six trees (`tree.sh`, `tree_tr.sh`; the recipe as it is, 2 slots, incumbent 6126)

| root row | rank by predicted size: power / exponential, of 30 | rank by the sum of shortfalls | solver's nodes | against production |
|---|---|---|---|---|
| r23 | 2 / 1 | 4 | **3,623** | 0.84 |
| r0 (production's: the QPB row lookahead) | 5 / 4 | 7 | 4,285 and 4,365 (two runs) | 1 |
| r1 | 3 / 2 | 5 | 4,566 | 1.06 |
| r4 | 1 / 3 | 2 | 4,616 | 1.07 |
| r6 | 4 / 17 | 1 | 4,833 and 4,875 | 1.12 |
| r2 | 30 / 29 | 30 | 10,165 | 2.35 |

- **The root row matters**: a factor 2.8 between the best and the worst tree run [m].
- **The score sorts good from bad**: its top four hold the best tree and nothing more than 7% over
  production; the 2.35x tree is last or next to last [m].
- **It does not order the top four**: of its three first picks one beats production by 16% and two lose
  by 6 and 7% [m]. Predicting each row from the other rows' children: within a factor 1.3 (power) to 1.46
  (exponential) of the tree.
- By the rule fixed before the trees: the best row's tree is at most 95% of production's only if the pick is
  r23 (the exponential's first, the power law's second); Spearman(score, nodes) over the six is +0.83
  (exponential) and +0.43 (power), against +0.03 for the sum of shortfalls [m, the fits of `models.py` on the
  five traced trees].
- The predictions written before the trees ran (17:55 UTC, commit 483bd294): r23 3,211 (3,623), r2 16,404
  (10,165), tho30 r25 333,945 (not closed: at least 21,095).

## tho30

| root row | rank by predicted size | by the sum | trace lines | |
|---|---|---|---|---|
| r5 (production) | 1 | 2 | 14,900 (14,815 solver's nodes) | `../grid/runs/tho30_prod` |
| r25 | 26 to 30 | 1 | **at least 21,095**, not closed | stopped by the two-hour limit on background commands with 6 units open and not restarted; `runs/tho30_r25_tr/nodes.partial.txt.gz` |
| r4 | 2 | 5 | **13,887** (13,812 solver's nodes): 0.93 of production | `runs/tho30_r4_tr` |

- Production's row is the size model's first and r4 its second (predicted 1.15 to 1.25 times r5); r4's tree is
  7% smaller than r5's [m]. As at nug30: the top picks are within some 15% of each other and the model
  does not order them. r25, first by the sum of shortfalls, is at least 42% worse [m].
- **The curve is not the same under every row here** [m]: r5's children have larger subtrees at the same
  shortfall than r25's (5:0, 3,198 short, 1,214 nodes; 25:2, 3,528 short, 335), and the pooled fit's residual is
  a factor 2 a child (nug30: 1.5). The first fit (r5's ten children alone, a narrow range) had twice the slope
  and predicted 334,000 for r25.

## The size score at depth 1 (`la_tree.sh`, `QAP_LA_SCORE`)

The size model holds one level down: over the depth-2 nodes of the production trees, ln(nodes under a node)
against its own shortfall has R^2 0.67 to 0.71 and the same slope in two tho30 trees (0.0011 and 0.0010;
nug30 0.031) [m]. `QAP_LA_SCORE=exp:B` scores a line in the lookahead by the sum of exp(B x shortfall).

| the lookahead at the depth-1 nodes (`QAP_LA=4:150:1:1`) | nug30 nodes | tho30 nodes |
|---|---|---|
| production (none) | 4,285 and 4,365 | 14,815 |
| the sum of shortfalls, the build of `../grid` (every column at a node with symmetry left) | | 11,100 |
| the sum of shortfalls, this build (a row's real children only) | 4,037 | **10,083** |
| the size score (exp:0.031, exp:0.0011), this build | 3,921 | 10,210 |

- **At depth 1 the two scores do the same** [m]: 3% apart at nug30 and 1% at tho30 (the other way), under
  the 5% set as the bar. They choose the same line at 6 of 9 and 8 of 10 nodes. The size score matters at the
  root, where rows differ widely and the sum picks badly, not here.
- **Restricting the lookahead to a row's real children is worth 9% at tho30** [m]: 10,083 against 11,100 with
  the same score. With symmetry left the old build scored a row over columns that are not its children.
- **The depth-1 lookahead is now measured on two grids** [m]: 6 to 10% fewer nodes at nug30, 32% at tho30.
  Its cost (3,581 and about 5,160 CPU s) is more than it saves on trees this small (`../budget/README.md`
  for the large ones).

## The open grids (scores from `../lsla/grids`)

| instance | production's root row | its rank by predicted size | best row | predicted tree, best / production's | calibration |
|---|---|---|---|---|---|
| tho40 | r26 | 15 of 40 | r2 | 0.68 (nodes), 0.71 (CPU); the same under the power law | the Nov 19 estimate under r26's 12 children: R^2 0.88 to 0.92, a factor 1.3 to 1.4 a child [m] |
| sko49 | r41 | 4 to 9 of 49, by the slope | r45 at every slope | 0.87 to 0.5 | none: its estimate is a depth-3 sample and gives no subtree under a depth-1 node [e] |
| wil50 | r21 | 2 of 50 at every slope | r45 at every slope | 0.98 to 0.84 | none (as sko49) |

- **tho40**: r2's children lie inside the range of shortfalls the fit saw (10,039 to 13,466 against 9,819 to
  13,366), so the form of the model does not matter there. The other rows the model puts within 0.85 of r26:
  r35, r16, r4, r30, r24, r0, r12. The worst row (r6) is predicted at 2.5 to 2.9 times r26.
- With nug30's experience (top picks within +-15% of each other in fact when predicted within 30%), r2 at
  0.7 is a bet, not a finding [e].
- **sko49**: r45 is first whatever the slope; how far ahead depends on a slope nothing here fixes.
- **wil50**: production's r21 and r45 are the top two at every slope, their children's shortfalls nearly the
  same cell for cell, and the third row is predicted 1.3 to 1.5 times r21. Nothing to gain at wil50's root:
  the QPB row lookahead already has it. The laptop's root study agrees (`../../wil50_feasibility/sampling`
  README, "Root lines": 45 best, 21 next).
- The slope: "0.87 to 0.5" and "0.98 to 0.84" are the model with ln(nodes) = k x shortfall / the mean
  shortfall for k from 3 to 16 (tho40's calibrated k is 12; nug30's about 4).

## What this adds up to

1. **Always worth its cost as insurance against a bad root row.** Hours of CPU against runs of thousands,
   and a bad row costs a factor of two or more (nug30 2.35x; tho40's worst predicted 2.5x).
2. **As a selector it needs a second stage.** The score makes a shortlist of three or four rows; choosing
   among them takes something that measures subtrees, which is the estimator: a reduced estimate per
   shortlisted row, or at tho40 one more full estimate with r2 against the Nov 19 one (its standard error is
   about 7% on CPU, so a 30% difference is about three standard errors of the comparison).
3. **What cognac would need**: a setting that forces the root row (`QAP_FORCE_LINE=0:rN` in the experiment
   build). The scoring runs offline and hands over one integer; the lookahead code can stay out.

## Not done

- A calibration for sko49 and wil50 (their estimates would have to be totalled under each depth-1 node).
- The root's row and the depth-1 lookahead together on one tree.
- Why the curve differs between rows at tho30 and not at nug30.

## Files

`root.patch`, `grid_root.par`; `tree.sh`, `tree_tr.sh` (a tree with the root row forced; `_tr` with the node
trace); `la_tree.sh` (a tree with the lookahead at depth 1 and a chosen score); `expscore.py` (one row's fit), `rootfit.py` (pooled, leave one row out), `models.py` (the forms
compared), `est_by_child.py` (an estimate's rest under each depth-1 node against the lookahead's shortfalls);
`runs/PROBLEM_root150/` (the scores: `children.txt`, `la.log`), `runs/PROBLEM_rN[_tr]/` (the trees).
