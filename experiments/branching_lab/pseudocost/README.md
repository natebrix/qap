# Pseudocosts for the DNN bound: can the run's own history predict a child's first solve? (Nov 19 2026)

**Answer: no, not well enough to build.** Trained only on what a run had already bounded, no
pseudocost table, linear model or gradient-boosted model ranks a parent's children better than the
child's own QPB bound. The child's QPB bound gives within-parent AUC 0.86-1.00 for "this child's first
solve will not fathom". The best learned model ties it at most, 0.98 against 0.97 on tai60b [m]. On the
line-truth nodes the learned line score has mean Spearman -0.00 (gradient-boosted trees) or +0.16
(pseudocost table), against the QPB lookahead's +0.34 and the children's first DNN solve's +0.80 [m].
That meets the kill rule on both counts.

What the history does learn is the tree's own neighbourhood. Re-solving the same subproblems (another run
of the same tree; random cross-validation folds whose siblings share all but one fix) looks good. On new
subproblems the advantage goes away.

One side finding is worth keeping. At the truth nodes, the parent's DNN marginal x_ij separates a line's
open children from its fathomed ones as well as or better than QPB's bound + U_ij: AUC 0.95 against 0.92
on tho40, 0.77 against 0.82 on tho30 [m]. It still cannot rank lines, because every line's marginals sum
to 1.

All numbers are **[m]** unless marked [e]. Scripts and their outputs are in this directory:
`build_table.py` -> `data/*.pkl` (not committed, `.gitignore`) and `build_table_out.txt`; `describe.py` ->
`describe_out.md`; `predict.py` -> `predict_out.md`; `overlap.py` -> `overlap_out.md`; `truth.py` ->
`truth_out.md`. `common.py` holds the features, the models and the metrics. Run them with
`OMP_NUM_THREADS=1 python3 <script>`. scikit-learn 1.9.1 was installed with pip for this.

## The data (step 1)

There is one row per (parent, child) in each trace (`build_table.py`). The parent of a child is the node
whose fixes are the child's less the last pair.

The orientation was checked against QAPSOLN (`build_table_out.txt`). A fix `r:c` assigns row r of the
instance's first matrix A to index c of the second matrix B, with objective sum A_ij B_p(i)p(j):
- On tai50b and tai60b, the best known permutation's chain runs to depth 47 and 57 in the traces as
  row -> column, and to depth 1 as its inverse. On these instances A is the symmetric clustered
  distance-like matrix and B holds the asymmetric flows.
- On tho40t.dat the best known permutation is the inverse of `QAPSOLN/tho40.sln` (the file is the
  transposed copy). On tho/nug, A holds the flows and B the grid distances.

Two older DNN traces (`QAP_DNN_TRACE`, Nov 8: tho30 and tho40, the run the tho40 truth nodes came from)
add 14,783 and 29,245 children. They have no QPB bound and no line, and their parents were rebuilt from
the sets of fixes.

| run | children | parents | open (first solve missed) | of the open: cut phase fathomed / branched | depths |
|---|---|---|---|---|---|
| tai60b-261117 | 331,147 | 7,495 | 6.2% | 64% / 36% | 1-56 |
| tai50b-261116 | 8,214 | 213 | 6.9% | 63% / 37% | 1-46 |
| tai50b-r9 | 6,000 | 143 | 7.1% | 67% / 33% | 1-23 |
| tho40ab-32 | 22,649 | 644 | 2.8% | - / 100% | 3-8 |
| tho40ab-0 | 22,438 | 637 | 2.8% | - / 100% | 3-8 |
| tho40 Nov 8 (DNN trace) | 29,245 | 959 | 3.3% | - | 2-8 |
| tho30 Nov 8 (DNN trace) | 14,783 | 605 | 4.0% | - | 2-8 |
| nug30-trace-a | 3,687 | 141 | 3.8% | - / 100% | 1-6 |

(`tai50b-r9`'s trace holds 6,001 lines, not the 6,629 nodes of the brief. `nug30-c7a`'s trace differs
from `nug30-trace-a`'s, but the two give the same table.)

**Censoring.** Targets were chosen so that the solver's stopping rules decide as little as possible:
- `open`: the first solve did not fathom (outcome C or B). This is the solver's own label, and the
  cutoff is fixed in every run used (the incumbent is the optimum + 1 from the start). A give-up that
  fires on a child that more iterations would have fathomed counts as open. That is production's
  behaviour, and the label the lookahead needs.
- `shortfall` = max(0, cutoff - first-solve bound). This is **exact (0) for every fathomed child**: the
  solve stops once it passes the cutoff, so the recorded bound is only a floor, but its positive
  shortfall is not censored. For the open children it is the bound where the solve stopped (iteration
  cap or the give-up's 5x-geometric extrapolation), so it overstates the converged shortfall. The
  ranking metrics use it only through Spearman, and every score is judged against the same values.
- The gain is used only capped at the cutoff: (min(child DNN, cutoff) - parent bound) = gap - shortfall.
  The uncapped gain of a fathomed child is a stopping artefact.
- Cut-phase bounds are never a target (a different quantity). They appear only as the parent's flag
  "its cut phase ran".
- No CPU, iteration or status field is a feature. P (popped) and E (enumerated) children are dropped.

The parent's final bound is often its cut-phase bound, close to the cutoff. 80% of tai60b's open
children have a first-solve bound *below* their parent's final bound. A raw "gain" is therefore mostly
negative noise, and the parent's gap is a poor normaliser. The shortfall is normalised by |cutoff|, and
models fit ylog = log1p(1e4 x shortfall / |cutoff|).

Variance explained, out of fold (5 folds grouped by parent; `describe_out.md`; R2 of ylog, total /
within parent). Within parent means the target less its parent's mean, which is what ranking a parent's
children needs.

| predictor | tai60b | tai50b-261116 | tai50b-r9 | tho40ab-32 | tho40 Nov 8 | tho30 Nov 8 | nug30 |
|---|---|---|---|---|---|---|---|
| between parents (ceiling of any parent-level feature) | 0.04 | 0.26 | 0.28 | 0.28 | 0.43 | 0.28 | 0.31 |
| row | 0.14 / 0.14 | 0.11 / 0.04 | 0.09 / 0.00 | 0.06 / 0.09 | 0.05 / 0.05 | 0.06 / -0.01 | 0.00 / 0.02 |
| column | 0.10 / 0.10 | 0.07 / 0.08 | 0.05 / 0.07 | 0.01 / 0.00 | 0.01 / -0.01 | 0.06 / 0.04 | 0.01 / 0.00 |
| row + column | 0.24 / 0.24 | 0.16 / 0.11 | 0.13 / 0.07 | 0.06 / 0.09 | 0.07 / 0.05 | 0.09 / 0.03 | 0.01 / 0.02 |
| row + column + pair (pseudocost table) | 0.62 / 0.63 | 0.40 / 0.44 | 0.34 / 0.37 | 0.11 / 0.15 | 0.08 / 0.06 | 0.31 / 0.31 | 0.08 / 0.11 |
| depth | 0.01 / 0 | 0.10 / 0 | 0.10 / 0 | 0.08 / 0 | 0.14 / 0 | 0.14 / 0 | 0.13 / 0 |
| parent's first-solve shortfall | 0.00 / 0 | -0.01 / 0 | -0.02 / 0 | 0.25 / 0 | 0.38 / 0 | 0.21 / 0 | 0.13 / 0 |
| child's QPB shortfall | 0.13 / 0.13 | 0.09 / 0.01 | 0.05 / -0.02 | 0.50 / 0.38 | - | - | 0.63 / 0.51 |
| structure (flow/distance sums, nnz, clusters, Gilmore-Lawler term, interaction with the fixed pairs) | 0.63 / 0.65 | 0.40 / 0.45 | 0.34 / 0.37 | 0.13 / 0.15 | 0.14 / 0.13 | 0.32 / 0.32 | 0.13 / 0.12 |
| agrees with the best known solution | 0.14 / 0.15 | 0.05 / 0.10 | 0.03 / 0.09 | 0 / 0 | 0 / 0 | - | 0.05 / 0.02 |
| all but the child's QPB (GBT) | 0.72 / 0.72 | 0.56 / 0.57 | 0.50 / 0.50 | 0.77 / 0.71 | 0.74 / 0.61 | 0.74 / 0.69 | 0.34 / 0.43 |
| all, with the child's QPB (GBT) | 0.72 / 0.73 | 0.58 / 0.58 | 0.51 / 0.51 | 0.85 / 0.80 | - | - | 0.80 / 0.73 |

(`open` as the target gives the same picture, `describe_out.md`.)

The capped gain (`gainrel`) is mostly the parent's gap: a fathomed child's capped gain equals its
parent's gap. The parent's shortfall explains 0.93-0.94 of its total variance on tho30 and tho40
(`describe_out.md`). Its within-parent part is the shortfall's with the sign flipped. Because so little
of it varies within a parent, its within-parent R2s are unstable and sometimes far below zero, so they
are not used.

How to read the table:
- On tai60b, 96% of the variance is within parents. The (row, column) pair identity alone explains 0.63
  of it, and so does the structure, which is a function of the pair. The tree is two mirrored chains with
  subtrees off them (`tai60b_postmortem`), so the same pair recurs under near-identical parents.
- Row alone or column alone explains little (0.00-0.14). Depth and the parent's shortfall explain
  nothing within a parent, by construction.
- The child's QPB bound is weak in R2 on tai (0.13 or less), where QPB is 13-29% under the cutoff, but it
  orders children well (AUC 0.97 below). R2 on a log shortfall rewards magnitudes, and ranking needs
  orders.
- These random-fold R2s are optimistic: a fold's test parents have siblings and cousins in its training
  folds. Step 2 is the honest test.

## Out of sample (step 2: `predict.py`, `overlap.py`)

Scores, with larger meaning more likely open:
- `qpb`: the child's QPB shortfall. Rule 9 has it on its two lines, and the trace records it for every
  child.
- `pc`: a pseudocost table, with row, column and pair effects shrunk by ridge (alpha 3).
- `lin`: ridge on the features.
- `gbt`: gradient-boosted trees (scikit-learn HistGB, 300 trees, row and column as categories).
- `gbt-w`: the same model fit to the target less its parent's mean.
- `+q`: with the child's QPB as a feature.

Splits:
- time: train on the children traced before the first 25% or 50% of the trace's lines, test on the
  parents traced after.
- depth: train on the depths at the top that hold 30% of the children, test below.

Metrics, mean over test parents with 3 or more children: AUC of open against fathomed, and Spearman with
the shortfall. `predict_out.md` also has R-precision and the recall and precision of a shortlist of each
parent's top quarter.

AUC within parent:

| test | split | qpb | pc | lin | lin+q | gbt | gbt+q | gbt-w+q |
|---|---|---|---|---|---|---|---|---|
| tai60b | time 25% | **0.97** | 0.90 | 0.94 | **0.97** | 0.92 | 0.92 | 0.93 |
| tai60b | time 50% | 0.97 | 0.86 | 0.93 | **0.98** | 0.93 | 0.94 | 0.93 |
| tai60b | depth <= 11 | **0.97** | 0.88 | 0.94 | **0.97** | 0.88 | 0.92 | 0.86 |
| tai50b-261116 | time 25% | 0.91 | 0.55 | **0.92** | **0.92** | 0.84 | 0.80 | 0.89 |
| tai50b-261116 | time 50% | **0.96** | 0.63 | 0.93 | 0.93 | **0.96** | 0.93 | 0.95 |
| tai50b-261116 | depth <= 6 | 0.93 | 0.55 | **0.94** | **0.94** | 0.81 | 0.80 | 0.81 |
| tai50b-r9 | time 25% | **0.86** | 0.60 | 0.85 | 0.84 | 0.65 | 0.54 | 0.66 |
| tai50b-r9 | time 50% | **0.90** | 0.71 | 0.85 | 0.86 | 0.81 | 0.82 | 0.85 |
| tai50b-r9 | depth <= 5 | 0.86 | 0.70 | 0.83 | 0.83 | 0.80 | 0.78 | **0.88** |
| tho40ab-32 | time 25% / 50% / depth | **0.98 / 0.99 / 1.00** | 0.88 / 0.97 / 0.95 | 0.91 / 0.89 / 0.85 | 0.90 / 0.92 / 0.89 | 0.81 / 0.96 / 0.83 | 0.84 / 0.96 / 0.88 | 0.90 / **0.99** / 0.74 |
| tho40ab-0 | time 25% / 50% / depth | **0.99 / 0.99 / 1.00** | 0.96 / 0.97 / 0.95 | 0.89 / 0.87 / 0.86 | 0.93 / 0.91 / 0.91 | 0.92 / 0.98 / 0.83 | 0.89 / 0.96 / 0.74 | 0.94 / **0.99** / 0.92 |
| nug30 | time 25% / 50% / depth | **0.99 / 1.00 / 0.99** | 0.80 / 0.88 / 0.78 | 0.92 / 0.93 / 0.92 | 0.96 / 0.97 / 0.94 | 0.88 / 0.86 / 0.82 | **0.99** / 0.99 / 0.87 | **0.99 / 1.00** / 0.92 |
| tho40 Nov 8 (no QPB) | time 25% / 50% / depth | - | 0.87 / 0.88 / 0.91 | 0.92 / 0.94 / 0.95 | - | 0.93 / 0.94 / 0.92 | - | - |
| tho30 Nov 8 (no QPB) | time 25% / 50% / depth | - | 0.85 / 0.94 / 0.92 | 0.90 / 0.92 / 0.84 | - | 0.86 / 0.94 / 0.84 | - | - |

Spearman with the shortfall follows the same pattern. qpb gives 0.39-0.41 on tai60b and 0.25-0.34 on
tai50b; the best model gives 0.42 (lin+q, tai60b 50%), and otherwise matches or falls below qpb
(`predict_out.md`). Two other results:
- The learned models are better at the very top of the list on tai50b. R-precision (the k best-scored
  children, k = the parent's open children) is 0.80 for gbt against 0.48 for qpb (tai50b-261116, time
  50%), and 0.62 against 0.32 on tai50b-r9. They are not better on tai60b (0.73 against 0.71), tho40
  (0.67-0.70 against 0.82-0.87) or nug30.
- Without QPB (the Nov 8 DNN traces), the models reach AUC 0.86-0.95. That is what history alone gives,
  short of the child's QPB.

**Across runs** (`overlap.py`). Some test children are subproblems that the training run also bounded,
with the same set of fixes. On those, a model can recall the node instead of predicting it:

| train -> test | test parents | parents with only new children | score | AUC, all | AUC, new parents | Spearman, new parents |
|---|---|---|---|---|---|---|
| tai50b-261116 -> tai50b-r9 | 143 | 66 | qpb | 0.84 | **0.92** | **0.28** |
| | | | gbt | **0.95** | 0.89 | 0.26 |
| tai50b-r9 -> tai50b-261116 | 205 | 128 | qpb | 0.89 | 0.95 | 0.33 |
| | | | gbt | **0.97** | **0.97** | **0.34** |
| tho40ab-0 -> tho40ab-32 (the same subtree) | 644 | 32 | qpb | 0.97 | 0.97 | 0.43 |
| | | | pc | 0.99 | 0.98 | 0.43 |
| tho40 Nov 8 -> tho40ab-32 (another subtree) | 644 | all | gbt | 0.96 (qpb within tho40ab-32: 0.98-1.00) | | |

Across the two tai50b runs, the model's apparent edge (0.95-0.97 against 0.84-0.89) comes from the
shared upper tree. On the parents whose children are all new, it is 0.89 against 0.92 and 0.97 against
0.95: a tie.

## Line ranking at the truth nodes (step 3: `truth.py`)

**Method.**
- Training sets:
  - tho30: the Nov 8 tho30 DNN trace, less the truth node's subtree. No tho30 node trace exists; the
    tho30 traces are under `experiments/tho40_estimate/traces/`.
  - tho40: tho40ab-32, as the brief asked. Also the Nov 8 tho40 trace less the three truth subtrees, and
    both together.
  - tai60b: its trace less the subtree of 44:45,43:58.
- A child is the subproblem F + (r, c), whatever line creates it, so the features never name the line.
- The truth nodes lack the children's QPB bounds, so the model is the no-QPB gbt (and pc).
- A line's score is the sum over its children of the predicted shortfall. The truth is each line's
  subtree nodes, as in `../../../../wt-cognac/experiments/line_signals/signals_out.md`.

Line-level results:

| signal | nodes | mean Spearman | best line ranked 1st | its top line / best, mean (worst) |
|---|---|---|---|---|
| marginal score (rules 8, 9), from `signals_out.md` | 7 | -0.07 | 1/7 | 1.52x (1.99x) |
| QPB lookahead over every line (rule 4), from `signals_out.md` | 7 | +0.34 | 1/7 | 1.67x (3.42x) |
| (a) the children's first DNN solve, from `signals_out.md` | 6 | +0.80 | 5/6 | 1.00x (1.00x) |
| **gbt on the traces: predicted shortfall sum** | 7 | **-0.00** | 1/7 | 1.95x (3.38x) |
| **pc, the pseudocost table: predicted shortfall sum** | 7 | **+0.16** | 2/7 | 1.52x (2.29x) |

By node, mean Spearman (`truth_out.md`):

| node group | gbt | pc | QPB lookahead |
|---|---|---|---|
| tho30 (3 nodes) | +0.62 | +0.80 | +0.46 |
| tho40 with tho40ab-32 (3 nodes) | -0.49 | -0.30 | +0.46 |
| tai60b 43:58 (1 node) | -0.40 | -0.40 | -0.40 |

Training tho40 on the Nov 8 trace instead, or on both traces, is no better: gbt -0.50 and -0.51, pc
-0.62 and -0.74.

How to read these results:
- tho30 is the model's good case: the gbt ranking at 5:1 is +1.00. Its training set holds the subtrees of
  the other depth-1 orbits, which differ from the truth node only in where facility 5 sits, so this is
  close to recall.
- On tho40 the model fails at every node. It puts row14 and row7, the truth's 7th and 8th lines, first
  and second at 24:7,26:1,39:0, because it predicts their children all fathom.
- The test is small: seven nodes, two instances and a third with four lines. At this size one node moves
  a mean Spearman by about 0.15-0.3. The kill rule is not close, though: -0.00 and +0.16 against +0.34.

**Child level** (the 232-296 candidate children of each node, against their cold first solves):

| training | score | Spearman within a line | AUC open vs fathomed |
|---|---|---|---|
| tho30 (3 nodes) | gbt | +0.77 | 0.90 |
| | QPB's bound + U_ij (the parent's dual) | +0.76 | 0.82 |
| | the DNN marginal x_ij | +0.71 | 0.77 |
| tho40, tho40ab-32 (3 nodes) | gbt | +0.28 | 0.64 |
| tho40, both traces | gbt | +0.42 | 0.78 |
| | QPB's bound + U_ij | +0.51 | 0.92 |
| | the DNN marginal x_ij | +0.56 | 0.95 |

The marginals and U, which the solver already has at the node, beat the learned model on tho40. They
match it on tho30.

The marginal tells which of a line's children stay open. It cannot say how many, because x is doubly
stochastic and each line sums to 1. That is why the marginal score fails at lines (-0.07) while x_ij
does well within a line.

## Verdict (step 4)

The kill rule was set before the runs. Out of sample, the model does not beat the child's QPB bound at
ranking a parent's children; at best it ties (0.98 against 0.97 on tai60b, 0.97 against 0.95 across the
tai50b runs on new subproblems). On the truth sets it does not beat the QPB lookahead's +0.34 (-0.00
gbt, +0.16 pc). **A pseudocost or learned predictor of the DNN first solve is not worth building.**

There are three reasons:
1. What a run has seen predicts a child well only when the child repeats it: the same pair under a
   near-identical parent, as along tai60b's chains. There the child's QPB bound already ranks the
   children (AUC 0.97).
2. At a node where a line must be chosen, the candidate lines' children are new subproblems. The learned
   predictor then falls to or below the cheap signals the node already has.
3. The DNN bound's gain from a fix has no stable per-(row, column) component that transfers. Row or
   column alone explains 0-14% of the variance.

So a reliability-branching rule ("stop looking ahead once the table is reliable") has nothing reliable
to switch to. The DNN lookahead (`line_signals`, Verdict 4: 4-6 lines a side, about 150 iterations,
at the top depths of the big instances) remains the only signal that ranks lines.

## What was not tested

- Features from the node's own DNN state:
  - the marginals and U as model inputs, which the traces do not carry;
  - the children's warm-start residuals;
  - a few warm iterations, which is the lookahead itself.
- A ranking loss (LambdaMART) instead of a regression on the log shortfall. The within-parent fit
  (`gbt-w`) did not help.
- Truth-set children's QPB bounds, which were not recorded. The line-level QPB baseline is rule 4's
  score from the dumps.
- tai50b truth nodes: none exist.
- The tai60b truth is the rough subtree model of `tai60b_lines`, for one node and four lines.
- Whether the shortfall of open children would rank differently if the give-up were off. The lookahead's
  own uncapped solves (`line_signals` section 3) would say. Here every score was judged against the same
  censored values.

## Cheapest next step, if any

None for pseudocosts. If the side finding is wanted, it can be used inside the DNN lookahead, not in
place of it. Among a candidate line's children, solve first those with the largest marginal x_ij (AUC
0.95 for "stays open" on tho40, 0.77 on tho30). Stop scoring a line once its open children's shortfall
sum exceeds the best line's so far. The cost is zero: the marginals exist at the node. The saving is the
lookahead's solves of the children that would fathom anyway [e]. It can be measured offline on the
`line_signals` dumps and the truth sets before any solver change.
