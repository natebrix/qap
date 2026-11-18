# Which signal picks the branching line? (Nov 18 2026, running)

Question (Nathan, after the tai60b line study, `../tai60b_lines/`): rule 2 uses QPB's dual matrix U
(how many children bound + U_ij fathoms, ties by the U sum); is there an equivalent that uses the DNN
bound? The tai50b and tai60b runs suggest there are savings in the line choice.

What was known: the DNN certificate's own child bounds (`DNNResult::child_lb`, lb + the reduced costs of
its outer assignment problem) are flat at convergence, by complementary slackness, and were measured
useless (`../wil50_feasibility/roadmap/certificate/REPORT.md`, item 3). The children's first DNN solve
ranks lines by their subtrees (Spearman 0.71-0.93, `../tho30_solve/line_truth/`,
`../tho40_line_truth/`); the marginal score does not.

## Data

The truth: the line-truth studies' 6 nodes (tho30 depth 1 x 3, tho40 depth 3 x 3), 8 candidate lines
each, every child's subtree solved to the end (b), and every child's cold first DNN solve, give-up
300 (a). And tai60b's 43:58 in the line study's state, 4 lines, ranked by its rough subtree model.

## Tools

- `experiment_build.patch`: the experiment build (a worktree of cognac 8df9573; never cognac's
  `latest/`): `../tai60b_lines/force_line.patch` plus
  - `QAP_LINE_SCORES=D` (and `QAP_LINE_SCORES_LA=1`): at a node of depth D, dump what the branching rules
    read: QPB's U, the bound, the cutoff, the DNN marginals, and rule 4's lookahead score of every line
    (each child bounded by QPB, on a copy of U);
  - `QAP_LA=K:CAP:DMIN:DMAX`: a DNN lookahead rule. At a node of depth DMIN-DMAX whose DNN bound missed,
    the candidates are the K rows and K columns with the largest marginal scores (K = 1: rule 9's two
    lines); every child of a candidate that U does not drop is solved warm from the node's state (the
    children's own start: Z and the conditioned Y) for CAP ADMM iterations, no give-up, no pool, no cut
    phase; the node branches on the line whose children's shortfall sum (cutoff - lb where positive) is
    least. `dnn-la:` lines on stderr give its choice, the candidates' scores and its CPU.
- `dump.sh`: the dumps at the truth nodes (`dumps/`); `signals.py` scores the lines by each signal
  against the truth; `warm.py run|analyze`: every child of the study's lines warm from the node's
  (uncut) state, 200 iterations, the bound every 10 (so any cap up to 200).
- `tai/`: tai35b and tai40b recipes, `infra/phase1/recipes/tai50b.par` with the instance's DNN depths and
  the bound lab's region files (`../bound_lab/tai_solves/`: comb for tai35b, combh for tai40b).

## Results

### 1. The cheap signals do not rank lines (`signals.py`, `signals_out.md`)

At the 7 nodes, Spearman of each signal against the truth (positive: agrees), whether the truth's best
line is the signal's first, and the truth nodes of the signal's first line over the best line's:

| signal | mean Spearman | best line 1st | in top 2 | in top 3 | its top line / best, mean (worst) |
|---|---|---|---|---|---|
| marginal score (rules 8, 9) | -0.07 | 1/7 | 1/7 | 2/7 | 1.52x (1.99x) |
| rule 2 (fathom count, then U sum) | +0.06 | 1/7 | 3/7 | 3/7 | 2.44x (4.71x) |
| U sum (rule 1) | +0.06 | 1/7 | 3/7 | 3/7 | 2.44x (4.71x) |
| QPB lookahead over every line (rule 4) | +0.34 | 1/7 | 4/7 | 5/7 | 1.67x (3.42x) |
| marginals near zero | -0.25 | 1/7 | 2/7 | 3/7 | 2.63x (8.01x) |
| (a) the children's first DNN solve | +0.80 | 5/6 | 5/6 | 6/6 | 1.00x (1.00x) |

- **Rule 2's count is 0 at every node.** QPB is 13% under the cutoff at tho30's and tho40's nodes, 29%
  at tai60b's 43:58 (431,951,821 against 608,215,054), so no child is fathomed by bound + U and rule 2 is
  its tie-break, the U sum, which says nothing (+0.06). That is QPB, not rule 2's idea: the DNN
  equivalent of U is flat (above), so the idea has no DNN form either.
- **tai60b 43:58**: every cheap signal ranks the four lines -0.40 (they favour r8, the worst). Among
  all 116 free lines, production's r1 is 13th by marginal score (the 6th row), 22nd by the QPB lookahead
  (14th row), 33rd by the U sum.

### 2. The DNN children decide, given a wide enough shortlist

Rule 9's shortlist widened to the K best rows and K best columns by marginal score, the line chosen by
(a): truth nodes of the chosen line over the best of the study's lines (6 nodes with full truth):

| node | rule 9 (QPB arbiter, K = 1) | (a), K = 1 | K = 2 | K = 3 | K = 4 | all the study's 8 lines |
|---|---|---|---|---|---|---|
| tho30 5:1 | 1.24x | 1.24x | 1.24x | 1.00x | 1.00x | 1.00x |
| tho30 5:2 | 1.47x | 1.47x | 1.47x | 1.00x | 1.00x | 1.00x |
| tho30 5:3 | 1.00x | 1.00x | 1.28x | 1.28x | 1.00x | 1.00x |
| tho40 24:7,26:1,39:0 | 1.42x | 1.42x | 1.84x | 1.00x | 1.00x | 1.00x |
| tho40 25:32,26:8,39:0 | 1.65x | 1.33x | 1.33x | 1.33x | 1.00x | 1.00x |
| tho40 26:3,31:7,37:0 | 1.65x | 1.00x | 1.00x | 1.00x | 1.00x | 1.00x |
| mean | 1.41x | 1.24x | 1.36x | 1.10x | 1.00x | 1.00x |

Rule 9 as it is leaves 1.41x on these nodes; its own two lines judged by the DNN, 1.24x; 3 to 4 lines a
side recover the best line. At tai60b's 43:58 the best of the four (r1) is the 6th row by marginal
score: there it takes K = 6.

### 3. A short warm solve is enough (`warm.py`, `warm_out.txt`; tho30's 3 nodes, all 8 lines)

Every child warm from the node's uncut state (Z and the conditioned Y), no give-up. Child by child, the
shortfall at 30 iterations correlates 0.92-0.98 with the cold 300-iteration (a) and 0.80-0.90 with the
child's subtree. By line:

| cap (iterations) | Spearman with the truth, 5:1 / 5:2 / 5:3 | top line / best | CPU s per line (29 children, m = 28) |
|---|---|---|---|
| 10 | +0.33 / +0.14 / +0.42 | 1.24x / 1.03x / 4.71x | 11 |
| 30 | +0.57 / +0.50 / +0.78 | 1.03x / 1.03x / 1.00x | 29 |
| 50 | +0.67 / +0.55 / +0.87 | 1.03x / 1.03x / 1.00x | 43 |
| 100 | +0.71 / +0.79 / +0.97 | 1.00x / 1.00x / 1.00x | 71 |
| 200 | +0.74 / +0.79 / +0.97 | 1.00x / 1.00x / 1.00x | 124 |
| cold 300, (a) | +0.71 / +0.71 / +0.93 | 1.00x / 1.00x / 1.00x | about 90 |

30 warm iterations pick within 3% of the best line at a third of (a)'s cost; 100 match (a) and better
it. A child's setup (allocation, the first full eigendecomposition) is a fixed part of the cost: about
7 CPU s per child at m = 34 with 30 iterations (tai35b's root, below).
