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
