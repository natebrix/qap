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

### 4. tai35b, the whole tree (`runs/tai35b_*`)

`infra/laptop/solve.sh tai35b` from the experiment build (coordinator local mode, 2 slots per arm, the
incumbent the optimum + 1, `tai/tai35b.par` = tai50b's production settings, rule 9 through depth 6 and
rule 2 below). The first three arms were stopped by the container restart at 06:45 UTC and resumed from
their stores. Nodes are the node trace's (one line per node, `qapctl export`): the coordinator's own
total varies with restarts and unit splits (972 and 923 for the same tree). **The production tree is
reproducible**: a second run, not interrupted (`tai35b_base2`), has the same 1,017 nodes at every depth;
its CPU is 7% higher on the loaded box, which is the CPU noise.

| arm | `QAP_LA` | nodes | CPU s | of which the lookahead | CPU s net of the lookahead | wall, 2 slots |
|---|---|---|---|---|---|---|
| production branching | - | 1,017 | 6,931 | - | 6,931 | 63 min |
| production branching, again | - | 1,017 | 7,390 | - | 7,390 | 61 min |
| rule 9's two lines, judged by the DNN, depths 0-6 | 1:30:0:6 | 1,046 | 8,682 | 1,672 (21 nodes) | 7,010 | 96 min |
| 4 rows and 4 columns, depths 0-2 | 4:30:0:2 | **780** | 8,900 | 2,803 (6 nodes) | **6,097** | 100 min |
| 4 rows and 4 columns, depth 2 only | 4:30:2:2 | 946 | 8,944 | 1,583 (3 nodes) | 7,361 | 78 min |

Every arm closed the tree at the optimum 283,315,445, found by the search.

- **Rule 9's own two lines judged by the DNN change nothing**: 21 branchings at depths 0-6 against
  production's 20, the same tree size. Both DNN arms branch the root on c26 where production takes r19
  (8 lines compared: c26's children 23M short in sum after 30 pooled iterations, r19's 43M, r18's 46M, the
  other five 140M or more).
- **Four lines a side at depths 0-2: 23% fewer nodes, 12-17% less CPU net of the lookahead**, all of it below
  depth 2: depths 3-5 branch 2, 2, 2 times against 4, 5, 2, and hold 225 nodes against 382. Its three
  depth-2 choices (c23, c23, c22) are not rule 9's lines.
- **At depth 2 alone it does little**: production's root and depth 1 (r19, then r18 three times), then
  the lookahead at production's three depth-2 nodes (c23, c6, c22, against rule 9's c5). 7% fewer nodes,
  all at depth 5 (93 against 155), and 7,361 CPU s net of the lookahead, within the noise of production's
  two runs. The three depth-2 nodes it branches are the depths 0-2 arm's three, the same pairs of fixes
  reached in the other order (r19 then r18, against c26 then r19 or r18); it picks the same line at two
  (c23, c22) and another at the third (c6 against c23), since the node's warm state depends on its
  parent, as the line study found for rule 9. The depths 0-2 arm's saving starts at the root: c26
  leaves 2 of its children to branch, where r19 leaves 3, so depth 2 holds 68 nodes against 102.
- **The lookahead costs more than it saves at this size**: 2,803 CPU s for 6 nodes (about 470 s a node:
  8 lines, about 260 pooled child solves of 30 iterations, 1.8 s each at m = 33), against 830-1,290 s saved,
  so +20-28% CPU in all. The cost scales with a node's children and the iteration cost, the saving with the
  subtree under the node: at tai60b's depth 2 the subtrees run to 1.7M CPU s and more, against about
  12,000 s for an 8-line lookahead there (464 children, 30 iterations at 0.86 s).

### 5. tai60b's 43:58, the lookahead over the line study's four lines (`dump.sh tai60b_la`)

The node in the line study's state (`../tai60b_lines/`), bounded as there, then the lookahead over r1, r45,
c26 and r8 (`QAP_LA_LINES`), 30 pooled iterations a child (`dumps/tai60b_la_cap30.log`):

| line | truth (rough subtree, `../tai60b_lines/`) | children open after 30 iterations | shortfall sum | rank |
|---|---|---|---|---|
| r1 | 19,600 (1.0x) | 47 of 58 | 236M | 3 |
| r45 | 22,900 (1.2x) | 53 | **162M** | 1 (chosen) |
| c26 | 38,800 (2.0x; rule 9's choice) | 58 | 191M | 2 |
| r8 | 67,100 (3.4x) | 45 | 236M | 4 |

- It picks r45, 1.2x the best line, where rule 9 from the same state picks c26 at 2.0x. Not r1.
- **30 iterations are too few at n = 60.** The full first solves of the line study fathom 52 of r1's
  58 children; after 30 warm pooled iterations 47 are still open, and the counts of open children
  rank the lines the wrong way round (r8 fewest). On tai35b (m = 33) 30 iterations already fathomed 21
  to 29 of a good line's 34 children. The cap has to grow with the instance.
- Cost: 6,364 CPU s for 229 pooled solves (28 s each, about 0.93 s an iteration), about 0.4% of the node's
  subtree in production (at least 1.67M CPU s). A cap of 100 to 150 iterations would cost 3-5x that and
  is still a small fraction.

### 6. tai40b, the whole tree (`runs/tai40b_*`)

The same set-up as tai35b (`tai/tai40b.par`, the incumbent the best known + 1, 2 slots):

| arm | `QAP_LA` | nodes | CPU s | of which the lookahead | CPU s net of the lookahead | wall, 2 slots |
|---|---|---|---|---|---|---|
| production branching | - | 731 | 7,379 | - | 7,379 | 89 min |
| 4 rows and 4 columns, depths 0-2 | 4:30:0:2 | 829 | 17,461 | 5,368 (6 nodes) | 12,093 | 213 min |

Both closed the tree at the best known 637,250,948, found by the search.

- **Production's tree is a path**: the root branches on r23, and at depths 1 to 19 one child of each
  node branches and the other 38 to 39 are fathomed. There is little to save.
- **The lookahead took the wrong root line**: c27, where production takes r23. At the root every child
  of all 8 candidate lines was still open after 30 pooled iterations (40 of 40), so the scores are the
  shortfalls of solves that have not converged: r23's children 419M short in sum, c27's 188M. Under c27,
  2 children branch at depth 1 and 3 at depth 2 (depths 2 and 3 hold 78 and 114 nodes, against 39 and 38),
  and those are the heavy nodes: 2.4 core-hours at depths 1-2 against 0.5, 1.2 of them the lookahead's.
- **+13% nodes, +64% CPU net of the lookahead, +137% in all.** As at tai60b's 43:58, 30 iterations are
  too few at m = 39: after the lookahead 37 to 40 of each candidate line's 38 to 40 children are open at
  every node it judged, except one depth-1 line (c16, 27 of 39).
