# Branching lines at wil50 clear misses (set up Oct 31 2026; to run on the laptop)

Question: at the depth-3 nodes the DNN bound misses by a lot (about 95% of the
estimated wil50 tree is below them, ../../NOTES.md section 4), which line
(facility = row, location = column) leaves the smallest subtree, and do the
rules we have pick it: the most eccentric location (rule 6), QPB lookahead
(rule 4), or a new candidate, the DNN primal's own marginals x_ij?

`lines.py` (see its docstring) bounds, for each parent of `parents.txt`, every
child of about 9 candidate lines (eccentric; the top 2 rows and 2 columns by U
sum, rule 4's candidates, with lookahead's own pick among them marked; the row
and column whose x_ij are most concentrated; 2 random rows and 2 random
columns), with a fixed DNN budget per child, and scores each line by
sum over children of 47 ** (gap / 130), the size of the subtree a child
leaves if the bound gains about 130 per level (wil50, ../../NOTES.md). It also
reports the rank correlation of x_ij (known before branching) and of the
child's QPB with the child's gap.

The parent run uses `dnntest -X`, which prints the marginals x_ij of the last
primal iterate (new, Oct 31).

## Running it

    VECLIB_MAXIMUM_THREADS=1 python3 experiments/wil50_feasibility/tuning/lines/lines.py \
        qaplib/wil50.dat --ub 48816 \
        --parents experiments/wil50_feasibility/tuning/lines/parents.txt \
        --work experiments/wil50_feasibility/tuning/lines/work --jobs 8

`--dry-run` lists the candidate lines (after the parent runs) without bounding
children. Cost: about 9 lines x 47 children x 3 parents = 1,270 child bounds
at 300 iterations each (about 150 s each on the 2.1 GHz cloud Xeon), plus the
child QPB bounds through the console: about 50 CPU-hours there, so roughly
2-3 hours on the M2 Pro with 8 workers. It caches every child in `--work`
and resumes after an interruption.

## Check on nug30 (cloud)

The nug30 depth-2 far miss `0:13,6:14` (75.7 short; the node of
experiments/dnn_bound/predictor/lines), `--child-iters 300 --gain 40
--branch 27`, 11 lines, 308 children (`nug30_check.txt`):

| line | selected by | children fathomed | cost |
|---|---|---|---|
| row 4 | xbar (most concentrated row) | 20/28 | 111 (best) |
| col 0 | ecc | 18/28 | 203 |
| col 5, col 29 | U (col 29: lookahead's pick) | 18/28 | 238 |
| random lines | | 11-16/28 | 221-420 |
| col 15 | xbar (most concentrated column) | 13/28 | 429 (worst) |

- **x_ij predicts a child's gap as well as the child's QPB does**, and it is
  free once the node's DNN bound has run: pooled over the 308 children the
  rank correlation with the gap is +0.745 for x_ij (larger x, harder child)
  and -0.750 for the child's QPB.
- **The line score matters.** Ranked by how well they order the 11 lines'
  costs (Spearman): children with x > 0.05 (fewer is better) +0.70; with
  x > 0.02 +0.51; sum of x^2 (the "xbar" selector above) +0.39; mean child
  QPB (lookahead's criterion, which costs a QPB solve per child) +0.28;
  children with x > 0.1 +0.25. The "most concentrated" selector found the
  best line as a row and the worst as a column.
- Eccentric and lookahead beat random (0.62x and 0.72x the random lines'
  mean cost) but were 1.8-2.2x worse than the best line.

One node and 11 lines: a hint, not a result. The wil50 run reports the same
score correlations per parent; if "children with x > 0.05" holds up there, a
branching rule that picks the line with the fewest such children would cost
nothing per node (the DNN bound runs before branching) and replace lookahead
at the DNN depths.
