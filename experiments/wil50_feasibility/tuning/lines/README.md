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

## wil50 result (laptop, Oct 31 2026)

Three depth-3 clear misses, 1,410 children: the xbar line (largest sum of
x_ij^2) is the best line at all three parents; eccentric is 2 to 4.7x worse,
lookahead 1.0 to 1.6x worse. Full tables and the line-score correlations are
in ../README.md ("Branching lines at wil50 clear misses"); log
`wil50_lines_m2pro.log`, child cache `work/`.

## Deep in the tail (Nov 6-7): does rule 9 still pick the best line at levels 4-5?

The study above scored lines at depth-3 parents only; the 329-405-short tail class
carries most of the estimated work, and its cost sits at levels 4-7. `parents_tail.txt`
(cloud) and `parents_tail2.txt` (laptop) list eight branched level-4 and level-5 nodes
of the node-1 and node-5 subtrees, each with the line the search itself branched on
(rule 9, third column, scored as `prod`). `--rules M,D,G` adds Fujii et al.'s rules and MatQapNB's rule G
(arXiv 2101.09629 section 5): the line whose children have the largest mean score,
the score being the child's objective at its barycentre (M) or the child's certified
lower bound from the parent's own certificate (D, our analogue of their dual rule)
or its Gilmore-Lawler bound (G, `matqapnb_rule_G.md`).
Production settings for the DNN runs (laptop, no AMX):

    DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 \
    DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 OMP_NUM_THREADS=1 \
    python3 experiments/wil50_feasibility/tuning/lines/lines.py qaplib/wil50.dat --ub 48816 \
        --parents experiments/wil50_feasibility/tuning/lines/parents_tail2.txt \
        --work experiments/wil50_feasibility/tuning/lines/work_tail2 --jobs 8 --random 1 \
        --rules M,D,G --dnn-flags "-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40" \
        > experiments/wil50_feasibility/tuning/lines/wil50_tail2_lines.log 2>&1

(the cloud adds `DNN_AMX=1 DNN_CHOLQR_PASSES=1` and runs `parents_tail.txt` into
`work_tail`).

### Results (Nov 7): rule 9 is within a few percent of the best line in the tail; branching is closed

Five parents scored [m]: `n1d4a` on the cloud (`work_tail/`, `wil50_tail_lines.log`), `n1d4c`,
`n1d4d`, `n1d5b`, `n5d4b` on the laptop (`work_tail2/`, `wil50_tail2_lines.log`). Each rule's line
cost against the best candidate line at that parent (1.00 = best; for the random baseline,
the mean of its two lines):

| rule | n1d4a | n1d4c | n1d4d | n1d5b | n5d4b | geometric mean | worst |
|---|---|---|---|---|---|---|---|
| production (rule 9) | 1.25 | 1.00 | 1.00 | 1.02 | 1.02 | **1.05** | 1.25 |
| better of the two xbar lines | 1.18 | 1.00 | 1.00 | 1.00 | 1.01 | 1.04 | 1.18 |
| M (barycentre objective) | 1.00 | 1.17 | 1.48 | 1.02 | 1.00 | 1.12 | 1.48 |
| G (Gilmore-Lawler) | 1.18 | 1.17 | 1.48 | 1.00 | 1.00 | 1.15 | 1.48 |
| QPB lookahead (rule 4) | 1.10 | 1.17 | 1.48 | 1.72 | 1.30 | 1.34 | 1.72 |
| eccentric (rule 6) | 2.77 | 2.86 | 3.29 | 2.90 | 1.62 | 2.61 | 3.29 |
| D (our certificate's child bounds) | 2.77 | 2.86 | 3.29 | 2.42 | 2.72 | 2.80 | 3.29 |
| random | 3.14 | 3.35 | 2.95 | 2.31 | 2.90 | 2.91 | 3.35 |

- **Rule 9 holds deep in the tail**: within 2% of the best line on four parents, 1.25x on
  one (`n1d4a`, where M's line was best), 1.05x on average. No rule beats it consistently:
  M and G win on two parents each and lose 1.2-1.5x on others; QPB lookahead is 1.1-1.7x
  worse; D and the eccentric rule are as bad as a random line. At depth 3 the xbar line was
  also the best (section above), so rule 9 is near the best line at every level measured.
- Choosing the better of rule 9's two xbar lines perfectly (its QPB lookahead picked the
  worse one on `n1d4a` and `n1d5b`) would give 1.04x against 1.05x: nothing to take.
- An oracle that always found the best of these lines would save about 5% per branching
  decision; over the three or four levels of a tail subtree that compounds to 1.15-1.2x at
  most, and no rule here comes close to the oracle. **Branching is closed for wil50.**
- D, the dual-bound rule, is consistently 2.4-3.3x off: no reason to spend more cycles on
  dual-bound branching.
- Not scored: the cloud's other three parents (`n1d4b`, `n1d5`, `n5d4`; on `n1d4b` rule 9, M, G
  and the xbar row all pick the same line). They were stopped at the cloud's background time
  limit and are not needed for the verdict; `work_tail/` resumes them if wanted.
- Caveats: five parents from the same deep region; the cost score is a model
  (47^(gap/130) per unfathomed child) and one far-short child can move a line's cost by
  10-20%, so the ranking among the top three lines is not resolved, only that nothing
  beats rule 9 by a margin that matters.

## Open (Nov 7): the root and depth-1 lines are untested, and they set the tail

Everything above is depth 3 and below. No line study has run at depths 0-2, and
that is where the remaining leverage is, for two reasons that do not hold deeper
down (`../../roadmap/tree_size.md` sections 2-3; nothing here was run):

- **What the top lines decide.** Nothing fathoms at depths 1-2 whatever the line,
  so the depth-3 frontier is 29,580 orbits regardless. The lines at depths 0-2
  decide the *distribution* of depth-3 gaps, and the two tail nodes among the 16
  sampled (329 and 405 short) carry 85-91% of every tree estimate. In the
  exponential model (subtree ∝ 47^(gap/130)) a 10% change in a tail node's gap is
  2.7-3.3x on its subtree. The 5% per decision left on the table at depths 4-5 is
  a different order of question.
- **Why rule 9 is unvalidated there.** The root is about 680 short and the
  relaxation is far from integral, so the marginal row and column it reads are
  diffuse; the parents where rule 9 was validated were 160-240 short. With
  `wil50t.sym` the root must branch on a facility (15 orbit children whatever the
  facility), so the root decision is only *which* facility; the same holds at
  depth 1 whenever the placed cell is in the grid's middle row.

**The test, when there is time** (not run; the budget says it is cheap):

- Root: all 50 facilities, 15 children each at about 250 s, about 50 CPU-hours
  (under 0.03% of the 6e8 CPU-s budget); score each facility by the lines cost
  (sum over its children of 47^(gap/130)) and compare rule 9's pick against the
  best. 750 independent depth-1 solves: an overnight job on the laptop at 8 jobs.
- Depth 1: 15 nodes; rule 9's two xbar lines plus M's line and one random line
  per node, about 2,700 child solves, about 2e5 CPU-s, a day on the cloud box.
  Exhaustive strong branching (49 lines x 15 nodes) is not needed.
- Both fit inside the sampling plan's stage A, which runs depths 0-2 under the
  production rule anyway.
- Caveat: child gap at depth 1 is a proxy for subtree size two levels down. In the
  tail study the same proxy's rank correlation with realised cost was +0.6 to +0.9,
  enough to see a factor of two, not a 5% difference.
