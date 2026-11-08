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

- **What the top lines decide.** *(Corrected Nov 8: with the CP cut phase at
  depth 2, stage A fathoms about half of the depth-2 nodes, so the frontier is
  about 15,100 depth-3 nodes, not 29,580; the fathomed ones were the near-cutoff,
  cheap ones.)* The lines at depths 0-2 decide the *distribution* of depth-3
  gaps and which depth-2 nodes fathom, and the two tail nodes among the 16
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

### First root result (Nov 7, laptop): rule 8's root pick is 5x worse than its runner-up

`root/`: the root solved 1,000 iterations with `-X` (`root_xbar.log`, lb 48,122.2,
691.9 short), then the 15 orbit children (J1 cells: grid rows 0-2 x columns 0-4) of
rule 8's three top facilities, cold, 1,500 iterations each, production flags, no
give-up (`run_children.sh`, 45 logs in `children/`, scored by `score.py`; 52 min
wall at 8 jobs, about 440 s a child). On `wil50.dat` a child is `-x cell:facility`.

| facility | rule 8 score (sum x^2) | flow total (rank of 50) | children fathomed | median gap | worst gap | cost = sum 47^(gap/130) |
|---|---|---|---|---|---|---|
| **45** | 0.125 (2nd) | 178 (50th) | 0 | 390 | 515 | **1.27e7** |
| 26 | 0.098 (3rd) | 189 (47th) | 0 | 494 | 552 | 6.33e7 (5.0x) |
| 25 | 0.134 (1st, **rule 8's pick**) | 184 (48th) | 0 | 495 | 559 | 6.72e7 (5.3x) |

- Rule 8's root choice (facility 25, by a 7% margin in sum x^2 over 45) costs
  5.3x the runner-up on the lines metric. The root marginals are flat (no entry
  above 0.18 against a uniform 0.02), and rule 8's ranking among the top lines
  does not predict their cost there. The question "is our root line the best
  possible" is therefore open in the direction that matters.
- The pattern: the three most concentrated facilities are among the four
  lowest-flow ones, and their mass sits on the four corner cells. The corner
  child (cell 0) is the *hardest* child for all three (gap 515-559); the
  easiest is the centre cell 24 (288-452). Facility 45's children gain 177-404
  on the root's 692, facility 25's 133-240.
- Not yet run: the two highest-flow facilities (40, flow 267; 31, 259) and a
  median facility as control, which test the opposite hypothesis (diffuse
  marginals but every child pins down something that matters). Same script:
  `run_children.sh 40 31 <median>`, about 50 minutes at 8 jobs.
- Caveats: cold children at 1,500 iterations (all `iter_limit`), so gaps are a
  little pessimistic for everyone alike; the lines cost is the model of
  `../README.md`, and 5x is far outside its noise.

**Facility 40, the highest-flow facility (flow 267), Nov 7, 18 min**: 15 children, 0
fathomed, gaps 190 to 643, median 504, cost 4.01e+08 (32x facility 45,
6.0x rule 8's pick). Its corner child (cell 0, gap 190) is the easiest child
seen at the root, but its centre children barely move the bound (cells 14, 23, 24:
619, 607, 643 on the root's 692), and the lines cost is set by the worst
children. So the opposite hypothesis fails too: a diffuse high-flow line is worse
than a concentrated low-flow one, by the spread of its children. The best root
line so far is facility 45, found by rule 8's score as runner-up, and nothing
here says how to find it a priori. Facility 31 and a median control not run.

**Facility 31, the second-highest-flow facility (flow 259), Nov 7, 18 min**: 15 children,
0 fathomed, gaps 302 to 570, median 506, cost 1.03e+08 (8.1x facility 45,
1.53x rule 8's pick, 0.26x facility 40). Corner child 302, centre children
(cells 14, 23, 24) 564, 566, 570. Median control not run.

**Facility 21, rule 8's fourth (flow 190, rank 46th), Nov 7, 18 min**: 15 children, 0
fathomed, gaps 276 to 513, median 403, cost 1.48e+07 (1.2x facility 45). Corner child
508, centre children (cells 14, 23, 24) 302, 320, 276.

### Root summary (Nov 7): six facilities, 90 children, about 2 hours at 8 jobs

| facility | why tested | rule 8 score | flow (rank) | median gap | worst gap | cost | vs best |
|---|---|---|---|---|---|---|---|
| 45 | rule 8's 2nd | 0.125 | 178 (50th) | 391 | 515 | 1.27e+07 | 1.0x |
| 21 | rule 8's 4th | 0.086 | 190 (46th) | 403 | 513 | 1.48e+07 | 1.2x |
| 26 | rule 8's 3rd | 0.098 | 189 (47th) | 494 | 552 | 6.33e+07 | 5.0x |
| 25 | rule 8's pick | 0.134 | 184 (48th) | 495 | 559 | 6.72e+07 | 5.3x |
| 31 | 2nd-highest flow | - | 259 (2nd) | 506 | 570 | 1.03e+08 | 8.1x |
| 40 | highest flow | 0.083 | 267 (1st) | 504 | 643 | 4.01e+08 | 31.5x |

- The root line matters by 5x or more between lines rule 8 ranks within 7% of each
  other, and the best line (45) is found by rule 8's score only as runner-up.
  Flow total is the wrong direction (the two heaviest facilities are the two worst
  lines). What distinguishes 45 is the narrow spread of its children; no root
  statistic tried here predicts that, child bounds do.
- Proposed fix (not built): at the root, and perhaps at depth 1, shortlist by rule 8
  (top 4-5 lines) and choose by the children's bounds (short cold solves, a few
  hundred iterations each; the full 1,500-iteration version above cost about 20
  CPU-hours per facility-set of six, once per run). Whether depth 1 has the same
  problem is the next test and is the cloud box's size (15 nodes x 4 lines x 45
  children).
- Caveats: cold children, 1,500 iterations, all `iter_limit`; the lines cost model
  of `../README.md`; one level measured, two tail nodes dominate the tree estimate
  and their lines are set at depths 0-2 together.
- Stage A (cloud, Nov 7-8, `sampling/`) happened to branch the root on facility 21
  (its root has only the QPB bound, so the rule there is not rule 8; rule 8's pick of
  25 arises only with DNN marginals at the root). Its depth-1 children at 300
  iterations rank as the 1,500-iteration solves here do (corner hardest, centre
  easiest), 333-570 short. The cloud is adding `--root-facility` to `stage_a.py`
  and running 45 and 25 on the 15 depth-1 nodes plus 200 random depth-2 nodes each,
  so the three root lines will be compared on the real depth-2 bounds, with the
  depth-1 line question answered by the same data.

### Calibration (Nov 8, laptop): 300 cold iterations rank the root lines as 1,500 do

The same 90 children (`children_i300/`, `children_i100/`), the same command with
`-i 300` and `-i 100`, scored with `score.py` against the 1,500-iteration ranking.
The order of the six facilities is identical at every iteration count; only the
cost scale moves, since the cold bound is still rising at 300 iterations (median
gap 457 against 391) and far from done at 100 (gap about 2,500, where the cost
model's numbers mean nothing but the order survives).

| facility | 1,500 iters: cost, vs best | 300 iters: cost, vs best | 100 iters: vs best | CPU per facility (1,500 / 300 / 100) |
|---|---|---|---|---|
| 45 | 1.27e+07, 1.0x | 6.93e+07, 1.0x | 1.0x | 6,637 / 1,277 / 327 s |
| 21 | 1.48e+07, 1.2x | 7.81e+07, 1.1x | 1.1x | 6,289 / 1,258 / 321 s |
| 26 | 6.33e+07, 5.0x | 3.01e+08, 4.3x | 1.8x | 5,976 / 1,275 / 322 s |
| 25 | 6.72e+07, 5.3x | 3.06e+08, 4.4x | 8.8x | 6,774 / 1,278 / 317 s |
| 31 | 1.03e+08, 8.1x | 4.85e+08, 7.0x | 14x | 6,331 / 1,274 / 317 s |
| 40 | 4.01e+08, 31.5x | 1.95e+09, 28x | 210x | 6,531 / 1,206 / 301 s |

So the proposed fix (shortlist by rule 8, choose by child bounds) can be run at
300 iterations for a fifth of the 1,500-iteration price: about 85 CPU s per
child, 1,300 CPU s per candidate line, 20 CPU-minutes for a shortlist of
rule 8's top 5 at one node of a 1,000-core run. The 100-iteration pass is a
further 4x cheaper and still orders these six correctly, but the two best lines
are only 7% apart there and the cost scale is meaningless, so 300 is the
setting to trust. One level, one instance, six lines; the cloud's depth-2
comparison (above) is the real test.

## Depth 1 (Nov 8, laptop): rule 8's pick is 2.4x the best line at node (cell 0, facility 45)

`depth1/`: the depth-1 node with facility 45 in cell 0 (the best root line's
corner child) solved 1,000 iterations with `-X` (`node_0_45_xbar.log`, lb
48,297.7, 165 s), rule 8's score (sum of squared marginals over a free line)
computed for all 98 free lines, cells and facilities alike (`shortlist.txt`), the
top 4 lines taken, and every child of each (49 per line, `-x 0:45,c:f`) solved
300 cold iterations with the Nov 7 child settings (`children/c{c}_f{f}.log`,
196 children, 94 min at 8 jobs, about 69 CPU s each). Cost is the lines model,
sum of 47^(gap/130); 300-iteration costs rank lines as 1,500-iteration ones do
(root calibration above) but are not comparable in scale to the root tables.

| line | rule 8 rank, score | median gap | min / max gap | cost | vs best |
|---|---|---|---|---|---|
| facility 21 | 2nd, 0.169 | 203 | 66 / 468 | 4.11e+06 | 1.0x |
| cell 40 | 4th, 0.149 | 311 | 98 / 488 | 4.52e+06 | 1.1x |
| facility 25 | 3rd, 0.166 | 362 | 319 / 487 | 9.23e+06 | 2.2x |
| cell 49 | **1st, 0.193** | 298 | 75 / 490 | 9.81e+06 | 2.4x |

- The root's pattern repeats one level down, smaller: rule 8's top-scored line
  is the worst of its own top 4 by child bounds, 2.4x the best (5.3x at the
  root). The spread between the four is 2.4x against 5x at the root, as expected
  if the line choice matters most where the subtrees are largest.
- Facility 21 is the best line here as it was the root's runner-up (1.2x there);
  facility 25, rule 8's root pick, is again a poor line. Cell lines are in the
  mix at depth 1 (cell 40 within 10% of the best), where the root study tested
  facilities only.
- Median gap and cost disagree (cell 40's median is 100 above facility 21's,
  its cost 10% above): the cost is set by the few hardest children, so a line
  with a narrow spread of children beats one with a better typical child, as at
  the root.
- Four lines are the shortlist, not the full 98; the best line at this node may
  be outside rule 8's top 4. Cost of the test: 196 children x 69 CPU s = 3.8
  CPU-hours, once per depth-1 node. At 300 iterations the shortlist-then-child-
  bounds rule costs about 1 CPU-hour per candidate line per node; at the 15
  depth-1 nodes of a wil50 run that is a rounding error against the tree
  (NOTES.md: 1e10 CPU-s), and the 2.4x here, if it holds, is worth more than
  any switch measured since the cut phase.
- Next: the cloud's stage-A comparison of root lines 21, 45 and 25 on real
  depth-2 bounds, and the same shortlist test on a second depth-1 node (an
  interior cell) to see whether 2.4x is typical.
