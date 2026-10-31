# wil50: cheap CPU levers before any cluster run (started Oct 30 2026, M2 Pro)

Two experiments from the discussion after the 100-core estimate
(../NOTES.md section 9). Both were cut short when the laptop session ended;
scripts and partial results are here so the cloud session can finish them.
All runs use the solver's mode (subspace tracker) with the Peaceman-Rachford
step and degree-4 filter: `dnntest -T 1 -d 4 -g 1.0 -G 1.0`.

## 1. Eliminated assignments as a bound-strength lever (`elim/`)

Aim: the DNN bound with eliminated assignments (DNNParams::elim, commit
9862021) is a strictly stronger relaxation of a node. If cheap eliminations
exist at wil50 depth 3 to 5, they move fathoming up and divide the tree.

Learned so far:

- **QPB's U matrix eliminates nothing at wil50 depth 3.** Console (`b`, `u`)
  on the 5 near-miss and 2 fathomed depth-3 sample nodes: parametric QPB
  46,928 to 47,081 against the cutoff 48,814 (gap about 1,870), largest U
  entry 222 to 282. Zero cells of 2,209 satisfy bound + U >= cutoff on any
  node (`elim/console_out.txt`, `elim/elims.txt`). The QPB route is closed.
- **The DNN certificate's own child bounds are the candidate instead.**
  `dnntest -C 2` prints, from the certificate's outer assignment problem,
  a certified bound for every child (lb + reduced cost), with no extra SDP
  solve. Children above the cutoff are exactly the assignments that can be
  eliminated. `elim/round.sh <id> <spec>` bounds a node with `-C 2`, forbids
  the children the certificate fathoms (`-E`), bounds again, and reports both
  bounds and both fathomed-children counts. It was started on the 7 nodes
  (`elim/order.txt`) with `-H 0 -i 2500` and killed before the first round
  converged; `elim/r1_*.log` are partial.

Open: how much does round 2 gain on the near misses (59 to 116 short), how
many more children fall, and does it iterate? A gain of 60 on a near miss is
half a level; the division-by-47 table in ../NOTES.md section 5 prices that.

## 2. Give-up threshold under PRSM (`giveup/`)

Aim: clear misses are most of the nodes and each costs the full give-up
budget (`-H 300`). PRSM's bound at 300 iterations is 60 better than ADMM's,
so the hopeless rule might fire at 150 or 200 with the same decisions, a
1.3 to 1.5x cut on the dominant node class.

Run: all 32 sample nodes (`giveup/all32.txt`: 16 depth-3, 16 depth-4
children of nodes 9 and 6) with `-H 150`, `-H 200`, `-H 300`, 8 at a time
(so CPU per node is about 2x the alone figures). Results in
`giveup/H150_*.txt`, `giveup/H200_*.txt`; the H=300 pass and 6 of the H=200
nodes did not run.

Learned so far, and it matters:

- **H=150 is unsafe.** It abandons 9 nodes that fathom: depth-3 node 12
  (given up at 170 iterations, lb 48,703; fathoms at 280 iterations with
  H=200) and depth-4 nodes 17, 18, 22, 24, 27, 29, 31, 32 (given up at 150
  with lb 48,644 to 48,732; 17 and 18 fathom at 210 and 270 iterations with
  H=200, and all 16 depth-4 children fathomed in the original -f/ADMM sample
  with H=300). A false give-up costs n-d children, about 46 nodes, so one
  in four is ruinous.
- **H=200 is unproven.** It fathoms everything that finished (node 12, and
  depth-4 nodes 17 to 21, 23, 25, 26, 28, 30), but the six slowest depth-4
  fathoms (22, 24, 27, 29, 31, 32) were still running when the session ended.
- The hopeless extrapolation was tuned on ADMM trajectories. Under PRSM the
  certified bound at 150 iterations sits 100 to 170 below its final value on
  nodes that fathom, so the rule's "recent gains" extrapolation is too
  pessimistic early. Either keep H=300 with PRSM, or re-fit the extrapolation.

## Suggested completion (cloud)

1. Finish `giveup/` at H=200 and H=300 for the 32 nodes (alone or 4 at a
   time); report false give-ups and CPU per class. If H=200 has any false
   give-up, look at the hopeless rule's extrapolation under PRSM rather
   than at H.
2. Run `elim/round.sh` on the 7 nodes to convergence; report round-1 and
   round-2 bounds, children fathomed by the certificate in each round, and
   whether a third round still gains. If it gains, the next step is wiring
   `child_lb`-based eliminations into bnb.C's DNN call (children the
   certificate fathoms need no bound at all, and the surviving children can
   inherit the parent's elimination set).
3. Only then the line-comparison experiment on wil50 depth-2 nodes
   (../NOTES.md discussion; port of experiments/dnn_bound/predictor/lines).


# Results (cloud, Oct 31 2026: Xeon 2.1 GHz, 4 cores; solver mode with PRSM and degree 4)

All runs `dnntest -T 1 -d 4 -g 1.0 -G 1.0 ... -u 48816` (cutoff 48,814.05),
at most 4 at a time. CPU figures are for this machine and are not comparable
with the laptop's.

## 1. Give-up threshold (`giveup/cloud/`, `run_cloud.sh`)

All 32 nodes of `all32.txt` (id = line number) with `-H 300` and `-H 200`,
`-i 2500`. Logs `cloud/H<H>_<id>.log`.

| | H = 300 | H = 200 |
|---|---|---|
| fathomed | 18 (depth-3 nodes 11, 12; all 16 depth-4 nodes) | the same 18, same iterations and bounds |
| CPU per fathomed node (mean) | 105 s | 105 s |
| given up | 14 (all other depth-3 nodes) | the same 14 |
| CPU per given-up node (mean) | 175 s | 95 s |
| certified lb at the give-up point | 48,359-48,752 (62-455 short; the near misses 6, 9, 14, 15, 16 at 360-480 iterations, 62-130 short) | 48,279-48,648 (166-535 short, all at 200 iterations) |
| false give-ups | - | **0** |
| total CPU, 32 nodes | 4,329 s | 3,205 s (-26%) |

**H = 200 makes no false give-up on this sample, but its margin is thin.**
`-v` traces of the four slowest fathoms (nodes 22, 24, 27, 29: 340-400
iterations; `cloud/trace/`) replayed through the rule (`replay.py`): it would
fire at iterations 130-180 on all four, so they are safe only from H = 180-190,
10-20 iterations below 200. And at H = 200 the near misses happen to give up
exactly at 200; at H = 220 or 250 they run on to 230-440 (`rules.py`). Varying
the rule's factor (5 = current, 10, 20, 40) does not help: anything that
protects the fathoms earlier than 180 stops giving up on near misses at all.

**Why: under PRSM every node's bound is still rising steeply before about 180
iterations** (the far misses gain about 1,300 between iterations 100 and 200),
so an extrapolation of recent gains cannot tell fathoms from misses that early.

**The gap to the cutoff separates them.** The gap of the running best shrinks
by a factor of 0.34-0.55 every 50 iterations on all four fathoms, right up to
the fathom; on the misses the 50-iteration ratio is 0.58-0.77 at iteration 200
and rises steadily toward 0.9-0.98 (`gaprule.py`). A rule "from iteration H,
give up when g(t)/g(t-50) > rho" with H = 150, rho = 0.65 never fires on the
fathoms (worst ratio 0.55) and gives up on the four traced misses at 190-230.
It is in the solver as an opt-in (`DNNParams::hopeless_gap_ratio`, `dnntest
-R rho`; off by default), and was run on the 32 nodes (`run_cloud_R.sh`,
`cloud/H150R065_<id>.log`):

| 32 nodes | H = 300 | H = 200 | H = 150, rho = 0.65 |
|---|---|---|---|
| fathomed | 18, mean 105 s | the same 18, 105 s | the same 18, 106 s |
| given up | 14, mean 175 s, at 300-480 iterations | 14, mean 95 s, all at 200 | 14, mean 98 s, at 170-230 |
| false give-ups | - | 0 | 0 |
| total CPU | 4,329 s | 3,205 s | 3,284 s |

The gap-ratio rule saves what H = 200 saves (-24% against H = 300), with no
false give-up, and its margin does not depend on when the check happens to
fall: on the traced fathoms the ratio never exceeds 0.55 against the threshold
0.65, and on misses it only rises. Clear misses go at 170-190, near misses at
210-230 (nodes 6, 9, 15, 16), where they are 134-170 short.

**Recommendation.** Do not go below H = 200 with the extrapolation rule under
PRSM. For runs with PRSM, H = 150 with the gap-ratio rule (`-R 0.65`; in the
solver, `DNNParams::hopeless_gap_ratio`, not yet exposed in bnbdrv's `-x`) is
the safer way to get the same 24-26%. Both rest on 32 nodes of one instance
(and four `-v` traces); a nug30 A/B on the laptop would check them where the
answer is known. On the dominant node class this is a constant factor; it
does not change the depth at which fathoming starts (../NOTES.md section 5).

## 2. Eliminations from the DNN certificate (`elim/`, `run_cloud.sh`)

`round.sh` on the 7 nodes of `order.txt`, `-H 0 -i 2500 -C 2`
(`r1_<id>.log`, `summary_<id>.txt`):

| node | round-1 certified lb | short of the cutoff | children the certificate fathoms | round 2 |
|---|---|---|---|---|
| 6 | 48,758.2 | 55.8 | 0 of 2,209 | nothing to eliminate |
| 9 | 48,770.6 | 43.4 | 0 of 2,209 | nothing to eliminate |
| 14 | 48,730.5 | 83.5 | 0 of 2,209 | nothing to eliminate |
| 15 | 48,753.7 | 60.4 | 0 of 2,209 | nothing to eliminate |
| 16 | 48,758.8 | 55.3 | 0 of 2,209 | nothing to eliminate |
| 11 | 48,815.4 | fathomed (290 iterations) | - | - |
| 12 | 48,814.9 | fathomed (280 iterations) | - | - |

The five near misses run to the 2,500-iteration limit with the primal and the
certified value within 2 of each other, so these are the converged DNN values:
43 to 84 short (the earlier ADMM sample had them 59 to 116 short at 720-800
iterations; PRSM to convergence gains 15-30). The certificate's child bounds
(lb + reduced cost of the outer assignment problem) are equal to lb for 90%
of the children; the largest gain is 2.8 to 6.5, against the 43 to 84 a child
needs; no child is within 38 of the cutoff. At a converged multiplier the
outer assignment problem is almost flat, as on nug30
(`experiments/dnn_bound/predictor`): the certificate cannot eliminate
assignments, and round 2 never ran. The answer to the question is no: this
adds nothing to the near misses.

## 3. Wiring child_lb eliminations into bnb.C

Not worth doing, given part 2: no child of any of these nodes is eliminated.
For the record, it would be: pass `child_bounds = 1` in bnb.C's DNN call
when the node is not fathomed; drop the children whose `child_lb` exceeds the
cutoff (no bound needed); and give each surviving child the parent's
elimination set (stored on the node, as the warm start's Z is) plus its own,
passed to the DNN call as `DNNParams::elim` (commit 9862021). The machinery
exists; the eliminations do not.

# Results (laptop, Oct 31 2026: Apple M2 Pro, Accelerate, 8 workers)

## nug30 A/B of the give-up rules (`giveup/nug30_ab/`)

`nug_deep.par --ub 6124`, PRSM and degree 4 (`QAP_DNN_GAMMA1=1.0
QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4`), same binary, runs back to back. The
gap-ratio rule reaches the solver through `QAP_DNN_GAP_RATIO` (commit d106a87).
Baseline is the Oct 30 run in experiments/eigensolver/README.md section 8.

| | `-x 2:7:1:300:1` (Oct 30) | `-x 2:7:1:200:1` | `-x 2:7:1:150:1`, gap ratio 0.65 |
|---|---|---|---|
| nodes | 4,555 | 4,529 | 1,040,545 |
| DNN nodes / fathomed | 4,545 / 4,381 | 4,519 / 4,356 | 6,367 / 6,132 |
| DNN CPU | 5,307 s | 4,876 s | 4,457 s |
| CPU | 5,334 s | 4,904 s | 4,634 s |
| wall | 13 min 21 s | 11 min 17 s | 10 min 26 s |

- **H = 200 is clean on nug30**: 8% less CPU, no new misses (DNN nodes
  4,519 against 4,545), no QPB tail. Together with the wil50 sample (0 false
  give-ups of 32, -26%), H = 200 is the new recommendation with PRSM.
- **The gap-ratio rule at 0.65 is not safe on nug30.** It gave up on about
  40 nodes that fathom under H = 300: 1,800 more DNN nodes and a million-node
  QPB tail. CPU still came out lower only because nug30's QPB tail is cheap;
  at wil50 each false give-up costs 47 DNN children. The ratio was fit on four
  wil50 traces without warm starts; the warm-started nug30 fathoms evidently
  flatten below 0.65 somewhere after iteration 150. Not recommended as is. A
  re-fit (later start, or a higher threshold) should be done on nug30 traces,
  where the truth is known, before it is tried on wil50 again.

## Branching lines at wil50 clear misses (`lines/`)

`lines.py` on the 3 depth-3 clear misses of `lines/parents.txt`, 10 lines
each, 1,410 child bounds at 300 iterations; 8.5 hours of wall time (about
3x the cloud's estimate: a child costs 70 s here, and the console QPB
runs are serial). Log `lines/wil50_lines_m2pro.log`, cache `lines/work/`.
Cost = sum over children of 47^(gap/130), an estimate of the subtree size.

| parent | short | best line (selector) | cost, best | ecc line | lookahead's line | xbar's line |
|---|---|---|---|---|---|---|
| 1 `34:40,28:22,14:29` | 405 | col 45 (xbar) | 8.5e4 | 4.68x best (8th of 10) | 1.58x | 1.00x |
| 3 `27:44,5:16,17:9` | 241 | col 45 (xbar) | 5,143 | 1.96x (5th) | 1.11x | 1.00x |
| 10 `39:27,30:48,25:49` | 158 | row 49 (U, xbar, lookahead) | 486 | 2.07x (4th) | 1.00x | 1.00x |

- **Eccentric location (rule 6) is a poor line at wil50 clear misses**: 2 to
  4.7x the best line's cost, and on parent 1 worse than two random lines.
  wil50's locations lie on a grid with a smooth flow matrix, so the most
  eccentric location has no special leverage.
- **The DNN primal marginals pick the best line at all three parents** (the
  row or column with the largest sum of x_ij^2), at no cost once the node's
  DNN bound has run. Lookahead (rule 4, one QPB per child) is as good on
  parent 10 and 1.1 to 1.6x worse on the others.
- **Line scores** (Spearman with line cost, per parent): sum of x^2 +0.93 /
  +0.73 / +0.84; children with x > 0.02 +0.93 / +0.75 / +0.92; mean child QPB
  (lookahead's criterion) +0.88 / +0.70 / +0.72. The nug30 favorite,
  children with x > 0.05, is -0.09 on parent 1; use sum of x^2 or the 0.02
  count.
- At the child level, pooled over 1,410 children, the child's QPB predicts
  its DNN gap far better than x_ij (Spearman -0.94 against +0.50), unlike on
  nug30 where they tied. x_ij works at the line level, which is what the
  branching rule needs.
- Column 45 was the xbar pick at all three parents (and the best line at
  two); location 45 may be structurally special in wil50. Worth checking that
  the rule is not just rediscovering one location.

**What it is worth.** About 95% of the estimated wil50 tree is below the
depth-3 clear misses (../NOTES.md section 4). Replacing the eccentric line
with the xbar line at those nodes divides their subtrees by 2 to 4.7x on
these three samples, and the same choice is available again at every DNN
level below. This is the largest lever found so far that costs nothing per
node; it needs a branching rule in branch.C that reads the DNN primal
marginals (dnn_bound would return the x_ij of the last iterate alongside the
warm start) and picks the line with the largest sum of x_ij^2, then an A/B
on nug28/nug30 (where ecc was tuned) and a re-sample of wil50 depth-4
nodes under the new rule.
