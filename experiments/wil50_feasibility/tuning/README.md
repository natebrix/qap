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

## Rule 8 in the solver: nug30 A/B (`lines/nug30_xbar_H200_j8_m2pro.log`)

Branching rule 8 (commit 622af94, `BRANCH_DNN_XBAR`): at a node where the DNN
bound ran and missed, branch on the free row or column with the largest sum
of squared primal marginals x_ij; lookahead (rule 4) where no DNN bound ran.
`PARAM/nug_deep_xbar.par` is `nug_deep.par` with rule 8 in the three upper
blocks. Same binary and settings as the H = 200 run above (PRSM, degree 4,
`-x 2:7:1:200:1 --ub 6124`, 8 workers).

| | `nug_deep.par` (rules 4/7/6) | `nug_deep_xbar.par` (rule 8) |
|---|---|---|
| nodes | 4,529 | 3,946 |
| DNN nodes / fathomed | 4,519 / 4,356 | 3,935 / 3,796 |
| DNN CPU | 4,876 s | 4,125 s |
| CPU | 4,904 s | 4,133 s (-16%) |
| wall | 11 min 17 s | 9 min 56 s |

Per depth (nodes / fathomed, from the checkpoints; depth 0-2 are identical
because depth-1 nodes branch by lookahead in both):

| depth | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| rules 4/7/6 | 1,552 / 1,467 | 2,247 / 2,228 | 470 / 469 | 25 / 25 | |
| rule 8 | 1,552 / 1,493 | 1,593 / 1,576 | 442 / 438 | 100 / 99 | 24 / 24 |

The depth-2 misses are the same 58 nodes; rule 8's choice of line for them
leaves 59 depth-3 misses instead of 85, so depth 4 has 1,593 nodes instead of
2,247. The depth-3 children also bound faster (1.10 s against 1.31 s per
node): children of a concentrated line are nearer to fathoming when the warm
start begins. The tail is a little deeper (one depth-6 miss), which is noise
at this level. Rule 8 is a 13% smaller DNN tree and 16% less CPU on nug30,
where rule 6 was tuned; on wil50 the line comparison above predicts a much
larger gain because eccentric is so poor there. The picks were rows and
columns both (`QAP_XBAR_DEBUG=1`).

### nug28: rule 8 loses (`lines/nug28_{ctl,xbar}_H200_j8_m2pro.log`)

Same settings, `--ub 5166`:

| | `nug_deep.par` | `nug_deep_xbar.par` |
|---|---|---|
| DNN nodes / fathomed | 1,674 / 1,613 | 2,928 / 2,813 |
| CPU | 1,130 s | 1,546 s (+37%) |
| depth 3 nodes / fathomed | 590 / 561 | 590 / 537 |
| depth 4 | 725 / 719 | 1,265 / 1,234 |

The same 25 depth-2 misses; rule 8's lines leave 53 depth-3 misses where
lookahead leaves 29. So on nug instances the single-line marginal rule is
high-variance: nug30 -16% CPU, nug28 +37%. That matches the nug30 line
check, where the marginal row was the best of 11 lines and the marginal
column the worst. On wil50 both marginal lines were the top two at all three
parents, so the rule's case there stands, but it should not replace lookahead
on nug-type instances as is.

Candidate fix, not yet built: use the marginals to shortlist (the top row and
top column by sum x^2, plus lookahead's own candidates) and let lookahead's
child QPB bounds choose among them. Lookahead costs 2(n-d) QPB solves per node,
negligible next to a DNN bound at these depths, and on wil50 the mean child
QPB ranked lines at +0.70 to +0.88, so it would separate the two marginal
lines. The decisive test for wil50 is a depth-4 re-sample under rule 8 (the
tree estimate in ../NOTES.md section 4), which does not depend on nug28.


# Two settings ported from Jean-Pierre Goux's nug30 solver (Nov 1 2026)

J.-P. Goux's solver (github.com/jpgoux/nug30, journal items D6 and D8)
reaches nug30 in 44 min 39 s on a 4-core 2.1 GHz Xeon with a set of changes
tuned on nug25 to leave the tree unchanged. Two of them are policies rather
than tunings and are now in `bnb.C`/`dnn.C` as opt-ins:

- `QAP_DNN_EXT=k` (D8): a node below the last DNN depth whose parent's DNN
  bound ran and did not fathom keeps the DNN bound, k levels down, warm-started
  as usual; under a fathomed or QPB-only node nothing changes. The flag
  travels with the parent's `WarmStart`, so not across a checkpoint restore.
  His motivating case: a depth-7 give-up at 6,116.7 produced 345,000 QPB nodes
  (99.8% of his tree); ours: 580,000 QPB nodes under two depth-7 misses in the
  laptop nug30 run, 3.5 million under a depth-5 miss on nug28 (../../eigensolver
  README section 8).
- `QAP_DNN_GU_ZONE=z` (D6; `dnntest -Z z[:d]`): a give-up decided while the
  certified bound is within z|cutoff| of the cutoff is deferred d = 100
  iterations and must then fire again. His lost fathoms under faster settings
  were all give-ups within 0.3% of the cutoff; so were ours under H = 150
  (section 1 above).

Functional checks (cloud Xeon, machine loaded by the regional-cut runs):

| check | without | with |
|---|---|---|
| nug20 `default.par -u 2570 -x 2:2:1:300:1`, PRSM + degree 4 | 84,592 nodes (DNN 97 / 95 fathomed), 36.1 s | `QAP_DNN_EXT=3`: **141 nodes** (DNN 133 / 131), 26.6 s: the 36 depth-3 children of the two misses all fathom |
| wil50 node 12 `38:40,40:27,18:29`, `dnntest -H 150` (a false give-up) | gives up at 170 at 48,703.4 | `-Z 0.003`: deferred at 170 to 270, **fathomed at 280** (48,814.9) |

`make test` and `make dnncheck` pass with the port. A/B in the solver:
`goux/run_ab.sh` (nug25 `-x 2:5:1:300:1` and nug28 `-x 2:5` with and without
`QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003`, 4 workers); nug30 on the laptop.

## Baseline with Goux's D6/D8 (`giveup/nug30_ab/nug30_prsm_H200_ext2_zone003_j8_m2pro.log`)

Binary at 01025c9 (QAP_DNN_EXT, QAP_DNN_GU_ZONE merged), default rules,
PRSM, degree 4, `-x 2:7:1:200:1`, `QAP_DNN_EXT=2 QAP_DNN_GU_ZONE=0.003`,
8 workers. Tree identical to the H = 200 run without the flags (4,519 DNN
nodes, 4,356 fathomed; the same per-depth counts): with H = 200 there is no
depth-7 miss for the extension to act on, and the give-up confirmation
changed no decision. CPU 5,088 s against 4,904 s; part of that is the
deferred give-ups' extra iterations, part a few minutes of nug20 tests that
shared the machine. This is the reference for the rule 9 runs below.

## Rule 9 A/B (`lines/nug2{8,30}_rule9_H200_d6d8_j8_m2pro.log`, control `lines/nug28_ctl_H200_d6d8_j8_m2pro.log`)

Rule 9 (commit ae5595f): lookahead restricted to rule 8's two lines, the
marginal row and the marginal column, chosen by the children's QPB bounds.
Same flags as the D6/D8 baseline.

| | nug28 rules 4/7/6 | nug28 rule 8 | nug28 rule 9 | nug30 rules 4/7/6 | nug30 rule 8 | nug30 rule 9 |
|---|---|---|---|---|---|---|
| DNN nodes / fathomed | 1,674 / 1,613 | 2,928 / 2,813 | 1,938 / 1,866 | 4,519 / 4,356 | 3,935 / 3,796 | 4,346 / 4,190 |
| depth-3 misses | 29 | 53 | 32 | 85 | 59 | 73 |
| CPU | 1,105 s | 1,546 s | 1,190 s (+8%) | 5,088 s | 4,133 s | 4,450 s (-13%) |

(rule 8 rows from the runs without D6/D8; the flags did not change the tree.)

Rule 9 hedges as intended: it repairs most of rule 8's nug28 loss (+8%
instead of +37%) and keeps most of its nug30 gain (-13% instead of -16%,
against the 5,088 s baseline), but it beats the default rules on neither
instance decisively. On nug instances full lookahead (use_best 30 at depth
<= 2 scores nearly every line) is already a good selector, and a shortlist of
two lines can only match it by luck. On wil50 the picture was the reverse:
lookahead's own line was 1.0 to 1.6x the best and the marginal lines were
the top two at all three parents. So: keep rules 4/7/6 for nug-type
instances; for wil50 the evidence favours rule 8 or 9, and the decisive
test is wil50 itself, a depth-4 (and deeper) re-sample under rule 8, rule
9 and eccentric from the same depth-3 misses, bounding every child with the
same budget as lines.py did. That is hours of laptop time or a cloud job.

Untried: adding the two marginal lines to lookahead's candidate set at the
depth >= 3 blocks (where use_best is 3 or 5), instead of replacing it.

## Goux's momentum (D2/D4) against our PRSM step: negative (Nov 1 2026)

Fast ADMM (Goldstein, O'Donoghue, Setzer, Baraniuk 2014) with Goux's three
safeguards is in `dnn.C` as `DNNParams::momentum` (`dnntest -M eta[:cap]`,
solver `QAP_DNN_MOMENTUM`): Nesterov momentum on (Y, Z), restart at the
current point when the combined residual grows by more than eta, back to the
last true iterate with the cap halved when a certificate falls below the
running best, plain iterations within 1e-3 |cutoff| of the cutoff. Certificates
stay valid (they are valid for any multiplier). Iterations to fathom, degree-4
filter, one core:

| node | ADMM (1.618, 0) | ADMM + momentum 1.05 | PRSM (1.0, 1.0) | PRSM + momentum 1.05 | PRSM + momentum, cap 0.25 |
|---|---|---|---|---|---|
| nug20 `1:7,19:11` (m = 18) | 260 | 200 | **160** | 230 | 200 |
| nug30 `4:10,6:8` (m = 28) | 490 | 370 | **310** | false give-up at 210 | 350 |

On plain ADMM the momentum reproduces Goux's gain (−23%). On top of the
Peaceman-Rachford step it loses at every cap tried (0.1 to 1.0), and the
certificate oscillation it causes made the extrapolation give-up rule abandon
a fathomable node. The two are substitutes, and PRSM is the better one
(310 against 370 iterations). Kept as an opt-in, off; not run on the 32
wil50 nodes. Goux's other iteration settings, a full eigendecomposition every
120 iterations (`-r 120`, `QAP_DNN_REFRESH`) and beta x 0.7 (`-B 0.21`,
`QAP_DNN_BETA`), use existing parameters and are on the all32 A/B queue
(`../solvers/accel/run_all32.sh`, configs `d4r120`, `d4B021`).
## nug30 refinements (`refine/nug30_xall_*_j8_m2pro.log`)

Four runs on the rule 8 setting, each changing one thing on top of the last
(binary at ae5595f, PRSM, `--ub 6124`, 8 workers, no D6/D8 flags). The
parameter file `PARAM/nug_deep_xbar_all.par` lowers the third block's gap
threshold from 0.05 to 0, so every node through depth 6 branches with rule 8;
with `nug_deep_xbar.par` about 40% of depth-3 nodes and 85% of depth-4
nodes fell through to the bottom block (rule 2). Baselines: default rules
H = 200 4,904 s; rule 8 4,133 s.

| run | change | DNN nodes / fathomed | depth-2 / 3 / 4 misses | CPU | wall |
|---|---|---|---|---|---|
| `xall_c4` | rule 8 for all nodes thru depth 6 | 3,936 / 3,797 | 58 / 59 / 17 | 4,008 s | 567 s |
| `xall_c3t2` | + `QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2` | 3,964 / 3,824 | 59 / 59 / 17 | 3,653 s | 526 s |
| `xall_c3t2_d1` | + DNN from depth 1 (`-x 1:7:1:200:1`) | 3,842 / 3,700 | 68 / 45 / 14 | **3,632 s** | **523 s** |
| `xall_c3t2_h150` | c3t2 with H = 150 | 4,044 / 3,901 | 61 / 59 / 17 (+1 at depth 6) | 3,800 s | 582 s |

- Opening rule 8 to the small-gap nodes changed no decision: the tree is the
  same as `nug_deep_xbar.par` (3,946 nodes), so none of the rule-2 nodes was
  a DNN miss. The 3% CPU difference is run-to-run noise.
- Degree 3 with the tight lam_hi (the accel recommendation, so far tested
  only on wil50 nodes) saves 9% of CPU in a full search at effectively equal
  strength: one extra depth-2 miss, otherwise the same per-depth counts.
- DNN from depth 1 costs 36 s on the nine depth-1 nodes and gives rule 8
  marginals for the depth-1 branching. It produces 10 more depth-2 misses
  (68 against 58) but their children fathom far better: 45 depth-3 misses
  instead of 59, and 1,215 depth-4 nodes instead of 1,593. Net 3,632 s, the
  lowest nug30 so far: -26% against today's default-rule baseline (4,904 s)
  and -57% against the Oct 29 laptop run (8,355 s, H = 300, no PRSM).
- H = 150 without the gap-ratio rule is unsafe on its own: two extra depth-2
  misses and a depth-6 miss (24 extra depth-7 nodes), +4% CPU. H = 200 stays.

Recommended nug30 setting: `nug_deep_xbar_all.par`, `-x 1:7:1:200:1`,
`QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2`.
The rule 8 caveat from the nug28 A/B stands: this setting has not been
re-run on nug28, where rule 8 lost 37%.

## Goux's D6/D8 in the solver: nug25 and nug28 A/B (`goux/run_ab.sh`, cloud Xeon 2.1 GHz, 4 workers, PRSM + degree 4)

| run | setting | nodes | DNN nodes / fathomed | total time (s) |
|---|---|---|---|---|
| nug25 base | `nug_deep.par -x 2:5:1:300:1` | 696 | 689 / 658 | 485 |
| nug25 + `QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003` | | 696 | 689 / 658 | 490 |
| nug28 base | `nug_deep.par -x 2:5` | **3,500,289** | 1,602 / 1,543 | 7,361 |
| nug28 + `QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003` | | **1,634** | 1,625 / 1,566 | 6,727 |

nug25: identical trees, +1% time (no last-level miss, no decision changed by
the confirmation). nug28 is the case D8 was ported for: under one depth-5
PRSM miss the baseline sends 3.5 million nodes to QPB alone (eigensolver
README section 8); with the extension those become 23 more DNN nodes, all
fathomed, and the run is 9% faster. Together with Fable's laptop nug30 run
(identical tree at H = 200, +4% from the deferred give-ups), the verdict:
`QAP_DNN_EXT` is a pure safety net that costs nothing where no last-level
miss exists and removes the QPB tails where one does; `QAP_DNN_GU_ZONE` is
insurance worth about 4%. Both belong in the production setting.
### + Goux's beta and refresh (`refine/nug30_xall_c3t2_d1_b021r120_j8_m2pro.log`)

The `xall_c3t2_d1` setting with `QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120`
(binary at 2b9d0fd), both changed in one run:

| | `xall_c3t2_d1` | + beta 0.21, refresh 120 |
|---|---|---|
| DNN nodes / fathomed | 3,842 / 3,700 | 3,820 / 3,679 |
| misses depth 2 / 3 / 4 / 5 / 6 | 68 / 45 / 14 / 5 / 1 | 68 / 46 / 14 / 4 / 0 |
| CPU | 3,632 s | **2,895 s** |
| wall (8 workers) | 523 s | 419 s |

The tree is the same to within a miss per depth, so the 20% is per-node
cost: fewer full eigendecompositions (every 120 iterations instead of 20)
and, with the smaller penalty, no loss of strength. The two were not
separated; Goux's own setting pairs them. nug30 is now 2,895 CPU s against
4,904 s this morning (default rules, H = 200) and 8,355 s on Oct 29.

Recommended nug30 setting: `nug_deep_xbar_all.par`, `-x 1:7:1:200:1`,
`QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2
QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120`. Still to do: separate beta from
refresh, re-run nug28 (rule 8 caveat), and try the same on the wil50 sample.

Addendum (Nov 1, cloud): Goux's `-r 120` and `-B 0.21` on the 32 wil50 nodes
(`../solvers/accel/README.md`): both weaker per iteration at m = 47 (−81 and
−79 at iteration 200 on the misses), `-r 120` nets −11% CPU on fathoms,
`-B 0.21` nets nothing; same fathoms. Good on nug25/nug30 (−20% on the
laptop), not on wil50.

## Cut phase on nug30, refined setting (`refine/nug30_refined_{base,cut,cut_d2}_j8_m2pro.log`)

The cloud's cut phase (`QAP_DNN_CUT_ZONE`, `stronger_bounds/` step 5) on the
refined nug30 setting above (rule 8 everywhere, `-x 1:7:1:200:1`, degree 3
tight, beta 0.21, refresh 120; binary at 61c9b6c), the cloud's nug30 cut
parameters `QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=5:6 QAP_DNN_CUT_SIDE=1
QAP_DNN_CUT_DELTA=20`, 8 workers. `QAP_DNN_CUT_MAXDEPTH` (new, bnb.C) limits
the phase to nodes of at most that depth. parallel.py's merged state does not
carry the cut counters, so cut fathoms = fathomed − DNN fathomed and cut CPU
= total − DNN time.

| | base | cuts at every depth | cuts at depth <= 2 |
|---|---|---|---|
| nodes | 3,821 | 2,086 | 2,782 |
| depth-2 misses (of 225) | 68 | 33 | 31 |
| depth-3 nodes / misses | 1,880 / 46 | 924 / 21 | 868 / 45 |
| depth-4 nodes | 1,242 | 567 | 1,215 |
| cut fathoms / cut CPU | – | 62 / ~1,840 s | 37 / ~640 s |
| DNN time | 2,937 s | 2,101 s | 2,273 s |
| **total CPU** | **2,943 s** | 3,940 s (+34%) | 2,977 s (+1%) |

* At depth 2 the cuts fathom about half the misses (35–37 of 68) at ~17 s a
  phase, but they fathom the near misses, whose children were the cheap
  ones: the depth-3 misses are unchanged (45 against 46) and the depth-4
  count barely moves. Children removed ≈ 600 s against ≈ 640 s of phases:
  break-even.
* At depth 3 the phases cost ~1,100 s and remove little: a loss.
* The contrast with the wil50 node-6 subtree (2.3x, `stronger_bounds/search/`)
  is per-node cost. A nug30 depth-3 child costs ~0.7 s, so a fathomed depth-2
  miss saves ~20 s and a phase must be cheaper than that; a wil50 depth-4
  child costs ~75 s. The cut phase is a wil50 tool; on nug30 under the
  refined setting it does not pay at any depth. nug30 stays at 2,895–2,943 s.

### Nov 3: the recommended setting re-run, plus warm Y, one-pass CholQR and D6/D8 (`refine/nug30_best_nov3*_j8_m2pro.log`)

Binary at a2836ef, 8 workers, `--ub 6124`, `nug_deep_xbar_all.par`,
`-x 1:7:1:200:1`, PRSM, degree 3 tight, beta 0.21, refresh 120 (the Nov 1
recommended setting), then the same with the cloud's Nov 3 per-node levers
that apply on Apple silicon (the AMX kernel is Intel-only and is skipped):

| run | added | nodes / DNN fathomed | tasks | CPU | wall |
|---|---|---|---|---|---|
| `nug30_best_nov3` | — | 3,821 / 3,679 | 300 | 2,870 s | 459 s (x6.25) |
| `nug30_best_nov3_warmy_cq1_ext` | `QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003` | 3,844 / 3,701 | 345 | **2,156 s** | **346 s** (x6.22) |

The baseline reproduces the Nov 1 tree exactly (3,821 nodes) at 2,870 s
against 2,895 s. The three levers together take 25% off the CPU and the
wall at a tree 23 nodes larger (the warm Y changes the iterates, so a few
give-up decisions move; EXT confirms the depth-7 children). The per-child
1.17x the cloud measured on wil50 children carries to nug30 as 1.33x
overall, so one-pass CholQR is adding on top; the three were not separated.
nug30 is now 2,156 CPU s against 8,355 s on Oct 29 (3.9x).

Recommended nug30 setting: `nug_deep_xbar_all.par`, `-x 1:7:1:200:1`,
`QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2
QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1
QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003`. Not yet re-run on nug28.
