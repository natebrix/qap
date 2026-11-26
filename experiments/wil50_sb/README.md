# wil50: what would DNN strong branching choose at the top of the tree?

Nov 25-26, 2026. Question: would strong branching (bounds one level down) at depths 0 and 1 be
likely to improve a fresh wil50 run under the recipe `infra/phase1/recipes/wil50_d2r10_g150xl2.par`?
Stage 1, here: fire it once at the root and at eight of the fifteen depth-1 nodes and see which line it
takes, by what margin, at what cost. No subtree was solved under a changed line: the gain below is a
prediction from a fitted model, not a measurement.

`par/sb.par` is the recipe plus `dnn.depths = 0:48`, `sb.depths = 0:1`, `sb.all = 0` (every row at the
root; 4 rows and 4 columns by marginal score at depth 1; 150 iterations a child; score the sum of
shortfall^P, P = 8 at the root, 1 at depth 1). `run.sh [NODE]` runs one firing (resumable), `queue.sh`
kept four going, `snap.sh` pushed the progress files while they ran. `scores.py` prints a line's score
from `runs/*.progress.txt`; `runs/*.sb.txt` is the solver's own trace of each decision.

## Lines chosen

| node | share of the estimate's rest | recipe's line | strong branching | margin over the next line | CPU h |
|---|---|---|---|---|---|
| root | | r21 | r45 | 661.6 against r21's 670.0 (8th-power norm) | 14.7 |
| 21:0 | 38% | r45 | r45 | 6,831 against c49's 8,878 | 6.8 |
| 21:10 | 37% | c0 | **r45** | 6,555 against c0's 8,226 | 6.7 |
| 21:1 | 8% | c0 | **r45** | 5,625 against c0's 6,652 | 6.7 |
| 21:20 | 7% | r45 | r45 | 3,429 against r26's 6,583 (rows only: symmetry left) | 2.1 |
| 21:11 | 4% | r45 | r45 | 4,637 against c0's 6,564 | 6.8 |
| 21:21 | | r45 | r45 | 2,215 against r26's 4,959 (rows only) | 2.1 |
| 21:2 | | not looked up | r45 | 3,860 against c0's 4,882 | 6.8 |
| 21:12 | | not looked up | r45 | 2,630 against c0's 4,064 | 6.7 |

59.4 CPU hours in all. Shares are of `wil50-est-261124`'s 10,145 CPU hours.

- At the root rows 45 and 21 are near mirror images: location by location their children's shortfalls
  agree to about 20 (sums 6,606 and 6,636, the hardest child 545 on both). The choice between them is
  a coin flip.
- Under r21, strong branching takes r45 at every depth-1 node fired. The recipe does too, except at
  21:10 and 21:1 (c0), which hold 45% of the estimated work.

## Is the changed line worth anything? (`predict.py`, `models.py`)

A child's subtree work modelled as a function of its 150-iteration shortfall (A exp(k s), or s^P), scale
and shape fitted to the estimate's shares at the five nodes where the recipe's line was scored
(21:0, 21:10, 21:1, 21:20, 21:11), then applied to r45 against c0:

| node | r45's work relative to c0's, over every model shape |
|---|---|
| 21:10 | 0.73 to 0.83 |
| 21:1 | 0.81 to 0.88 |

That is 7 to 11% of the tree (700 to 1,100 CPU hours), against at most about 107 CPU hours for the firings
at the root and all fifteen depth-1 nodes. Reasons not to lean on it:

- the models reproduce the five shares to a factor 1.3 (best) to 2.1 (the plain sum), and the estimate
  itself resolves 30 to 40%: the predicted gain is smaller than either;
- r45 wins on the easy children, not the hard ones: at 21:10 its median child is 102 short against
  c0's 168, the five hardest 308 to 355 against 314 to 355.

A hint about the score: the shares are fitted far better by a tail-weighted score (P = 12 to 16, or
k = 0.04 to 0.05: a factor 1.3) than by the sum used at depth 1 (P = 1: a factor 2.1); 21:0 and 21:1
have nearly equal sums on the recipe's line and 38% against 8% of the tree. Under the tail-weighted
fits c9 and c40 come out 20 to 40% ahead of r45 at 21:0 and 21:10. Five points and the noisiest
children: a hypothesis (`sb.power` at depth 1 well above 1), not a result.

One child at 21:11 (row 25) ended 769 short, far below its parent's bound: a failed solve, which
makes r25's figures there meaningless under the tail-weighted models. r25 is not a contender.

## What would settle it

A second estimate run under the recipe with strong branching at depths 0 and 1 (about 750 CPU hours;
it resolves 30 to 40% differences, so a 10% gain would not show in the total, only in the 21:10 and
21:1 subtrees compared like for like), after deciding the depth-1 exponent, since one run tests one
setting.
