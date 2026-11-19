# The landscape signal: does a child's best completion predict whether it branches, and which line to take? (Nov 19 2026)

branching_lab agent "landscape". Brief: test the bound lab's survive rule ("a child branches iff the best
completion a pairwise-exchange descent finds for it costs less than the incumbent plus the bound's slack")
on the production trees, and whether counting the children it predicts to survive ranks lines.

## Answer

**At the nodes that hold the best known solution (the chains), the rule holds in production: it predicted
the taken line's surviving children at tai60b's 21 chain nodes of depths 0-10 with no false positive and 8
misses in 110 [m]. Among every free line it prefers production's line at depths 0-2 and 5-6, and it strongly
prefers other lines at depths 7-8, where rule 2's c54 and c2 leave 9 and 7 survivors and lines exist that it
says leave 1 or 2 (the path child alone, or with one more child).** The depth-7 branching is 35% of tai60b's
nodes, so if the prediction holds there the saving is large. That c54 is a bad line there is robust to the
search's starting point. Which line is best is not: c47 holds at 2, while r31 moves from 1 to 5. That is a prediction about lines nobody has
bounded; one cheap pair of arm runs below tests it. Off the chains (inside the subtrees of nodes whose own best
completion is above the cutoff), the rule's threshold is wrong: recall 0.09 to 0.40 [m]. There the children
branch at excesses up to 3 to 10 times the parent's shortfall, and more search effort does not fix it. Its
ranking of children within a parent still beats the children's QPB bounds (0.80 against 0.55 on tai60b). On
the grid truth sets (tho30, tho40) the signal ranks lines worse than the children's first DNN solve (mean
top line 1.3x to 1.5x the best, against 1.00x), as the brief expected. It costs 1 to 2.3 CPU s per node for
all lines with descents alone, and 31 to 69 s with the 30-kick search used here (n = 50, 60) [m], so it can
serve as a free shortlist for the DNN lookahead. It needs a near-optimal incumbent permutation, not just its
value. Without the best known solution, an ILS from random starts found one 13,483 above it on tai60b in
25 s (see "Without the optimum").

## Conventions (checked, `check.py`)

A permutation p is 0-based with p[i] the location of facility i (QAPLIB's `.sln` minus 1), and the cost is
sum A[i,k] B[p(i),p(k)] with A the instance file's first matrix. A trace fix `r:v` means p[r] = v. Checked: on
every production trace, the nodes consistent with p0 under this reading are all B (or E at the leaf), down to
the leaf (tai60b 58 nodes, tai50b 48, tai35b 33, tai40b 38). Under the inverse reading only 2 are, and they
include D and C nodes. `rN` branches on facility N (children N:v), `cN` on location N (children r:N).
tai60b's mirror is p0 with facilities 43 and 44 exchanged (43:45, 44:58), 2,786 above the best known.

Definitions used throughout. The **excess** of a node is the cost of the best completion found minus the best
known. The **slack** of a parent is cutoff - lb, its shortfall after its own cut phase: the trace's `cut` field,
else `dnn`. The rule says a child survives iff excess <= slack. A child **survives** when its trace outcome is
B (or E). **x** = excess / slack.

## 1. The survive rule in the production trees

Data: every child of every branched node in tai60b-261117 (331,153 children), tai50b-261116 (8,224),
tai50b-r9 (6,000), and line_signals' tai35b_base (1,016) and tai40b_base (730) [m `excess.py`, `rule.py`,
`out/rule_*.md`]. The search: one 2-exchange descent from p0 adjusted to the child's fixes (one exchange per
fix), one from the parent's best completion, and for tai60b one from the mirror (`ex_k`). Each descent is
C with Taillard's delta table, O(n^2) per move (`ls.c`).

**Where the parent is decides it** [m `split.py`, `out/split.md`]:

| run | parents | children | survive | rule, slack = parent's shortfall: precision / recall | "gap" rule (slack + parent's excess): P / R |
|---|---|---|---|---|---|
| tai35b chain | 32 | 575 | 37 | 1.00 / 0.95 | 1.00 / 0.95 |
| tai35b off-chain | 14 | 441 | 9 | 0.69 / 1.00 | 0.50 / 1.00 |
| tai40b chain (the whole tree) | 37 | 730 | 37 | 0.95 / 1.00 | 0.95 / 1.00 |
| tai50b-261116 chain | 47 | 1,183 | 78 | 0.92 / 0.85 | 0.92 / 0.85 |
| tai50b-261116 off-chain | 167 | 7,041 | 136 | 0.89 / 0.31 | 0.33 / 0.80 |
| tai50b-r9 chain | 27 | 1,073 | 55 | 0.92 / 0.82 | 0.92 / 0.82 |
| tai50b-r9 off-chain | 116 | 4,927 | 88 | 0.88 / 0.40 | 0.31 / 0.92 |
| tai60b chain | 194 | 4,602 | 347 | 0.94 / 0.83 | 0.92 / 0.89 |
| tai60b off-chain | 7,302 | 326,551 | 7,149 | 0.96 / 0.09 | 0.43 / 0.71 |

"Chain" means the parent's own best completion is within 3,000 of the best known: it holds p0, the mirror, or
another solution that close. The bound lab's 270-of-286 was measured on dive nodes next to the path, which are
of this kind. Calibration at tai60b's chain parents of depth <= 11 [m]:

| x = excess / parent's shortfall | <= 0.5 | 0.5-0.8 | 0.8-1.0 | 1.0-1.2 | 1.2-1.5 | 1.5-2 | 2-3 | > 3 |
|---|---|---|---|---|---|---|---|---|
| children | 47 | 36 | 13 | 26 | 20 | 29 | 37 | 942 |
| branched | 47 | 36 | 13 | 10 | 2 | 1 | 0 | 1 |

That is a sharp step at x = 1 with a soft edge to about 1.3. At tai60b's 43:58, the line r1 production took
has six survivors at x = 0 (the path child), 0.22, 0.46, 0.50, 0.61 and 0.72. The cut phase fathoms the next
six, at x = 1.01 to 1.22, and the first solve fathoms the next one, at 1.41.

**Off the chains** the parent's shortfall understates the bound's gap at the child. The parent's best
completion is above the cutoff, and its children branch at x up to 3 to 10. Adding the parent's excess to the
slack (the "gap" rule) recovers recall 0.71 to 0.92 at precision 0.31 to 0.43. Neither is a usable threshold.

**Slack variants, all children** (`out/rule_*.md`). The on-path shortfall at the child's depth (`pathd`)
predicts a little better than the parent's shortfall on tai60b (F1 0.40 against 0.21) and tai50b (0.69 against
0.65). At the chain nodes the two are the same number. The best multiplier on `pathd` is 1.25 to 1.5 (F1 0.69
on tai50b, 0.62 on tai60b).

**The kill rule: is it better than a trivial rule?** The trivial rules get, per parent, the true number k of
survivors (an oracle the rule does not have), then take the k children with the smallest excess, the smallest
QPB, or, for reference only, the smallest first DNN bound. Fraction of survivors found [m]:

| run | k smallest excess | k smallest QPB | k smallest first DNN (not available in advance) | the rule, F1 (slack = parent's shortfall) |
|---|---|---|---|---|
| tai35b | 0.98 | 0.78 | 1.00 | 0.94 |
| tai40b | 1.00 | 0.62 | 1.00 | 0.97 |
| tai50b-261116 | 0.84 | 0.41 | 0.90 | 0.65 |
| tai50b-r9 | 0.86 | 0.30 | 0.83 | 0.69 |
| tai60b | 0.80 | 0.55 | 0.86 | 0.21 |

The excess ranks children nearly as well as their own first DNN solve and far better than QPB, so the idea
survives the kill rule. The threshold does not survive off the chains. At the chain nodes, which are where
line choice is worth money (step 2), both rank and threshold hold.

**Sensitivity to the search's effort** [m `effort.py`, `effort_score.py`, `out/effort*_*.md`]. The sample was
every B and C child plus 200 to 400 D children, reweighted. The search was an ILS of K kicks, each 6 random
exchanges then a descent, from the best descent. Precision / recall with slack `pathd`:

| search | tai50b | tai60b (800 B/C + 400 D sampled) | tai35b |
|---|---|---|---|
| descents (ex_k) | 0.92 / 0.55 | 0.94 / 0.22 | 0.94 / 0.96 |
| ILS 30 kicks | 0.88 / 0.61 | 0.93 / 0.25 | 0.92 / 0.98 |
| ILS 300 | 0.83 / 0.77 | 0.93 / 0.28 | 0.90 / 0.98 |
| ILS 1,000 | 0.79 / 0.84 | | |
| ILS 3,000 | 0.78 / 0.84 | | |
| no knowledge of p0: ILS 300 from a random start | 0.84 / 0.48 | 0.91 / 0.10 | 0.92 / 0.76 |
| no knowledge: ILS 3,000 from a random start | 0.79 / 0.70 | | |

On tai50b a descent overstates many children's best completion, and recall climbs to 0.84 by 1,000 kicks,
then stops. On tai60b, whose errors are almost all off the chains, effort barely moves recall: those misses are
the bound, not the search. At the chain nodes 30 kicks change the predicted count of 30% of the lines
(363 of 1,210 on tai60b, 373 of 990 on tai50b), so descents alone are not enough for the line choice.

## 2. Does it rank lines?

### Truth sets [m `truth.py`, `linescore.py`, `out/linescore.md`]

The excess of every child of every candidate line uses ILS 300 from the best known solution and all its images
under the grid's 4 automorphisms (tho), or from p0 and the mirror (tai60b). For tho30, QAPSOLN has no
permutation, so one was found by `ts.c`, a robust taboo search, in 2 seeds (`truth/tho30.perm`). The node's
slack s: tho, a cold `dnntest` solve with give-up 300 at each truth node, from the existing
`/home/user/qap/latest` binary (`truth/*_nodes_dnntest.txt`: 3,537 to 3,875 on tho30, 2,565 to 3,088 on tho40).
tai60b 43:58 uses the line study's arm state, 1.76M. The signals, smaller meaning a better line: predicted
survivors (ex <= s); the same over the node's own best completion (ex - ex_node <= s, where ex_node is the min
over its children); and the chain size model below. In the form of `line_signals/signals_out.md` (7 nodes):

| signal | mean Spearman | best line 1st | in top 2 | in top 3 | its top line / best, mean (worst) |
|---|---|---|---|---|---|
| landscape: predicted survivors, ILS 300 (ex <= s) | +0.58 (defined at 4 of 7) | 1/7 | 2/7 | 2/7 | 1.43x (2.25x) |
| landscape: predicted survivors, descents only | +0.70 (3 of 7) | 2/7 | 3/7 | 3/7 | 1.67x (2.72x) |
| landscape: survivors over the node's own best (ex - ex_node <= s) | +0.53 | 2/7 | 3/7 | 3/7 | 1.30x (1.70x) |
| landscape: chain size model, sum P(x) x^c | +0.53 | 1/7 | 2/7 | 4/7 | 1.47x (2.28x) |
| (line_signals) marginal score | -0.07 | 1/7 | 1/7 | 2/7 | 1.52x (1.99x) |
| (line_signals) children's first DNN solve (a) | +0.80 | 5/6 | 5/6 | 6/6 | 1.00x (1.00x) |

- **tho40 (depth 3)**: every child's excess is above the node's slack (the node's own best completion is
  5,020 to 6,396 above the optimum, about twice s), so the plain count is 0 for every line. Off the optimum the
  rule needs the node's own best completion subtracted, and even then the ranking is weak: Spearman
  +0.21 to +0.88, top line 1.49x to 1.70x.
- **tho30 (depth 1)**: the count ranks lines +0.44 to +0.86, with a top line within 3% of the best. But
  precision child by child is 0.64 to 0.85, and at the 1.95-tolerance even objective of tho the threshold sits
  among many near-equal completions.
- **tai60b 43:58, the line study's four lines**: child by child, recall 1.00 and precision 0.76 (48 branching
  children, 15 false positives, all at x 0.7 to 1.0). The count puts r1, the truth's best, first: r1 12,
  r8 14, r45 17, c26 20. Against the branched children (6, 12, 10, 20) the order is right except r8 and r45.
  Against the rough subtree model it is Spearman +0.40, because r8's subtree is the largest through its path
  child, which every line has once.
- **tai60b 43:58, all 116 free lines** (`all_lines.py`, `out/all_lines.md`, ILS 100). At production's slack
  (1.11M), **r1, production's line, is first of 116 by the count (6 predicted, 6 realised) and by the model**.
  At the arms' slack (1.76M) it is 21st by count. The signal depends on how far the node's own bound has got.

### The production chains [m `chain.py`, `chainscore.py`, `out/chain_*.md`]

At every chain node of depths 0-10, every child of every free line was searched (m^2 children: descents from p0
and the mirror, then ILS 30). s is the node's own shortfall after its cut phase, which the node knows when it
branches. "mod" is the chain size model (`model.py` with CHAIN=1, `out/model_chain_fit.txt`), fitted on the
children of chain nodes of depth <= 15 of all five traces:
P(branch | x) = logistic(0.18 - 4.36 log x), and log subtree = a_parent - 1.07 log x within a parent. A line's
cost is sum P(x) max(x, 0.05)^-1.07. "off" is the same without the chain children, whose subtree (the rest of
the chain) every line has.

tai60b, the path chain (the mirror chain is the same to within one child at each depth,
`out/chain_tai60b-261117_mirror.md`):

| depth | taken | taken: predicted / realised survivors | taken's rank of all lines, by count / by model | least count: line, count | signal's line (least off-chain cost) and taken / it |
|---|---|---|---|---|---|
| 0 | r44 | 4 / 5 | 2 / 1 of 120 | c58, 2 | r44, 1.0 |
| 1 | r43 | 1 / 2 | 1 / 1 of 118 | r43, 1 | r43, 1.0 |
| 2 | r1 | 6 / 6 | 1 / 1 of 116 | r1, 6 | r1, 1.0 |
| 3 | r6 | 8 / 8 | 18 / 11 of 114 | c47, 4 | r51, 1.4 |
| 4 | r7 | 6 / 7 | 10 / 2 of 112 | c38, 3 | c19, 1.4 |
| 5 | r0 | 2 / 2 | 1 / 1 of 110 | r0, 2 | r0, 1.0 |
| 6 | r10 | 3 / 3 | 2 / 2 of 108 | c47, 2 | r31, 1.1 |
| 7 | c54 (rule 2) | 9 / 10 | 77 / 49 of 106 | r31, 1 | c10, 11.0 |
| 8 | c2 (rule 2) | 7 / 7 | 49 / 46 of 104 | r31, 1 | r31, 10.4 |
| 9 | r8 (rule 2) | 4 / 4 | 25 / 15 of 102 | r31, 1 | r31, 5.3 |
| 10 | r4 (rule 2) | 3 / 4 | 14 / 6 of 100 | c10, 1 | c10, 2.5 |

The taken lines' children, child by child: tp 53, fp 0, fn 5 of 605 children (58 branched). On the mirror
chain: 49 / 0 / 3 of 545.

tai50b-261116 (`out/chain_tai50b-261116.md`): predicted / realised survivors of the taken line 11/10, 2/2, 1/1,
**12/6** (depth 3), 6/5, 3/4, 4/4, 2/2, 1/2, 2/4, 1/1. Child by child, tp 33, fp 12, fn 8 of 495: less clean
than tai60b's, and wrong by 2x at depth 3. The signal again agrees with rule 9's lines at depths 0-2 (ranks 1-2)
and prefers other lines at depths 3-9, mostly c37 or r17, with 1 or 2 predicted survivors against production's
2 to 6.

**Where the signal prefers another line, and how far to trust it.**
- **tai60b depth 7, `10:39`: c54 against r31 or c10.** This is the branching that holds 35% of the tree's nodes
  and 24% of its CPU (the postmortem). Under c54 the signal counts 9 children at x <= 1, at 0, 0.16, 0.26,
  0.41, 0.57, 0.61, 0.66, 0.66 and 0.79, and 10 branched. Under r31 it counts the path child alone, with four
  children just above the step (x 1.03, 1.06, 1.11, 1.18) and three at 1.40 to 1.55. The calibration above
  says about 40% of the first four and 10% of the rest branch, so **about 3 survivors against c54's 10 [e]**.
  The model puts c54's off-chain cost at 11x r31's or c10's. If it held at both chains, most of c54's
  off-chain subtrees (114,669 nodes, 873 CPU-h, at both chains) would go, less what r31's 2 or so survivors
  hold. That is of the order of a third of tai60b's tree [e]. Trust: the threshold is calibrated on these very
  chain nodes, with 0 false positives in 1,150 taken-line children. But a line's count rests on its
  borderline children, and r31's sit at x 1.0-1.2, exactly the soft edge, while c54's survivors are well below
  it. With a start that does not know p0 (below), r31's count goes to 5, while c47 stays at 2 (x 0, 0.22,
  then 1.23 and up) and c54 at 9. **That c54 is among the worst lines at this node is robust (rank 77 of
  106 by count with the p0 start; 9 predicted survivors under both starts); which line is best is not. c47 is the safer prediction.**
- **tai60b depth 3, `1:33`: r6 against c47 or r51.** The model's saving here is small (1.4x off-chain): lines
  with 4 predicted survivors exist, but their survivors sit at lower x (bigger subtrees). Do not trust a
  saving here without the lookahead.
- **tai50b**: at depth 3 the signal mispredicted production's own line by 2x (12 against 6), so its
  preferences there (c37) need a check before they mean anything.

## 3. Cost [m `timing.py`, `out/timing.txt`; this shared 4-core box, single-threaded]

All m^2 children of all free lines at a depth-1 chain node (a child r:v is shared by row r's line and column
v's, so m^2 searches, not 2m^2):

| | tai50b (m = 49, 2,401 children) | tai60b (m = 59, 3,481 children) |
|---|---|---|
| one descent per child from the adjusted incumbent | 1.0 CPU s | 2.3 CPU s |
| + ILS 30 kicks | 31 CPU s | 69 CPU s |
| + ILS 100 kicks | 104 CPU s | 223 CPU s |
| (line_signals) DNN lookahead, 8 lines x 30 warm iterations | about 470 CPU s at m = 33 | |
| (line_signals) DNN lookahead, 4 lines x 150 iterations at 43:58 | | 16,576 CPU s |

It is cheap because of the delta table: a descent is O(n^2) per move after one O(n^3) set-up. A child's
start differs from the incumbent by one exchange per fix, so its descent takes about 10 moves, against 70 from
a random start. The ILS loop is Python around C, so its numbers are 2-3x what a C loop would cost [e]. A
further cut, not done: start each child's table from the node's own table updated by one exchange, O(n^2)
instead of O(n^3).

## Without the optimum ("what the signal needs")

The signal needs (a) a cutoff close to the optimum, as every bound does, and (b) **a near-optimal
permutation** from which each child's best completion is a few exchanges away. Random-start local search
does not give (b) child by child: on tai60b an ILS of 300 kicks from a random start reaches recall 0.10
(table above). But one good permutation at the root is enough, and it is cheap to find. `noknow.py`, an ILS
from 4 random starts of 10,000 kicks each (8 s apiece), found tai60b's 13,483 above the best known (16
positions differ from p0; the other three seeds 2.6M to 9.5M above) and tai50b's best known itself (3 of 4
seeds) [m `logs/noknow.log`]. **The chain analysis with that permutation as the only start**, in place of p0
and the mirror (`START=ils`, `ilscompare.py`, `out/ilscompare.md`) [m]: the taken lines' children are tp 51,
fp 0, fn 7 (against 53 / 0 / 5). The predicted survivors of the taken line are the same at 9 of 11 depths
(depth 0: 3 against 4; depth 10: 2 against 3). The line counts correlate +0.72 to +0.97 with the p0-start
counts, and the least count is within one of p0's at every depth. But single lines move: at depth 7, c54
stays at 9, c47 at 2, while r31 goes from 1 to 5 and c10 from 2 to 3. Counts that rest on children just
above the step depend on the start; a line whose predicted survivors are well separated (c54's, c47's) does
not.
In production the run starts from the best known value + 1 and the best known permutation is in QAPSOLN,
so (b) is free. A run that improves its incumbent should rescore with the new one.

## 4. A rule

**LS-shortlist** (as the shortlist of line_signals' DNN lookahead; or alone at chain nodes where it is clear):

- *Where*: at a node whose DNN bound missed and whose own best completion, found by the signal, is within its
  slack: ex_node = min over the children's excess <= s. These are the chain-like nodes, where the threshold is
  calibrated. Elsewhere keep rule 9 or rule 2. Depth <= 10 on tai60b and tai50b, where the subtrees are large.
- *Inputs*: the incumbent permutation p* (QAPSOLN, or the run's best), its value z, the node's fixes, its bound
  after the cut phase lb_N, and s = cutoff - lb_N. Nothing about the children's bounds.
- *Search*: for each free cell (r, v): start = p* with the node's fixes and r:v forced by one exchange each.
  Also start from the node's cached best completion (its parent's search), adjusted the same way. A 2-exchange
  descent, then 30 kicks of 6 random exchanges with a descent each, keeping improvements. Then
  ex(r, v) = best - z.
- *Score* per line: cnt = #children with ex <= s. Tie-break by
  off = sum over children with ex > 3,000 of P(x) max(x, 0.05)^-1.07, x = ex / s, with
  P(x) = 1 / (1 + exp(-(0.18 - 4.36 ln x))) (`out/model_chain.json`). Then by rule 9's marginal score.
- *Alone*: take the least (cnt, off). *As a shortlist*: the K = 4 lines with the least (cnt, off), plus rule
  9's two lines, go to `QAP_LA_LINES` with a cap of 150 iterations (line_signals Verdict 4). The shortlist costs
  about 70 CPU s at n = 60 against 2-4 core-hours for the lookahead itself.

**The cheapest test, before any whole tree**: one arm pair at tai60b's depth-7 chain node `10:39`, the
prediction with most money on it. Bound the children of c47 (2 predicted survivors under both starts; r31 as a
third arm if a core is free) and of c54 (production's line: 10 branched in the trace, 9 predicted), from one
parent state. That is the line study's design at depth 7 (`tai60b_lines/run.sh`,
experiment build with `QAP_FORCE_LINE`): a parent run that bounds the depth-6 chain node `44:45,43:58,1:33,6:24,7:19,0:32` with `-F`,
`QAP_FORCE_LINE=6:r10` and `-d 7 -o parent.state`. `cutstate.py` keeps `10:39` (as `run.sh arm` does for
43:58). Then the arms `QAP_FORCE_LINE=7:c47` and `7:c54`, each `-r` from that state with `-d 8`, under
`infra/phase1/recipes/tai60b.par` and `-u 608215055`. It costs about 53 children
x 45 to 1,000 CPU s, 2 to 4 CPU-h an arm [e]. Kill rule: if c47 leaves 5 or more branching children, or as many as
c54's arm from the same state, the depth-7 saving is not there and the rule is a shortlist at best. If
it passes, the whole-tree test is tai35b through the coordinator's local mode (about 7,000 CPU s an arm),
with `QAP_LA_LINES` fed from the LS-shortlist. That needs the shortlist in C (the descents above are ~100
lines of `ls.c`), or a precomputed `QAP_FORCE_LINE` per depth, which only works on a path-like tree such as
tai40b's.

## What was not tested

- Any bound of a line production did not take, so every "the signal prefers X" is a prediction (section 2).
  The truth sets are the only realised comparison, and on tai they are one node with four lines.
- The off-chain regime as a line ranker: the threshold fails there, and no truth set has tai off-chain nodes.
- A slack that the node does not know in advance (the children's own first solve) was deliberately not used.
- Whole trees, and a C implementation inside the solver.

## Files

Scripts (all in this directory; rerun in this order): `ls.c` / `lsq.py` (descents, conventions, traces),
`check.py`, `excess.py RUN` (per-node excess, `out/RUN.tsv`; tai60b with `-1 norc`), `rule.py RUN`,
`split.py`, `effort.py` / `effort_score.py`, `model.py` (CHAIN=1 for the chain model), `truth.py SET`,
`linescore.py`, `all_lines.py`, `chain.py RUN 11 30` (ONLY=mirror, START=ils), `chainscore.py`, `timing.py`,
`noknow.py`, `ilscompare.py`, `ts.c` / `tsq.py` (robust taboo search, used only for tho30's optimum),
`autom.py` (grid automorphisms). The over-2 MB outputs are listed in `.gitignore` and rebuilt by these
scripts: `out/tai60b-261117.tsv`, 59 MB, 10 min; the `out/*.pkl` files; and `chain/*.tsv`, 20 min a run.

Note, fixed during the work: the first version of `excess.py`, `truth.py` and `chain.py` built tai60b's mirror
wrongly (it exchanged the facilities at locations 43 and 44); `effort.py`'s tai60b sample also ran with it.
For the excesses this makes no difference that matters: p0 adjusted to any node under 44:58 already has the
mirror's two entries, and `ex_par` carries the mirror chain's completions down. The mirror-chain run
(`ONLY=mirror`) and the ILS-start run used the corrected one (`lsq.mirror`).
