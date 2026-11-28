# line_truth: what a branching line costs, read from wil50's own tree

Nov 27-28, cloud container (4 cores), session branch `claude/cool-babbage-ilnogj`. The question: at the nodes
that carry wil50's work, is the line the run branches on the cheapest one, and can a strong-branching firing
find a cheaper one? The evidence is the published tree of the live run `wil50-261126` (qap `coord/tree.py`,
fetched through the site's proxy by `fetch_tree.py`, never from the bucket) and firings run here at its nodes.

**The result.** At wil50's four heaviest depth-2 nodes the run branches on row 31 (rule 10: every line, the size
score on the children's QPB bounds; it takes r31 at all 162 depth-2 nodes). A DNN strong-branching firing at
each (`sb.depths = 2:2`, four rows and four columns by marginal score, r31 added) ranks r31 last of the nine
lines by its own score and takes a column, c9 or c49. Priced by the size law of the run's own subtrees, the
firing's line costs **0.57 of row 31's** over the four nodes, and the cheapest line (c9 at every node) 0.49.
[m: firings; price: model, checked three ways below]

## The pieces

| file | what it does |
|---|---|
| `fetch_tree.py RUN OUTDIR` | a public `--trace` run's tree, shard by shard, through the site's proxy (resumable: use a new directory for a fresh copy) |
| `tree.py` | the shards as one table; which subtrees are closed |
| `where.py`, `where.txt` | the size law by depth (ln of the CPU under a node against its shortfall), and what a gate on the shortfall would fire on |
| `first_look.py`, `first_look.txt` | depth-1 firings (`../wil50_sb`) against the run: a child's 150-iteration shortfall against the work under it |
| `score_forms.py`, `score_forms.txt` | which function of the children's shortfalls orders nodes by work (plain sum, powers, exp, largest) |
| `truth_powers.py`, `truth_powers.txt` | the same forms at the six tho30/tho40 nodes where every line's subtree was solved |
| `lines_price.py SHARDDIR d1\|d2`, `lines_price_d1.txt`, `lines_price_d2.txt` | every line a firing solved in full, priced by the law fitted on the run's own line |
| `price_check.py`, `price_check.txt` | that pricing at the truth nodes: fitted on one line, does it predict the others? |
| `line_effect.py`, `line_effect.txt` | in the live tree, does the work at equal shortfall depend on the line a child came from? |
| `d2/run.sh NODE`, `d2/par/d2.par`, `d2/runs/` | the depth-2 firings (9 lines x 48 children, 150 iterations, no early exit; 5.9 CPU h each) |
| `d2/price.py`, `d2/price.txt` | the same lines priced from the tree's first-solve shortfalls (an earlier form; kept for the 10-unit offset below) |
| `d2/kids.sh`, `d2/kids.py`, `d2/kids/` | children of c9, c49 and r31 at one node solved to the end: the direct check |

## The size law, and the price of a line

Under a node of this tree the work grows exponentially in the node's shortfall (cutoff less its bound): among
branched nodes at depth 2, ln(CPU h) = a + 0.025 s, a factor 12 per 100 units, R2 0.97 (`where.txt`); the slope
falls with depth (0.021 at depth 3, 0.017 at depth 6). A firing gives each child of each candidate line a
shortfall after 150 iterations. Where the run took one of those lines, its children have both a firing
shortfall and a subtree in the run, so the law can be fitted in the firing's own units: over r31's 178 children
at the four nodes, ln h = -3.27 + 0.0262 s, rms factor 1.82 per child, R2 0.86 (tree of Nov 28, 12:00 AM PT).
A line's price is the sum over its open children of exp(a + k s).

A firing's shortfall is not the run's: for the same child it reads 11 units below the shortfall after the run's
first DNN solve (rms 12; `d2/price.txt`), and 28 to 49 above the node's final one. Price lines from one source.

## Depth 2: the four heaviest nodes

All four are under root r21 and r45; bounds 48,489 to 48,521 against the cutoff 48,816.05 (shortfalls 295 to
327, relative gap 0.11 to 0.12). `lines_price_d2.txt`:

| node | the run's r31: price, CPU h under it so far (subtrees open) | the firing's line (plain sum): price / r31's | the cheapest line: price / r31's |
|---|---|---|---|
| 21:0,45:40 | 466 h, 547 h (37 of 48) | c9: 0.47 | c9: 0.47 |
| 21:10,45:0 | 896 h, 449 h (37 of 43) | c49: 0.60 | c9: 0.50 |
| 21:1,45:0 | 353 h, 362 h (32 of 48) | c49: 0.64 | c9: 0.48 |
| 21:20,45:0 | 626 h, 1,093 h (34 of 48) | c49: 0.57 | c9: 0.49 |
| the four | 2,340 h, 2,452 h | **0.57** | **0.49** |

- r31 is the line with the most left to do by the plain sum at every node (sum of shortfalls 8,072 to 9,787
  against 5,039 to 7,225 for c9 and c49) though its largest child is among the smallest of any line's. Its marginal
  score is the lowest of the nine (0.06 against 0.13 to 0.26): rule 9's shortlist would not hold it. [m]
- Rows r25, r26, r0 price at 0.63 to 0.98 of r31; r13 above it at every node (1.18 to 1.27), c20
  level with it (0.96, 1.05). Six or seven of the
  other eight lines are cheaper than the run's.
- Powers 1 and 2 take c49 at three nodes and c9 at one; powers 4 and 8 take c9 at all four. Here the
  tail-weighted scores choose the cheaper line (0.49 against 0.57); at the truth nodes powers above 2 lose
  (`truth_powers.txt`). Not settled; the gain is in firing at all.
- The slope is the soft spot: open subtrees are still growing and the heaviest are the most open, so k may be
  too small, and a larger k weighs a line's largest child more. The firing's line against r31 had k been
  otherwise: 0.56 at k = 0.015, 0.55 at 0.02, 0.59 at 0.03, 0.66 at 0.04, 0.77 at 0.05, 1.09 at 0.07. The
  cheapest line (c9) stays at 0.46 to 0.62 up to 0.05. No fit in this tree gives a slope above 0.03.

## Does the pricing hold? Three checks

1. **Where the truth is known** (`price_check.txt`): at the tho30 and tho40 truth nodes, fit the law on one
   line's children alone, price the other seven, compare with their solved subtrees. Over 252 (fitted line,
   other line) pairs the predicted ratio is 0.97 of the true one at the median, rms factor 1.33. Of the 88
   pairs priced under 0.7 of the fitted line, 94% are truly cheaper and the true ratio's median is 0.50. Those
   subtrees are hundreds to thousands of nodes; wil50's are far larger. [m]
2. **Whether the law depends on the line** (`line_effect.txt`): the run takes one line at depth 2, but many
   below. Over the children at depths 4, 5 and 6 with closed subtrees (4,390, 4,027 and 2,071, under thousands
   of parents) the work at equal shortfall differs by the parent's line by a standard deviation of 4%, 10% and
   4%; the extremes are x0.89 and x1.16. Children of c9 and c49, the lines in question, sit at x0.90 to x1.04.
   Rows against columns: x0.90 to x0.98. And r31's children are not unusually cheap: at 150 units short a
   closed subtree costs 1.8 h at depth 3 (all r31's), 1.6 at depth 4, 1.2 at depth 5, 0.9 at depth 6. [m]
3. **Children of the other lines, solved** (`d2/kids/`): see below.

## Depth 1, for the record

`lines_price_d1.txt`: at the five depth-1 nodes fired in `../wil50_sb`, the run's r45 is the cheapest line at
three and within 0.83 of it at the other two (c9). The depth-1 exponent is worth at most 8% there.

## What is not shown

- No line but r31 has been solved at these nodes: 0.57 is a price, not a measurement. The checks say the price
  is the right kind of number; an end-to-end run with strong branching at depth 2 is the measurement.
- Four nodes, the heaviest (a quarter of the run's CPU so far at the time of the firings). Whether nodes 200 to 250
  short behave the same is not known.
- What rule 9 would take at depth 2 of wil50 was not run. Rule 10 at depth 2 came to wil50 from sko42, where it
  halved rule 9's nodes on five subtrees (`../branching_lab/ucheap`); it was not measured on wil50.
