# lpbound: the LP bound without the PSD condition (Kurt Anstreicher's suggestion, Nov 20 2026)

Kurt's question: drop R >= 0 from the DNN relaxation and solve the LP. Is it closer to DNN than QPB is, do its duals
help branching, and can PDLP get it fast?

## Answer (in progress: tai60b depths 7-8 and one child-LP check still running)

**No on both counts at the sizes that matter.** The LP is exactly the Adams-Johnson level-1 RLT LP. It is far from
the DNN bound wherever the solver branches: at tai35b's root and depth-2 nodes it is *below QPB* (104-112M against
QPB's 151-182M and DNN's 261-282M), and at tai60b's depth-17 nodes it is within 15M of QPB on either side, with DNN
60-80M higher [m]. Its duals give flat child bounds: at every node tried, no child of any line reaches the cutoff
(the LP is 11-64% short of it, the largest reduced-cost gain at most a third of what is needed), because at the LP
optimum 85-95% of the x_ij are fractional and complementary slackness sets their reduced costs to zero. What is
left is a tie-break whose line rankings are no better than chance against the bounded truth. PDLP does solve it on
this box: m = 33 in 3-4 minutes, m = 43 in 8-24 minutes, under 1 GB, to 1e-4; the cost is not what kills it.

## What the LP is (T0)

The DNN relaxation of the solver (`docs/dnn_bound/dnn_bound.tex` eq. (dnn), cognac) with R >= 0 dropped:
Y = V R V^T with R a free symmetric matrix is the same as range(Y) in range(V), i.e. T Y = 0 for the assignment map
T on y = (1, x). Column (kl) of T Y = 0 says sum_j Y_(ij),(kl) = Y_0,(kl) = x_kl for every i, and
sum_i Y_(ij),(kl) = x_kl for every j; column 0 says X is doubly stochastic; the gangster zeros reduce the i = k
equation to Y_(kl),(kl) = x_kl; symmetry of Y is the complementary constraint y_ijkl = y_klij; 0 <= Y <= 1 gives
y >= 0. That is (1.9) of the thesis (section 1.5.2), the Adams-Johnson LP, with nothing added and nothing missing.
The same identification was already made in `experiments/wil50_feasibility/roadmap/bound_strength.md` (item 2: the
certificate's polyhedral set intersected with Y = Y' is the AJ RLT-1 LP).

Checked numerically [m `t0check.py`, `out/t0check.txt`]: the literal LP (all of Y's upper triangle, T Y = 0,
gangster zeros, the box) and the compact one (`ajlp.py`: one variable per unordered compatible pair, m^2 (m-1)^2 / 2
of them, plus x) agree to the unit on three nodes of nug12, had12 and tai12b (m = 9, 9, 10). On the same nodes every
child's own LP is at least the bound `ajlp.child_bounds` derives for it from the parent's dual (min difference 0.0),
so the child bounds below are valid. Against published values: nug12 522.89 and nug15 1040.99 round up to
RRD95's 523 and 1041 (Resende, Ramakrishnan and Drezner's IPM solve, QAPLIB's lower-bound table).

### Roots and nodes (T0)

QPB is the solver's (`bnbdrv`, `PARAM/default.par`, the node trace's field), DNN is `dnntest -H 0 -i 20000` without
a cutoff (the first solve's certified value, no cut phase; 20,000 iterations or convergence), LP by HiGHS IPM
(exact, small) or PDLP to 1e-4 with its valid dual bound (`t0.py` -> `out/t0_*.tsv`; `t2.py` for the tai nodes).

| instance, node | m | GLB | QPB | **LP** | DNN | opt / cutoff | LP's share of QPB -> DNN |
|---|---|---|---|---|---|---|---|
| nug12 root | 12 | 493 | 479.4 | **522.9** | 568.0 | 578 | 49% |
| nug15 root | 15 | 963 | 990.6 | **1041.0** | 1140.6 | 1150 | 34% |
| nug18 root | 18 | 1554 | 1698.1 | **1663.0** | 1893.5 | 1930 | below QPB |
| had12 root | 12 | 1536 | 1589.6 | **1621.5** | 1652.0 | 1652 | 51% |
| had14 root | 14 | 2492 | 2626.1 | **2666.1** | 2724.0 | 2724 | 41% |
| had16 root | 16 | 3358 | 3589.4 | **3560.2** | 3720.0 | 3720 | below QPB |
| had18 root | 18 | 4776 | 5135.6 | **5087.9** | 5357.5 | 5358 | below QPB |
| tai12b root | 12 | 9.79M | 15.78M | **31.70M** | 39.46M | 39.46M | 68% |
| tai15b root | 15 | 11.24M | -124.5M | **51.56M** | 51.70M | 51.77M | 99.9% |
| tai35b root | 35 | 30.9M | 150.7M | **104.4M** | 261.5M | 283.3M | below QPB |
| tai35b 19:26,18:12 | 33 | | 180.7M | **108.5M** | 280.7M | 283.3M | below QPB |
| tai35b 19:12,18:26 | 33 | | 181.1M | **108.6M** | 282.1M | 283.3M | below QPB |
| tai35b 19:26,18:8 | 33 | | 182.2M | **111.9M** | 281.4M | 283.3M | below QPB |
| tai60b d17 ...22:31 | 43 | | 527.7M | **541.9M** | 608.0M | 608.2M | 18% |
| tai60b d17 ...8:49 | 43 | | 549.5M | **525.0M** | 607.6M | 608.2M | below QPB |

[m; the tai35b and tai60b DNN values are the traces' first-solve bounds, production settings.] The LP sits between
QPB and DNN at the roots of the small nug, had and tai-b instances (for tai15b, where QPB is negative, it nearly
reaches DNN), but from n = 16-18 on it falls below QPB on nug and had, and at the tai35b and tai60b nodes where the
branching decisions are made it is at or below QPB and 60-170M below DNN. So "the LP may be closer to DNN than QPB
is" is true of small instances only. Kurt's stop rule (no closer on every instance) did not fire, because of the
small instances, so T1 and T2 were run.

## Do the duals rank children and lines (T1, T2)

For every child x_ij = 1 of a node, from any row dual u of the parent's LP (all variables boxed in [0,1]):
bound_ij = b'u + d_ij + sum of min(0, d) over the variables still compatible with x_ij = 1 (`ajlp.child_bounds`).
With an exact dual this is the LP + d_ij of the brief (plus nothing); with PDLP's approximate dual it stays valid.
Rule-2 analogue: per line, the count of children whose bound reaches the cutoff, ties by the sum of the gains
(bound - LP); also tried, the sum of the 3 smallest gains. Truth: the bounded children of each line (linetest,
lsla, line_signals traces; `truth.py`).

### The reduced costs are flat (the T2 stop rule)

| node | m | LP short of cutoff | x_ij > 1e-3 at the LP optimum | child gain: median / largest | children reaching the cutoff |
|---|---|---|---|---|---|
| tai35b root | 35 | 180.1M | 1,062 of 1,225 | 0.13M / 6.2M | 0 of 1,225 |
| tai35b 19:26,18:12 | 33 | 174.8M | 1,026 of 1,089 | 80 / 26,686 | 0 |
| tai35b 19:12,18:26 | 33 | 174.7M | 1,022 of 1,089 | 86 / 34,335 | 0 |
| tai35b 19:26,18:8 | 33 | 171.5M | 1,016 of 1,089 | 317 / 40,222 | 0 |
| tai60b d17 22:31 | 43 | 66.3M | (not saved) | 26 / 0.48M | 0 of 1,849 |
| tai60b d17 8:49 | 43 | 83.2M | 1,428 of 1,849 | 2,713 / 26.2M | 0 of 1,849 |

[m `out/t2_*.txt`, PDLP 1e-4.] No child of any line reaches the cutoff at any node, so every line's count is 0 and
the rule is decided entirely by the tie-break. This is the same failure as the DNN certificate's child bounds
(dnn_bound.tex section 3: median 0.1, largest 12.1 against 14.4 needed) and for the same reason, now with a proof:
at an exact LP optimum, d_ij = 0 wherever x_ij > 0, and the LP's x is fractional almost everywhere (tai35b depth 2:
max x_ij 0.12). The small exact solves show it with exact duals [m `tsmall.py`, `out/small_summary.txt`]: at nug12's
root 138 of 144 x_ij are positive; the simplex (vertex) dual gives a positive gain to 1 child (largest 6.6) while
the children's own LPs are tens of units higher; the IPM (central) dual gives 6 children a gain, the largest 0.13.

### The tie-break against the truth

| node | lines with truth (branching children) | LP rule's order of them (rank among all lines) | right? |
|---|---|---|---|
| tai35b root | c26 2, r19 3 | r19 #1, c26 #5 | wrong |
| tai35b 19:26,18:12 | c23 1, c5 2, c6 3, r26 3 | c5 #1, c23 #25, c6 #34, r26 #46 | top pick 2nd best |
| tai35b 19:12,18:26 | c23 1, c5 2 | c5 #1, c23 #43 | wrong |
| tai35b 19:26,18:8 | c5 0, c22 0 | c5 #3, c22 #56 | tie in truth |
| tai60b d17 22:31 | r13 2, c54 3, c26 7 | c54 #2, c26 #10, r13 #80 | wrong (best last) |
| tai60b d17 8:49 | r2 0, r51 5, c30 12 | r2 #2, r51 #41, c30 #85 | right |

[m `out/t2_*_lines.tsv`.] The sum of gains is dominated by a few far children (gains of millions on lines whose
other children gain tens), so it behaves like rule 2's tie-break: it picks rule 2's own line at both depth-17 nodes
where we know it (r2 at 8:49 #2; c54 at 22:31 #2). The sum of the three smallest gains does no better (`*_lines.tsv`).

### Children (T1)

| node | children | Spearman with child DNN: parent LP + gain | child's own QPB | AUC (non-branching above branching): LP bound / QPB |
|---|---|---|---|---|
| tai35b root | 69 | 0.65 | 0.67 | 0.94 / 0.75 |
| tai35b 19:26,18:12 | 129 | -0.10 | 0.59 | 0.44 / 0.46 |
| tai35b 19:12,18:26 | 66 | -0.04 | 0.76 | 0.69 / 0.55 |
| tai35b 19:26,18:8 | 66 | -0.08 | 0.66 | - |
| tai60b d17 22:31 | 127 | 0.14 | 0.81 | 0.69 / 0.87 |
| tai60b d17 8:49 | 113 | 0.56 | 0.71 | 0.52 / 0.84 |

The child's own QPB (which costs a QPB solve per child, as the solver already does) predicts the child's DNN bound
better than the parent's LP dual everywhere but tai35b's root. QPB + U (rule 2's prediction from the parent alone)
was not computed here; the comparison above is against the stronger child QPB.

Small exact nodes (all m^2 children bounded, `out/small_summary.txt`): at nug12's root the IPM dual's gains, though
at most 0.13, rank the children's DNN with Spearman 0.79 and the vertex dual's not at all (-0.00); the PDLP dual
0.16. **Exact and approximate duals disagree on the ranking** (Kurt's two reasons pulling against each other); only
the interior dual carries order, and none carries magnitude. The child's own exact LP ranks them 0.47, its QPB 0.56.

## Cost (T3) (in progress)

| size | LP size | method | wall / CPU | memory | note |
|---|---|---|---|---|---|
| m = 12 | 8.9k vars | IPM / simplex / PDLP 1e-4 | 3 s / 124 s / 0.3 s | <0.1 GB | |
| m = 33-35 | 0.57-0.73M vars | PDLP 1e-4, 2 threads | 130-230 s | <1 GB | box shared |
| m = 43 | 1.63M vars, 155k rows, 6.7M nnz | PDLP 1e-2 / 1e-3 / 1e-4, 4 threads | 18 / 52 / 468-520 s wall | 0.9 GB | 1e-2's dual bound is useless (negative); 1e-3's 120M low |
| m = 43 | | PDLP 1e-4, 3 threads, box shared | 1,421 s wall, 2,482 s CPU | 0.9 GB | 8:49 |

For comparison: the node's DNN first solve plus cut phase at tai60b's depth 17 costs about 600 s (the parents'
CPU field, 252-622 s), so even at the cheapest accuracy that gives a usable dual the LP costs as much as the node's
own DNN bound, and it is a much weaker bound. At tai60b depth 7 the 150-iteration ADMM lookahead costs 5,100 s and
ranks lines right; the LP's duals do not rank them at any cost.

Hahn-Grant: thesis 1.5.3 describes it as a Hungarian-style dual ascent for (1.10), the dual of this same LP, giving
"bounds that are close to the optimal solution of (1.10)" cheaply. The solver already had one: the DNN_CERTX
experiment (commit 594cd3f1, removed Nov 9; `roadmap/certificate/`, `pipeline_ab/`) runs Hahn-Grant skew polishing on
L + Z'; with Z' = 0 it is a Hahn-Grant ascent on exactly this LP. Trying it would mean restoring that code path with
Z = 0. Not done: T2 says the exact duals are not worth approximating.

## What was not tested

- QPB + U per child (rule 2's own prediction); compared against the child's QPB instead.
- tai60b depths 7 and 8 (m = 52-53) and the third depth-17 node: running.
- The children's own LPs at m = 42 (a few, running); at the small nodes they rank children's DNN no better than QPB.
- Adding the LP's cuts to the DNN (the LP as a cut generator) or a stronger LP (RLT-2): out of scope.

## Files

`ajlp.py` (LP builder, HiGHS and PDLP solves, valid child bounds, GLB), `t0check.py` (equivalence and validity),
`t0.py` (T0 rows), `tsmall.py` (exact duals on small nodes), `truth.py` (bounded-children sets from the traces),
`t2.py` (T1/T2 at a truth node), `scale.py`, `childlp.py`, `queue_*.sh`; outputs in `out/`.
