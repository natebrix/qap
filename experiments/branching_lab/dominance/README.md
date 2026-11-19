# Exchange dominance: dropping nodes that hold no optimal solution (Nov 19 2026)

**Answer.** The pair exchange test is correct, exact and almost free: for a node with rows i and k fixed, the
minimum over all its completions of cost(p) - cost(p with i and k exchanged) is one linear assignment problem,
and on every instance the solver runs (one matrix symmetric) that LAP is a sort. If the minimum is >= 1, the node
holds no optimal solution and can be dropped. Replayed on the traces, it would have dropped **24.5% of tai60b's
nodes and 11.5% of its CPU (413 of 3,597 core-hours)** [m]. It drops nothing material on tai50b (1.5% of nodes,
0.04% of CPU), tai50b-r9 (4.7%, 0.4%), tai35b (8%, 0.8%), tho40ab-32 (0.8%, 0.4%) or nug30 (0) [m]. It does
**not** touch the mirror chain `44:58` itself: along it, the (44,43) exchange's minimum stays between -1.45M and
-145k down to depth 46, where the chain dies [m]. What it drops on tai60b is a second kind of near-miss subtree: the
children of a chain node where two rows of the chain swap their columns (rows 2 and 5 at depth 15, 8 and 4 at
depth 11, ...). These come in equal shares from the path chain (27% of its nodes) and the mirror chain (22%).
The one-fixed test adds nothing (0.2% of tai60b's nodes, all leaves). Since the share is above 10% of tai60b's
nodes, this README ends with a specification for `latest/bnb.C`.

Caveats: savings come from a replay with the run's own cutoff (608,215,054.05); a run that finds its incumbent late
saves the same nodes, since the test does not use the incumbent. The saving is 11.5% of tai60b's CPU, not 24.5%,
because the dropped subtrees are deep (their nodes are cheap). On the four other trees it is noise. Only
known optima can be checked on the real trees (the QAPSOLN solution, and its images under the grid automorphisms
for nug30 and tho40), and no node holding one is ever flagged.

Everything is in this directory; per-node replay files (5 to 15 MB) are in the scratchpad (`replay.py` rewrites them).

## 1. The argument

**Notation** (checked against `QAPSOLN` in `check_fixes.py`, `path_check.py`). cost(p) = sum_{u,v} A[u,v] B[p(u),p(v)]
with A the instance file's first matrix. A trace fix `r:c` means p(r) = c. On tai-b the rows index the
distance matrix A, which is symmetric, and B (flows) is asymmetric. On nug30, B is the grid; on tho40, A is.

**Pair exchange.** Let a node N fix i at a and k at b, and let p' be p with p'(i) = b, p'(k) = a. Only the terms
touching i or k change, and with no symmetry assumed:

    Delta(p) = cost(p) - cost(p') = D0 + sum_{l != i,k} D_l(p(l))
    D_l(c) = (A[i,l] - A[k,l]) (B[a,c] - B[b,c]) + (A[l,i] - A[l,k]) (B[c,a] - B[c,b])
    D0     = (A[i,i] - A[k,k]) (B[a,a] - B[b,b]) + (A[i,k] - A[k,i]) (B[a,b] - B[b,a])
             (+ C[i,a] + C[k,b] - C[i,b] - C[k,a] if there is a linear term C)

The fixed rows l contribute constants. So min over the completions of N of Delta is the sum of those constants
and a LAP over the free rows and free columns with cost D_l(c). The sketch was right in substance: its abar = A + A^T
form is this formula when A is symmetric, and the diagonal and asymmetric terms above are what it left out.
**If either matrix is symmetric, the LAP matrix has rank one.** With A symmetric, D_l(c) = x_l y_c, where
x_l = A[i,l] - A[k,l] and y_c = (B[a,c] - B[b,c]) + (B[c,a] - B[c,b]). The minimum is then sum x_(j) y_[j]: x sorted up
against y sorted down (rearrangement inequality). With B symmetric, the roles swap. No LAP is needed.

**One-fixed form.** k is fixed at a and i is free. For each free column c, the completions with p(i) = c form the
node N + {i: c}, where the pair formula applies to (i at c, k at a). N is dominated by (k, i) if that minimum is >= 1
for every c: one sort per c (`dom.py`, `onefixed_fast`).

**What is claimed.** If min Delta >= 1, then every p in N has an exchange image p' with cost(p') < cost(p). So no p in N
is optimal: N contains no optimal solution (strict dominance). Nothing is said about where p' is, so **the twin
need not be a node of the tree.** It may already be fathomed, dropped, or never created. The argument only needs
p' to exist as a permutation.

**Why the proof stays valid next to bound fathoming.** The invariant a branch-and-bound proof keeps is: *either the
incumbent is optimal, or some open node contains an optimal solution*.
- Fathoming N by a bound (bound >= cutoff = incumbent - DELTA) keeps it. If N held an optimum p*, then
  cost(p*) >= bound > incumbent - 1, so the incumbent is optimal.
- Strict dominance keeps it trivially, because it removes only nodes that contain no optimum. This holds whatever the
  incumbent is and whether or not it is optimal.
- So the run ends with an optimal incumbent, and "no solution below the cutoff" claims still hold: if no solution
  was found below a cutoff L, then the global minimum, which lies in no dominated node, lies in an explored or
  bound-fathomed one, so it is >= L. What is lost is only the enumeration of *all* solutions below the cutoff, which
  the solver does not claim.
- Monotonicity: a descendant of a dominated node is dominated by the same pair, since it has fewer completions and
  the same fixes. So testing a node when it is popped is enough. On every replay, the nodes flagged equal the nodes
  under the highest flagged ones [m].

**Ties and non-strict dominance.** With min Delta >= 0, two twins can drop each other. esc08a's nodes {0:0, 1:1} and
{0:1, 1:0} both have best min Delta = 0, both hold optima, and a non-strict test drops both (`nonstrict.py`,
`nonstrict.out`: 336 of 392 small esc08a nodes holding an optimum have best min exactly 0 [m]). A non-strict test is
valid only with a strict total order to break ties. Lexicographic order is one: for the pair (i < k), p' is
lexicographically smaller iff b < a, a property of the node. That order can conflict with the grid symmetry
reduction, which keeps an orbit representative chosen by other rules: the symmetry reduction may keep only the
lexicographically larger image, and the tie-break test then drops it. So: strict only.

**Integrality.** The tested instances have integer data, so Delta is an integer and "strict" is ">= 1". With A
and B both symmetric and zero diagonals (nug, tho, esc), every cost is even, so Delta is even and ">= 1" is the same
as ">= 2". This changes nothing, since no Delta lies between 0 and 2. For tai-b (B asymmetric) costs are not
always even, so the test is ">= 1". In floating point, the solver should evaluate the test exactly. The data are
integers or halves after `symmetrize_QAP`, and the sums are < 2^53. It then tests min > 0.5.

**Grid symmetry reduction (nug, tho, wil).** The solver keeps one child per orbit of the node's stabiliser. That keeps
the invariant: if N holds an optimum p*, a kept child holds an image g p*, and g p* is optimal because the automorphism
preserves cost. Strict dominance never drops a node holding *any* optimum, so it never drops the one holding g p*.
The two compose in any order. This holds for strict dominance only (see ties).

**QPB child elimination.** branch() drops children whose predicted bound (parent bound + U_ij) is above the cutoff.
That is a bound fathom and keeps the invariant as above.

## 2. Validation by brute force [m]

`validate.py` → `validate.out` covers 37 instances:
- nug05-08, esc08a-f (up to 17,280 optima), law07, tai10a, tai10b (B asymmetric), had10, scr10, rou10,
  lipa10a/b (A asymmetric), car10ga/pa: all n! permutations.
- tai12a/b, nug12, had12, scr12, chr12a/b: completions of three 3-fix nodes each, with the optimum value from QAPSOLN.
- 12 random n = 7 or 8 instances with entries in {0, 1, 2}: both matrices asymmetric, or one symmetric, all with
  nonzero diagonals, many ties and up to 13 optima.

The nodes were drawn half at random and half from near-optimal permutations, a quarter of them with a mirror-like
exchange. For each node it checks:
- **Exactness**: for every pair, the sort form and the general LAP form each equal the brute-force minimum of
  cost(p) - cost(swap(p)) over all completions. For the one-fixed test, its value equals the brute-force maximum
  over (k, i) of the minimum over completions.
- **The claim**: for every flagged node, min cost over N > optimum value. This covers every optimal solution without
  listing them.
- **Stronger**: every completion's exchange image is strictly cheaper.

| nodes | pair evaluations | inexact (sort / general) | pair-dominated | counterexamples | images not strictly cheaper | one-fixed evaluations | inexact | one-fixed dominated | counterexamples |
|---|---|---|---|---|---|---|---|---|---|
| 8,100 | 92,304 | 0 / 0 | 18,686 | **0** | 0 | 8,002 | 0 | 2,739 | **0** |

At n = 60:
- `xcheck.py` → `xcheck.out`: the vectorised one-fixed test equals the validated loop version on 200 nodes of tai35b,
  nug30, tho40, tai12b and chr12a.
- `realcheck.py` → `realcheck.out`: on tai60b's flagged nodes, the sort value equals scipy's LAP value, the LAP
  minimiser evaluated directly gives exactly that Delta, and sampled completions (random, and exchange-descent local
  optima inside the node) never fall below it. This covered 80 flagged nodes, 40 of them branched, with 100 random
  and 3 local-search completions each: 0 mismatches and 0 violations (the closest was 554 above the minimum).

## 3. Replay on the real trees [m]

`replay.py` computes, for every trace node, the best pair value (max over pairs of min Delta) and the one-fixed value.
`analyze.py` rebuilds the tree: the fixes are listed in fixing order, so the parent is the prefix. It then takes the
highest flagged node on each branch and adds up the subtrees under them (outputs `analyze_<run>.out`). A flagged
node's own bound is counted as saved.

| run | nodes | CPU h | pair test: nodes saved | CPU saved | branched nodes flagged | one-fixed: nodes saved | known-optimum nodes flagged |
|---|---|---|---|---|---|---|---|
| tai60b-261117 | 331,154 | 3,597 | **81,265 (24.5%)** | **413 h (11.5%)** | 1,672 of 7,496 | 536 (0.16%) | 0 of 58 |
| tai50b-261116 | 8,225 | 132.7 | 124 (1.5%) | 0.06 h (0.04%) | 0 of 214 | 3 | 0 of 48 |
| tai50b-r9 | 6,001 lines | 115.3 | 280 (4.7%) | 0.42 h (0.36%) | 6 of 144 | 1 | 0 of 19 |
| tai35b (line_signals/runs/tai35b_base) | 1,017 | 1.9 | 82 (8.1%) | 0.01 h (0.75%) | 0 of 46 | 4 | 0 of 33 |
| tho40ab-32 (depth-2 subtree) | 22,650 | 21.9 | 174 (0.8%) | 0.10 h (0.44%) | 0 of 644 | 0 | (no optimum in it) |
| nug30-trace-a | 3,688 | 0.6 | 1 | 0 | 0 of 141 | 0 | 0 of 6 (4 optima images) |

- **Validity on real data**: no node containing a known optimum is flagged in any run. That is the QAPSOLN solution,
  and for nug30 and tho40 also its images under the grid automorphisms (`autos.py`). On tai60b that is all 58
  nodes of the best known solution's path.
- On every run but tai60b, the flagged nodes are almost all ones the solver fathomed with its first bound anyway
  (D, Q, P). Only the bound's CPU at those nodes is saved.
- tai50b-r9's trace has 6,001 lines against the 6,629 nodes counted for that run, so the trace may be incomplete.

**tai60b, by depth of the highest flagged node** (pair test; `analyze_tai60b.out`):

| depth | highest flagged nodes | nodes under them | CPU h |
|---|---|---|---|
| 8 | 13 | 13 | 0.05 |
| 11 | 161 | 15,647 | 115.6 |
| 12 | 18 | 993 | 7.6 |
| 13 | 418 | 418 | 2.2 |
| 14 | 542 | 1,487 | 8.6 |
| 15 | 3,135 | 50,975 | 253.1 |
| 16 | 974 | 3,074 | 12.6 |
| 17 | 2,358 | 2,358 | 2.7 |
| 18 | 4,538 | 5,706 | 10.7 |
| 23-54 | 550 | 594 | 0.1 |

Of the 12,707 highest flagged nodes, 114 were branched; the rest are leaves the bound fathomed.

**Where tai60b's savings are** (`tops.py` → `tops_tai60b.out`). They come from 114 branched nodes at depths 11 to
18 that hang directly off the two chains. Each is a chain node's sibling in which two of the chain's rows have
exchanged their columns, and its twin is the chain node itself. Examples:
- The pair (2,5) at depth 15 (min Delta about +3,250 on the path chain, +3,600 on the mirror chain): 140 of these
  highest nodes per chain, 25k and 23k nodes.
- The pair (8,4) at depth 11 (min +6,087 and +7,013): about 230 per chain, 8.7k and 8.2k nodes.

All 50 such nodes with subtrees of 300 or more have their twin in the trace, branched. The test drops 27.2% of
the path chain's subtree (215 h, 13.5% of its CPU) and 22.2% of the mirror chain's (198 h, 9.9%).

**The mirror chain `44:58`** (`mirror.py` → `mirror_tai60b.out`):
- The mirror of the best known solution (rows 43 and 44 exchanged) costs 2,786 more, as the postmortem says.
- The (44,43) exchange's min Delta over the mirror node's completions is -1,453,988 at depth 2, -1.17M at depth 14,
  -315k at depth 26 and -145,064 at depth 46, the chain's last node. So the test never fires on the chain itself.
  Neither does any other pair (`analyze_tai60b.out`, chain listing: best pair -1.45M at depth 2, -27k at depth 46).
- The chain's own off-chain subtrees lose 22% of their nodes, as above.

**Why the mirror survives.** Rows 43 and 44 are an isolated two-location cluster, so x_l = A[43,l] - A[44,l] is small
but not zero. The flows y_c of facilities 45 and 58 differ by a lot. The LAP pairs the largest |x| with the
opposite-sign largest |y|, so Delta goes very negative on completions that are useless otherwise.
`run_restricted.py` → `restricted_tai50b.out`: on tai50b's branched nodes, the LAP minimiser of the best pair costs a
median **57% above the cutoff** [m]. The test is defeated by completions no proof needs to consider.

## 4. Restricting to useful completions and multi-exchange tests [m]

**Restricted test** (`restricted.py`, validated by `validate_restricted.py` → `validate_restricted.out`: 720 nodes,
377 flagged, 0 failures).
- The test is min over completions of Delta(p) + mu (G(p) - UB), with G the Gilmore-Lawler linearisation of cost
  (G <= cost). If it is >= 1 for some mu >= 0, every completion below UB is strictly dominated, and dropping is as
  valid as fathoming.
- On tai50b's 203 branched nodes, with the 8 best pairs and mu in {0 ... 10}: **0 flagged**, and the best mu was
  always 0. GLB is 30-80% below the cutoff at these depths, far too weak to restrict anything.
- The same idea with the solver's own bound would mean bounding the node's QAP plus a linear term t·D for each
  candidate pair: a DNN solve per pair. It interpolates between dominance and the bound and is not cheap. Not tried.

**More than one exchange** (`cycles.py` → `cycles_tai50b.out`). A 3-cycle, or two disjoint exchanges, of fixed rows
still gives a Delta linear in the free rows: one LAP.
- Tested on 80 tai50b branched nodes: every 3-cycle of the 12 rows in the 8 best pairs, and every double exchange
  among those pairs. **0 flagged.** Every such value was below the best single pair; 3-cycles were typically 2 to
  3.5 times more negative.

**One-fixed test**: never fires on a branched node in any run.

## 5. Cost [m]

`replay.py` times each test in Python and numpy, against the solver's CPU per node in the trace:

| tai60b depth | pair test ms | one-fixed ms | solver s per node |
|---|---|---|---|
| 5 | 0.26 | 1.35 | 141 |
| 10 | 0.32 | 2.26 | 76 |
| 15 | 0.32 | 2.42 | 17 |
| 25 | 0.69 | 4.11 | 11 |
| 40 | 2.1 | 5.3 | 0.93 |
| 55 | 9.2 | 4.0 | 0.10 |

Mean over the tree: pair 0.40 ms (max 33 ms) against 39 s per node, a cost of 1e-5 of the run [m]. In C the work per
node is at most d(d-1)/2 pairs × (O(d) constant + an O(m log m) sort), under 10^6 operations at n = 60 [e].

## What was not tested

- The test inside a real solver run. The replay assumes the rest of the tree unchanged. That holds when the cutoff
  is fixed from the start, as in all these runs; with an incumbent found late it saves the same or more.
- tai40b: no trace found. wil50 and sko traces were not replayed. On tho and nug, the dominated nodes are fathomed
  anyway.
- The restricted test with the DNN bound in place of GLB (see 4).
- Lookahead: testing the children of a line for exchange dominance when choosing the line. A line whose children are
  dominated in large numbers would be a good line. Not measured.

## Specification for `latest/bnb.C` (exchange dominance, pair form)

**Where.** In the node loop of `bnb()`, right after the pop-time discard (the `std::max(node.bound,
node.predicted_bound) > s->getObjective() - DELTA` block, about line 566) and before `pick_param` and any bound:

    if (exdom_on && node.nfix >= 2 && exchange_dominated(q, node, &pi, &pk, &val)) {
      stats->anotherDominated(node.depth);           // a new counter, beside anotherDiscard
      if (trace) node_trace(trace, node, 'X', NAN, NAN, NAN, s->getObjective() - DELTA, stats->ier_start_time);
      continue;
    }

Testing at pop covers the children of every branching and nodes restored from a state file. Because of
monotonicity, a dominated node's descendants are never created. Optionally, run the same test on `kid[]` after
`branch()` and before `push`, so that dominated children never take stack or state-file space.

**Data.**
- `q->A`, `q->B` (and `q->C` if nonzero) as the objective uses them, after `symmetrize_QAP`, so both are symmetric.
  The symmetrised matrix holds halves; the objective is unchanged, and so is Delta.
- The node's fixes: rows `node.history[0..nfix-1]`, columns `node.p[row]`, in original indices, as `node_trace`
  prints them.
- The free rows and the free columns.
- If any perturbation of the objective is on (EVB3), use the unperturbed objective, or disable the test.

**Computation** (A and B symmetric). For each pair of fixed rows (i at a, k at b):

    const = (A[i][i]-A[k][k])*(B[a][a]-B[b][b]) + (C[i][a]+C[k][b]-C[i][b]-C[k][a])
          + sum over fixed l != i,k of 2*(A[i][l]-A[k][l])*(B[a][p(l)]-B[b][p(l)])
    x_l = A[i][l]-A[k][l] for free rows l, sorted ascending
    y_c = 2*(B[a][c]-B[b][c]) for free columns c, sorted descending
    min = const + sum_j x[j]*y[j]
    dominated if min > 0.5 (integer data; 0 if values are not integers, or skip the test then)

Accumulate in long double or int64 on the doubled values, exactly as `LowerSum`. Stop at the first dominated pair.
Integrality: `set_fathom_tolerance` already knows whether objective values are integers. If they are not, use
min > 0 with an exact or safely rounded sum.

**Cost.** At most nfix(nfix-1)/2 pairs, each O(nfix + m log m). The pair's y depends only on (a, b) and x only on
(i, k). Under a millisecond per node at n = 60 [e], against 10 to 140 s of bound per node at tai60b's heavy depths.

**Validation for the implementer.**
1. A debug mode that writes each popped node's best pair value. Compare it with `replay.py`'s output on the
   tai35b and tai50b traces: the values must match exactly.
2. A real subtree. The path chain's depth-11 node of tai60b
   (`44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,24:2,8:29,4:23`, 21,745 nodes, 86.6 CPU h in the run) should lose
   9,371 nodes and 42.4 CPU h [m, replay]. Its depth-15 node (`...,2:49,5:57,9:55,12:0`, 12,726 nodes, 41.4 h)
   should lose 1,417 nodes and 5.3 h.
3. In both, the subtree's result (no solution below the cutoff, or the same optimum) must be unchanged.

Trace outcome `X` is new: `lib/qap/` on the website and `coord/` parsers need to accept it.
