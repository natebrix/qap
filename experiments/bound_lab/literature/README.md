# Literature sweep 2021-2026 for the bound lab (agent `literature`, Nov 15 2026)

Question (`../COMMON.md`): is there anything published since 2021 that would materially
strengthen the certified DNN bound or speed up its computation for tho40, wil50, sko42,
sko49 and tai50b? No solver was run. Earlier surveys in this repo
(`wil50_feasibility/roadmap/admm/alternatives.md`, `new_cuts/literature/REPORT.md` and
`NOTES.md`, `solvers/research_*.md`, `stronger_bounds/research_bound_quality.md`) worked
from search snippets; this pass read the papers.

Marks: **[primary]** read in this session in the paper's PDF or its arXiv abstract;
**[snippet]** from a search-result summary only; **[repo]** taken from this repo's earlier
notes and not re-read; **[m]** computed here (`structure.py`, `structure.log`);
**[e]** my estimate.

Access: arxiv.org PDFs and the arXiv API, optimization-online.org PDFs, the Semantic
Scholar API and one NHR project page were readable through the shell. **Not reachable:**
`qaplib.mgi.polymtl.ca` (the current QAPLIB site: no DNS from this session; its figures
below are [snippet]), web.archive.org, the OpenAlex API (shared anonymous budget used up),
ScienceDirect PDFs. The coral.ise.lehigh.edu QAPLIB mirror is readable but frozen at April
2012. To close the QAPLIB gap, allow `qaplib.mgi.polymtl.ca` in the session's network
settings or open its `news.html` by hand.

## Bottom line

1. **Nothing published solves tho40, sko49, wil50 or tai50b.** sko42 is the exception:
   Fujii et al. solved it in 2021 (6.0e6 nodes, about 737,000 core-hours). Solving tho40
   would be a first as far as any reachable source shows. A certified root bound above
   226,522 (tho40), 22,654 (sko49), 48,126 (wil50) or 431,245,612 (tai50b) would be a
   new best published lower bound.
2. **Cuts: the regional-triangle idea was not found in print**, and the partition
   literature holds nothing stronger than what is closed here. What the literature does
   add is a prediction that fits today's tho40 gate test: triangle-type cuts on a
   DNN partition relaxation gain much more on sparse graphs than on dense ones.
3. **No relaxation stronger than the DNN has been computed on a QAP at n >= 30** since
   Graham et al. 2022. The 2024 exactness papers are theory with n = 10 experiments.
   One older fact matters for tai50b and was misread before: on tai-b, RLT2 beats the DNN
   up to n = 30 (tai30b: 2.62% gap against 6.04%), so that family has 3-wise slack the
   DNN misses (Q2c).
4. **No solver beats a warm-started PRSM on a high-rank DNN.** The low-rank solvers'
   own authors say QAP is the wrong class for them. The newest GPU work on first-order
   SDP (polynomial-filter projections with warm-started subspace tracking, in half
   precision) is the design `latest/dnn.C` already has on the CPU.
5. **Branch and bound: one transferable item**, Fujii et al.'s estimates by branching
   rule, which differ by 25x to 19,000x between rules on the same instance, in opposite
   directions on tho40 and tai50b.

## Q1. Status of the open instances

**Short answer.** As of the newest reachable sources (an NHR project description dated
April 15 2025 and searches through November 2026): tho40, sko49, sko56, wil50, tai35a,
tai40a, tai50a, tai50b, tai60a/b, tai80a/b, tho150 and wil100 have no published proof.
sko42 and tai30a were solved by Fujii, Ito, Kim, Kojima, Shinano and Toh (2021), tai35b
and tai40b by Gonçalves, Pessoa et al. (2015, RLT2 on GPUs) and again by Fujii et al.
The same group has been running tai50b: the April 2025 text reports "a significant
increase in the lower bounds" and expects that "dedicating more time ... will likely
reach ... the exact solution", with no numbers. No follow-up paper exists on arXiv, and
the 2021 report's citing papers (Semantic Scholar, 5) contain no exact QAP result. The
group's later papers are on tai256c (1.25% gap, not solved).

Fujii et al.'s solved runs (their Table 4; ISM HPE SGI 8600, Xeon Gold 6154, branching
rule M, Newton-bracketing on the Lagrangian DNN) [primary]. Core-hours are wall time
times cores, computed here.

| instance | nodes | wall s | cores | core-hours | core-s per node |
|---|---|---|---|---|---|
| nug30 | 26,181 | 3.14e3 | 1,728 | 1.5e3 | 207 |
| tai30a | 34,000,579 | 5.81e5 | 1,728 | 2.8e5 | 30 |
| tai35b | 2,620,547 | 2.49e5 | 1,728 | 1.2e5 | 164 |
| tai40b | 278,465 | 1.05e5 | 1,728 | 5.0e4 | 652 |
| sko42 | 6,019,419 | 5.12e5 | 5,184 | 7.4e5 | 441 |

Their node estimates for the open instances (their Table 5; a sampling estimator that
they found to be within 3x of the generated count on the solved ones) [primary]:

| instance | rule M | rule P | rule D | their reading |
|---|---|---|---|---|
| tho40 | 7.1e11 | 3.7e7 | 9.8e7 | "could be solved" with P |
| wil50 | 3.3e7 | 2.1e7 | 3.9e7 | "could be solved" with P |
| sko49 | 1.6e8 | 1.3e9 | 1.2e9 | "we may challenge sko49" with M |
| tai50b | 5.5e7 | 1.4e9 | 2.1e11 | "could be solved" with M |
| tai35a | 2.1e11 | 3.9e10 | 1.3e10 | "too difficult" |
| sko42 | 3.6e6 (6.0e6 generated) | 3.1e6 | 2.8e6 | solved with M |

"Could be solved" is their projection from these estimates, not a result; a search
engine summary misreads it as solved. Their tho40 estimate with rule P (3.7e7) is close
to this repo's 3.4e7.

Best published lower bounds (all are the same DNN relaxation at the root; none has cuts):

| instance | best known | Graham et al. 2022 rPRSM [primary] | Kim-Kojima-Toh Newton-bracketing [primary] | BBCPOP / SDPNAL+ (Ito et al. 2019) [primary] | QAPLIB current page [snippet] | gap of the best |
|---|---|---|---|---|---|---|
| tho40 | 240,516 | 226,522 | 226,516 | 226,490 / 226,482 | 226,490 (BBCPOP) | 5.82% |
| sko42 (solved) | 15,812 | 15,336 | 15,334 | - | - | 3.01% |
| sko49 | 23,386 | 22,654 | 22,653 | 22,650 / 22,650 | 22,650 (BBCPOP) | 3.13% |
| wil50 | 48,816 | 48,126 | 48,123 | 48,121 / 48,120 | 48,121 (BBCPOP) | 1.41% |
| tai50b | 458,821,517 | - | 431,245,612 | 431,090,745 / 431,074,160 | 431,090,700 (BBCPOP) | 6.01% |
| tai35a / tai40a | 2,422,002 / 3,139,370 | 2,216,648 / 2,843,314 | - | 2,216,540 / 2,843,199 | - | 8.5% / 9.4% |
| tai50a / tai60a | 4,938,796 / 7,205,962 | 4,390,982 / 6,326,350 | - / 6,326,094 | 4,390,743 / - | - | 11.1% / 12.2% |
| tai60b / tai80a / tai80b | 608,215,054 / 13,499,184 / 818,415,043 | - | 592,829,061 / 11,657,055 / 786,900,529 | - | - | 2.5% / 13.6% / 3.9% |
| sko56 | 34,458 | 33,390 | 33,389 | - | - | 3.10% |
| tho150, wil100 | 8,133,398 / 273,038 | - | - | not in the arXiv version | not seen | 2012 mirror: 7,620,628 (TDB) / 264,442 (SDRMS-SUM) |

Other bounds are weaker at these sizes: RLT2 dual ascent on 95 GPUs gives tho40 213,372
and sko42 14,741 (Date-Nagi 2019, Table 5) [primary]; Mittelmann-Peng's matrix-splitting
SDP gives wil50 46,526, sko49 21,449, sko42 14,393 [primary]. Gueye-Michelon 2015 (an LP
in "distance variables" for Manhattan instances) reports improving "seven best lower
bounds of QAPLIB instances for which no optimality proofs exist" [primary, abstract
only; the full text was not reachable, and the 2012 list it improved on was far below
the DNN on sko100/wil100/tho150, so it may still hold the record there].

**For us.** (a) sko42 is a re-solve, at roughly 1,500 times less CPU than the first proof
if the 500 core-hour estimate holds. (b) tho40 would be a first. (c) Any certified
whole-instance bound from the cut phase above the Graham et al. column is a new record,
and it costs one root run per instance (idea 4 below). The published values are ceilings
of dual bounds from inexact iterates, not interval-certified; ours would be the first
certified ones as well.

| source | year | id | verified |
|---|---|---|---|
| Fujii, Ito, Kim, Kojima, Shinano, Toh, *Solving challenging large scale QAPs* | 2021 | arXiv 2101.09629, ZIB 21-02 | [primary] Tables 4, 5, section 5, remarks (i)-(v) |
| Shinano, *A massively parallel solver for QAPs*, NHR project bem00050 | PDF dated 2025-04-15 | zulassung.hlrn.de/kurzbeschreibungen/bem00050.pdf | [primary] tai50b in progress, no numbers |
| Kim, Kojima, Toh, *A Newton-bracketing method for a simple conic optimization problem* | 2019/2021 | arXiv 1905.12840 | [primary] Table 2 |
| Ito, Kim, Kojima, Takeda, Toh, *BBCPOP* | 2018/2019 | arXiv 1804.00761 | [primary] QAP table |
| Graham, Hu, Im, Li, Wolkowicz, rPRSM | 2020/2022 | arXiv 2006.01529 | [primary] Tables 4.2, 4.3 |
| Date, Nagi, RLT2 on GPU clusters | 2017/2019 | arXiv 1710.03732 | [primary] Tables 5-7 (tai-b rows: see Q2c) |
| Mittelmann, Peng, Hamming and Manhattan distance matrices | 2010 | SIAM J. Optim. 20(6); optimization-online 2008/05/1980 | [primary] |
| Gueye, Michelon, O(n^2)-variable linear formulation, Manhattan | 2015 | doi 10.1007/s13675-014-0033-4 | abstract [primary], results [not read] |
| Fujii, Kim, Kojima, Mittelmann, Shinano, tai256c | 2022-2024 | arXiv 2210.15963, 2401.09439 | [primary] abstract, section 4 |
| QAPLIB, current site | - | qaplib.mgi.polymtl.ca | [snippet] only; 2012 mirror [primary] |
| Gonçalves, Pessoa et al., RLT2 on CPUs and GPUs | 2015 | arXiv 1510.02065 | as cited by Fujii et al. [primary there] |

## Q2. Cuts

**(a) Has the idea been published?** Not found. No paper imposes Boolean-quadric or
triangle inequalities on region-membership moments, taken from a row/column threshold
decomposition of a grid distance, inside the QAP SDP/DNN relaxation. Searches: the arXiv
API (QAP with semidefinite, DNN, cuts, grid, Manhattan, linear arrangement), web
searches, and citation passes over Oliveira-Wolkowicz-Xu 2018 and de Meijer et al. 2023.
The nearest prior art, none of which is the same thing:

- Mittelmann-Peng 2010 use the grid's structure through a different fact: alpha E - X B X'
  is positive semidefinite for every permutation, which gives a small "matrix splitting"
  SDP. No thresholds, no cuts, and a bound 3% below the DNN (wil50 46,526) [primary].
  This repo already notes that the fact is implied by the lifted PSD block [repo].
- Gueye-Michelon 2015: an LP on distance variables with metric inequalities for
  Manhattan instances [primary, abstract].
- The cut decomposition itself is classical for one-dimensional layout: minimum linear
  arrangement and cutwidth are sums or maxima of threshold cuts. SDP work there lifts
  ordering variables, not assignments: Gaar, Puges and Wiegele 2023 on cutwidth
  (3-dicycle equations, triangle inequalities, inequalities from the squared linear
  ordering polytope; final gaps 20% to 63%) [primary]; Rendl, Sotirov and Truden on
  bandwidth through vertex partitions [primary, abstract]; Hungerländer and Rendl 2013
  on single-row layout (matrix of order C(n,2)+1, optimal up to 42 departments)
  [snippet].
- The DNN compressed onto one region is the vector-lifting DNN relaxation of a graph
  bisection, and the regional triangles are the "BQP inequalities" de Meijer, Sotirov,
  Wiegele and Zhao add to it by CP-ADMM [primary]. That paper is the closest relative,
  for graph partition only.

So a write-up can call the construction new with the usual caution; the related-work
section is the four items above.

**(b) Families stronger than triangles that are not on the closed list.** None found.
What each line of work measured:

| family | source | measured strength | status here |
|---|---|---|---|
| triangle / BQP inequalities on DNN partition relaxations | de Meijer et al. 2023, Table 1 | 2-equipartition, DNN to DNN+cuts against the upper bound, with the share of that gap closed (computed here): G124,2.5 (density 2%) 7.29 to 12.01 of 13, 83%; G124,10 (8%) 152.9 to 170.9 of 178, 72%; G124,20 (16%) 418.7 to 440.0 of 449, 70%; G250,2.5 (1%) 15.2 to 28.3 of 29, 95%; G250,5 (2%) 81.5 to 105.0 of 114, 72%; G250,10 (4%) 303.0 to 330.4 of 357, 51%. "The improvement ... gets more significant as the graphs get sparser" [primary]. All their graphs are at most 16% dense | in hand (the cut phase) |
| independent-set (pigeonhole, Chopra-Rao clique) inequalities for k-equipartition | same, section 5 | "the triangle inequalities contribute more to the strength of the bound than the independent set inequalities" [primary] | closed Nov 6 (a sum of line boxes) [repo] |
| bisection knapsack tree / knapsack walk inequalities | Armbruster, Fügenschuh, Helmberg, Martin 2012 | with the SDP relaxation "the best results are obtained when the cycle inequalities are separated exclusively"; knapsack trees improve the SDP bound but less than cycles and at higher cost; knapsack walks "seem not to help at all" [primary] | not worth building |
| pentagonal, heptagonal (5- and 7-clique) | BiqBin 2022, MADAM 2021 | used after triangles saturate (10n triangles + 300 + 200 per round) [primary]; no isolated gain reported | 5-point and local-7 closed on wil50 [repo] |
| exact subgraph constraints, order 5 and 7 | Gaar-Rendl 2019/2020, Table 2 | max-cut, deviation from the optimum in % with orders 3 / 5 / 7: density 10%: 0.40 / 0.00 / 0.00 and 0.13 / 0.00 / 0.00; density 50%: 0.90 / 0.51 / 0.39 and 3.91 / 1.41 / 0.85; density 90%: 0.58 / 0.38 / 0.31 and 8.06 / 5.66 / 5.09 [primary] | order 7 on one region = local-7, closed on wil50 [repo] |
| box-cone grouping of bound-factor RLT cuts | Charbonneau, Hu, Sotirov 2026 | max-cut triangles only inside MADAM's tree: median whole-search time 6.5% lower on 199 graphs, median node count 2.33x higher; "should not be interpreted as certified" [primary] | read Nov 7; our families lack the shared-slack form [repo] |
| Chvátal-Gomory cuts for integer SDPs; Lagrangian duals of mixed-integer SDPs | de Meijer-Sotirov 2024, 2025 | QTSP, max-k-cut; no QAP [primary, abstracts] | the single-threshold integer subproblem is already saturated on wil50 [repo] |
| squared linear ordering polytope inequalities | Gaar, Puges, Wiegele 2023 | cutwidth gaps stay at 20-63% [primary] | need an ordering lift we do not form |
| RLT against disjunctive cuts | Hof, Walter 2025 | RLT dominates disjunctive programming on cardinality equations with right-hand side 1, "applications in the strength comparison for the QAP" [primary, abstract] | agrees with the Nov 4 verdict on disjunctive cuts [repo] |

The one pattern that transfers is the density dependence, which appears independently in
de Meijer et al. (triangles: 95% of the gap at 1% density down to 51-70% at 4-16%) and in
Gaar-Rendl (orders 5 and 7 on top of triangles: all of the remaining gap at 10% density,
43-78% at 50%, 30-47% at 90%). wil50's flows are 90% dense, sko49's 69%, tho40's 40% [m],
which is denser than anything in de Meijer et al.'s tables. This is the literature's
version of today's gate test, and it is the reason the wil50 verdicts on the 5-point and
local-7 families should not be assumed for tho40 (idea 6).

**(c) Clustered distances and sparse or heavy-tailed flows.** No paper exploits either
for a QAP bound at this size. Bravo Ferreira, Khoo and Singer 2018 build an SDP with one
PSD block of order O(n) per edge of the sparse graph; it is strictly weaker than the full
SDP by construction ("we sacrifice positive semidefiniteness of Q") and is tested on chr
and esc [primary].

**RLT2 beats the DNN on tai-b up to n = 30, which it does on no other family at that
size.** Date and Nagi's dual ascent (2,000 iterations) closes tai20b and tai25b at the
root and leaves 2.62% on tai30b (620,444,000), where the DNN leaves 6.04% (BBCPOP
598,629,428; SDPNAL+ 598,979,464) [primary, both tables]. At n = 35 and 40 the same
2,000 iterations fall behind the DNN (tai35b 5.79% against 4.79%, tai40b 5.69% against
4.39%) [primary], which may be the dual ascent not converging rather than the relaxation.
The earlier survey's rule "RLT2 is below the DNN from n = 24" [repo] holds for nug, tai-a,
tho and sko, not for tai-b. So on tai-b there is 3-wise consistency that the DNN misses
and that is worth more than half of its gap at n = 30. The triple families closed on
wil50 (violations of 0.5-4% of mass) [repo] are the local form of exactly that, and the
cluster triangles are its aggregate.

What the data says [m, `structure.log`; the `regions` agent's `logs/structure.txt` has
the same clustering for all tai-b]:

- tai50b's distances are 5 clusters of 16, 15, 12, 4 and 3 locations: mean distance 71
  inside a cluster (max 152), 1,070 between clusters (min 356), and the first MDS axis
  carries 97% of the Gram matrix, so the clusters lie nearly on a line.
- In the best known solution **94.6% of the objective is between clusters**; the
  within-cluster pairs are 24.5% of the pairs, carry 51.5% of the flow and 5.4% of the
  cost. The root DNN gap is 6.0%.
- The flows are 56% zeros and heavy-tailed: the top 10% of the nonzero entries carry 61%
  of the weight (19-26% on tho40, sko and wil).

So tai50b is close to a 5-part graph partition problem with part sizes 16/15/12/4/3 on a
sparse heavy-tailed graph, which is the regime where the partition literature's cuts are
strongest. The transfer: regions = the clusters and their unions (prefixes along the
first MDS axis are the analogue of grid thresholds), the same Boolean triangles on
(p, W), separation weighted towards the heavy flow pairs. The `regions` agent is testing
exactly this; idea 1 below adds what the literature says to try in which order.

| source | year | id | verified |
|---|---|---|---|
| de Meijer, Sotirov, Wiegele, Zhao, *Partitioning through projections* | 2022/2023 | arXiv 2205.06788 | [primary] sections 4, 5, Table 1 |
| Armbruster, Fügenschuh, Helmberg, Martin, LP and SDP branch-and-cut for bisection | 2011/2012 | optimization-online 2011/03/2952; MPC 4 | [primary] section 5.2, 5.4 |
| Gaar, Rendl, exact subgraph constraints (bundle) | 2019 | arXiv 1902.05345 | [primary] Tables 2, 3 |
| Gusmeroli et al., BiqBin; Hrga, Povh, MADAM | 2020-2022 | arXiv 2009.06240, 2010.07839 | [primary] |
| Charbonneau, Hu, Sotirov, box cones | 2026 | arXiv 2610.41051 | [primary] |
| Gaar, Puges, Wiegele, cutwidth | 2023 | arXiv 2301.03900 | [primary] abstract, section 4 |
| Rendl, Sotirov, Truden, bandwidth | 2019/2021 | arXiv 1904.06715 | [primary] abstract |
| Bravo Ferreira, Khoo, Singer, QAP with a sparse graph | 2017/2018 | arXiv 1703.09339 | [primary] abstract, section 2 |
| Hof, Walter, geometry of RLT | 2025 | arXiv 2511.13805 | [primary] abstract |
| de Meijer, Sotirov, CG cuts; Lagrangian duality for MISDP | 2022, 2025 | arXiv 2201.10224, 2501.11397 | [primary] abstracts |
| Hungerländer, Rendl, single-row layout; Amaral, Letchford | 2013 | COAP 55; Math. Prog. | [snippet] |
| Sotirov, *Graph bisection revisited*; Wiegele, Zhao, extended ADMM for graph partition | 2018, 2022 | Ann. OR 265; arXiv 2105.09075 | [repo], not re-read |

## Q3. Relaxations stronger than the facially reduced DNN

**Short answer: nothing new with a computed QAP bound at n >= 30.** The earlier survey's
conclusion stands (RLT3 stops at n = 25; RLT2 is below the DNN from n = 24), with the
tai-b exception of Q2c: RLT2 is above the DNN on tai20b, tai25b and tai30b.

| line | what it is | QAP numbers | verdict |
|---|---|---|---|
| Ling 2024; Chen, Soh 2024 | sufficient conditions for the SDP relaxation to be exact (signal plus small noise; linear inequalities on the data) | experiments at n = 10 only, no QAPLIB [primary] | theory; QAPLIB instances are far from the regime |
| Hou, Tang, Toh, RiNNAL+ 2025 | proof that the DNN relaxation of order n+l+1 equals SDP-RLT of order n+1 | QAP only as a reformulation example (RNNAL: chr, esc32, kra30 at solver tolerance) [primary] | the same relaxation |
| Sinjorgo, Sotirov, Vera 2025 | an intermediate Lasserre level for the stability number: degree-2 monomials chosen from the level-1 solution, matrix order up to 2,500, solved by ADMM | stable set only; "best-known bounds" there [primary] | the QAP analogue is a block on products of assignments for chosen facility tuples: the per-region moment block and the triples, both closed on wil50 [repo]; on tai50b the selection rule would be the heavy flow pairs (idea 1) |
| de Meijer, Sotirov 2023, 2025 | compact mixed-integer SDP formulation of the QAP; Lagrangian dual hierarchy | max-k-cut only [primary, abstracts] | no QAP bound |
| Hu, Xu 2026 (MIQCR-FR) | exact convex reformulation from a facially reduced SDP, then Gurobi | 431 small instances of six families [primary, abstract] | a different route to exactness; small n |
| Kim, Kojima 2021; Qu, Sun, Xu 2025 | cuts proved by a conic solve; DNN cuts at KKT points of continuous QPs | QUBO n <= 250; continuous QPs [primary, abstracts] | closed Nov 4 [repo] |
| Hu, Sotirov, Wolkowicz 2023 | facial plus symmetry reduction | symmetric instances only | no help below the root [repo] |
| Massoulié et al. 2025-26; Valdivia, Tyagi 2023 | convex relaxations for random graph alignment | random models | not applicable |

| source | id | verified |
|---|---|---|
| Ling, *On the exactness of SDP relaxation for QAP* | arXiv 2408.05942 | [primary] |
| Chen, Soh, *Exactness conditions for semidefinite relaxations of the QAP* | arXiv 2409.08802 | [primary] |
| Hou, Tang, Toh, RNNAL, RiNNAL+ | arXiv 2502.13849, 2507.13776 | [primary] |
| Sinjorgo, Sotirov, Vera | arXiv 2506.08648 | [primary] |
| de Meijer, Sotirov, *On integrality in SDP for discrete optimization* | arXiv 2306.09865 | [primary] abstract |
| Hu, Xu, MIQCR-FR | arXiv 2610.11154 | [primary] abstract |
| Qu, Sun, Xu | arXiv 2510.02948 | [primary] abstract |

## Q4. Solvers for large DNN relaxations with certified dual bounds

**Short answer: nothing beats 100-300 ms x 20-600 warm iterations.** Published costs
for this relaxation, cold to 1e-5, per iteration (time over iterations, computed here):

| method | tho40 | wil50 | sko49 |
|---|---|---|---|
| rPRSM, MATLAB (Graham et al.) [primary] | 6,200 it, 2,323 s, 0.37 s/it | 4,700 it, 6,133 s, 1.3 s/it | 7,900 it, 14,235 s, 1.8 s/it |
| Newton-bracketing (Kim et al.) [primary] | 10,065 APG it, 5.0e3 s, 0.50 s/it | 8,650, 1.3e4 s, 1.5 s/it | 9,170, 1.2e4 s, 1.3 s/it |
| BBCPOP (Ito et al.) [primary] | 4,826 it, 1.8e4 s | 5,453 it, 6.7e4 s | 3,383 it, 3.6e4 s |
| SDPNAL+ (same table) [primary] | 1.4e4 s | 7.0e4 s | 5.9e4 s |

| candidate | measured, closest class | could it beat the PRSM here? |
|---|---|---|
| RNNAL / RiNNAL+ (low-rank Riemannian ALM) | the authors: "the solutions to the QAP instances are typically of high rank. Consequently, SDPNAL+ is anticipated to exhibit superior computational speed"; esc32a (order 1,024) 1,630 s against SDPNAL+ 309 s; kra30a (order 900) 974 s against 252 s. On a QMST DNN of order 1,225 whose solution has rank near n/2: RiNNAL+ 484 s, SDPNAL+ 1,987 s, **PRSM 440 s** [primary] | no; at high rank a PRSM is level with RiNNAL+ and 4.5x faster than SDPNAL+ in their own table |
| SDPF+ (Tang, Toh 2026) | "outperform other solvers on large-scale SDPs whose optimal solutions exhibit low-rank structure" [primary, abstract] | no (rank 0.7 K) |
| cuHALLaR, LoRADS / cuLoRADS | Burer-Monteiro on GPUs; stable set, matrix completion, phase retrieval, max-cut at 1e6-1e8 [primary, abstracts]; no DNN, no nonnegativity | no |
| McADMM (Lian et al. 2026) | multi-clique decomposition for sparse aggregate patterns [primary, abstract] | no; the QAP DNN is dense, chordal splitting closed [repo] |
| PDHG for SDP (Ding, Yang 2026; Prakhya, Wang 2026; Jiang 2026) | SDPLIB-scale examples; projections by Frank-Wolfe or Lanczos [primary, abstracts] | no; the PSD projection is the same cost or inexact, no DNN results |
| optimal-metric ADMM (Ran, Vlaski, Dai 2024); Halpern-accelerated preconditioned ADMM (Sun et al. 2024) | QCQP relaxations; convex QP [primary, abstracts] | unlikely; equilibration and Halpern were measured negative here [repo] |
| factorization-free PSD projection (Kang, Han, Groudiev, Yang 2025) | composite polynomial filter, 22 GEMMs, relative error 1e-3 in half precision, about 10x faster than cuSOLVER's eigendecomposition; a 10,000 x 10,000 projection in 55 ms (fp16) or 400 ms (fp32) on a B200 [primary] | this is our Chebyshev-filtered projection on a GPU; it supports the GPU port and offers nothing for the CPU code |
| spectral deflation with warm-started subspace tracking (Sun, Kang, Yang 2026) | lower KKT residuals for the same projection time on four large SDPs in fp16 [primary] | already in `dnn.C` (tracked subspace) |
| single-precision eigendecompositions in ADMM (Sinjorgo et al. 2025) | 600 iterations in 114 s against 205 s, bound difference at most 0.045 [primary] | already used (fp32/bf16) [repo] |
| ADMM theory (Kang, Jiang, Yang 2025) | local linear convergence under strict complementarity alone; slow cases are near-violations of it [primary, abstract] | diagnostic only: here rank R* (about 1,500) plus dual rank (about 600) is close to K = 2,117 [repo], the favourable case |
| warm start in SDP branch and bound | MADAM: within a node only, "passing cutting planes and the output of the bounding routine from parent to children" is future work; BiqBin: bundle, no inheritance [primary]. SOS-SDP (Piccialli, Sudoso, Wiegele 2022): children inherit the parent's last-round cuts, cutting-plane rounds capped at 50 at the root and 30 below [primary]. Chalkis et al. 2023: a warm-startable projection method for QUBO bounds [primary, abstract] | we already inherit Y, Z and the pool; nobody reports more |

GPU arithmetic [e]: 22 GEMMs of order 2,117 in fp16 are about 4e11 flops, a few
milliseconds on a current GPU if the small matrices are batched, against 88-330 ms per
iteration on a CPU core. This agrees with the roadmap's 3-8 ms and is the only large
solver-side lever the literature supports.

| source | id | verified |
|---|---|---|
| Hou, Tang, Toh, RNNAL / RiNNAL+ | arXiv 2502.13849 (section 6.1, Table 8), 2507.13776 (Tables 11, 12) | [primary] |
| Tang, Toh, SDPF+ | arXiv 2605.17089 | [primary] abstract |
| Aguirre et al., cuHALLaR; Han et al., cuLoRADS | arXiv 2505.13719, 2407.15049 | [primary] abstracts |
| Lian, Shaikh, Hou, Xie, Toh, McADMM | arXiv 2610.01503 | [primary] abstract |
| Ding, Yang; Prakhya, Wang; Jiang (PDHG for SDP) | arXiv 2611.03638, 2610.37998, 2607.08035 | [primary] abstracts |
| Kang, Han, Groudiev, Yang; Sun, Kang, Yang | arXiv 2507.09165, 2609.21102 | [primary] |
| Kang, Jiang, Yang | arXiv 2503.20142 | [primary] abstract |
| Piccialli, Sudoso, Wiegele, SOS-SDP | arXiv 2104.11542 | [primary] section 4 |
| Chalkis, Kleinert, Sofranac | arXiv 2312.15026 | [primary] abstract |

## Q5. Branch-and-bound mechanics

**Short answer: little, and one item worth a test.**

- **Branching rules (Fujii et al., section 5)** [primary]. All three are polytomic (one
  facility to every free location, or one location to every free facility) and pick the
  line whose children have the largest mean score. Rule M: the score is the child's
  objective at the uniform point x_ij = 1/(n-r-1), no bound needed. Rule P: the
  objective of the parent's approximate DNN solution projected onto the child's linear
  manifold. Rule D: a child lower bound from the parent's dual restricted to the child.
  The estimates in the Q1 table differ by 19,000x on tho40 (P best), 25x and 3,800x on
  tai50b (M best) and 8x on sko49 (M best). Rules P and D are relatives of our rules 8/9
  and of the child bounds from the certificate; nothing here resembles M. See idea 2.
- **Bound inheritance.** Rule D is the published form of "child bounds from the
  certificate", closed here (flat reduced costs) [repo]. Their estimate for it is the
  worst of the three on tai50b and tai40b.
- **Polytomic against binary.** Every exact QAP code since 2000 branches polytomically.
  Binary branching appears only after converting to a binary problem: the tai256c work
  rewrites the instance as a cardinality-constrained BQOP and uses orbital branching
  with isomorphism pruning (LB gap 1.48% to 1.25%, not solved) [primary]. Nothing
  measured for a general QAP.
- **Group fathoming of children.** Not in the literature. The nearest device is BiqBin's
  gate: measure at the root the gain `diff` of the expensive bound over the cheap one,
  and at other nodes run the cut phase only if the cheap bound is within `diff + 1` of
  pruning; "more B&B nodes than necessary, but the performance ... improves" [primary].
  Our near-miss rule is the same device; the constant could be the measured gain of the
  cut phase per depth (on tho40 today: +468 to +636 at depth 4).
- **Ramp-up.** Rule M's independence of bounds lets them enumerate the tree to a fixed
  depth in advance (50 x 49 x 48 subproblems for n = 50) and hand those out [primary];
  our coordinator's frontier does the same job.
- **Cut inheritance.** SOS-SDP passes the parent's last-round inequalities to the
  children and lowers the round cap from 50 to 30 [primary]; `pool.on` is the same.

## Ideas worth an experiment (ranked; none is on the closed list)

1. **tai50b: cluster regions in the cut phase, in the partition literature's order.**
   (The `regions` agent is on it; this is what the literature adds.) Regions: the five
   clusters, their unions, and prefixes along the first MDS axis. Order of families:
   (i) the Boolean triangles on (p, W), which de Meijer et al. find carry most of the
   gain; (ii) only if (i) pays, count ("box") inequalities on the two small clusters
   (3 and 4 locations), where the count of a facility set inside the cluster is an
   integer in 0..3 or 0..4, the small-count regime where boxes closed nug12's gap
   [repo]; (iii) the capacitated pigeonhole inequalities across the five clusters last
   (weaker than triangles in the k-partition experiments). Separate with the heavy flow
   pairs first (top 10% of entries = 61% of the weight).
   *Expected:* on sparse 2-equipartition instances the same cuts close 51-95% of the
   DNN gap [primary]; RLT2 closes 56% of the DNN gap on tai30b [primary]; tai50b has
   94.6% of its cost at the partition level [m] and a 6.0% gap. A third to a half of the
   gap at the root is plausible [e], with no evidence yet on what survives at depth 6-9.
   *Cheapest tests:* (a) the Nov 4 triple gate (`wil50_feasibility/new_cuts/literature/small/lmp_scale.py`,
   exact LP separation on the heaviest-flow facility triples) on a tai30b root iterate:
   on wil50 it found 0.027, on nug30 0.28; a value near nug's says the 3-wise slack is
   there. (b) tai30b and tai35b roots, then one tai50b depth-3 node, `dnncuts -regions`
   with `top:5` unions; kill if less than 20% of the gap to the cutoff is closed at twice
   the plain iterations.
2. **Rule M as a branching rule on tai50b and sko49.** Static, a few lines in `bnb.C`:
   score each child by its objective at the uniform point, take the row or column with
   the largest mean. Fujii et al.'s estimates favour it by 25x (tai50b) and 8x (sko49)
   over their DNN-solution rule and disfavour it by 19,000x on tho40, so it must be
   per-instance.
   *Expected:* anywhere from none to an order of magnitude [e]; their rule P is not our
   rule 9 and their estimator is good to about 3x.
   *Cheapest test:* the existing Knuth dives (`estimates/tai50b/run_dives.sh`) with
   rule M at the DNN depths, one night; for sko49, stage A's estimate with rule M.
   Kill if the pooled estimate is not at least 2x lower.
3. **Enter the cut phase early where the parent says the node will miss.** BiqBin's gate
   read in reverse: when the parent gave up by more than the children can recover in a
   plain solve, start the CP phase from the inherited pool after a short plain solve
   instead of a full one. Today's gate-test node needed 300 plain iterations and then
   120 in the phase.
   *Expected:* up to half of the iterations on nodes that reach the phase [e]; the
   tho40 agents can size it from their logs.
   *Cheapest test:* the gate-test node with `-I 100` and `-I 50` in place of `-I 300`;
   total iterations to the certified fathom.
4. **Root records.** One whole-instance run with the cut phase per instance, to stall:
   tho40 against 226,522, sko49 against 22,654, wil50 against 48,126, tai50b (plain or
   with cluster regions) against 431,245,612. No effect on the tree; it answers "is it a
   record" with a certified number and is the natural first table of a paper.
   *Cost:* a few hours of one core each.
5. **tai50b: branch on clusters before locations.** A child "facility f in cluster c"
   (5 children) is a face of the DNN (zeros on the rows of f outside c) and replaces up
   to 16 location children whose distances differ by 71 on a scale of 1,070 [m]. There
   is no QAP precedent (the BQOP conversion of tai256c is the nearest), and region
   branching was closed on wil50 on a far miss [repo], so this is a regime argument
   only: level 9 holds 64.6% of the pooled dive CPU, and in the largest dive all of its
   level-9 nodes fathom [repo, `estimates/README.md`].
   *Expected:* unknown; bounded by how much the union bound loses against the weakest
   location child.
   *Cheapest test:* at sampled level-8 parents whose children all fathomed, the DNN
   bound of the five cluster faces (needs forbidden assignments in `dnntest`, in a
   worktree). Kill if fewer than half of the cluster faces fathom.
6. **tho40: re-gate the 5-point and order-7 families after the triangle phase.** They
   were closed on wil50 (90% dense flows). Gaar-Rendl's orders 5 and 7 close all of what
   triangles leave at 10% density, 43-78% at 50% and 30-47% at 90% [primary, max-cut];
   tho40 is 40% dense.
   Only worth doing on nodes the CP phase fails to fathom.
   *Expected:* 5-15% of the residual after triangles [e].
   *Cheapest test:* the existing 5-point separation on the CP-phase iterate of three
   such nodes; kill if the largest violation is below 0.05 in triangle units.

Not ideas, but confirmed: keep the PRSM; the GPU port is the solver-side lever; no
relaxation above the DNN is affordable; bisection-specific knapsack inequalities are
not worth building.

## Files

- `README.md`: this report.
- `structure.py`, `structure.log`: sparsity and tails of the open instances' matrices,
  tai50b's clusters and the split of its best known objective [m].
- Downloaded PDFs and their text were kept in the session scratchpad, not here.
