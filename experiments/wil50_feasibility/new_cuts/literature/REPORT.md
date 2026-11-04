# Relaxations and inequality families beyond the DNN, for wil50

Agent 1 (literature), Nov 4, 05:26-12:05 UTC, core 1. Saved by the
coordinator from the agent's hand-back message (the agent could not write
report files). The box-cut values on nodes 13 and d2 and the nug12
cutting-plane log were re-read; note that the nug12 "DNN + 1-box" value
578.011 is an SCS objective at tolerance 1e-5 (slightly above the optimum
578), not a certified bound. `NOTES.md` is the full running log.

## Summary

1. **Only two published relaxations beat the DNN on QAPLIB instances, and
   neither is usable at n = 47.** Level-3 RLT (RLT3) closes 52-100% of the DNN
   gap at n <= 25 but needs O(n^8) memory (384 GB at n = 25, about 10 TB at n
   = 50). SDP+RLT2 has only been computed on highly symmetric instances and
   needs n^2 PSD blocks of order n^2. Level-2 RLT (RLT2) is not stronger than
   the DNN at our sizes: in units of the DNN gap it falls from +1.0 at n <= 15
   to -0.6 at nug30, -0.94 at tho40 and -1.25 at sko42. Nothing computed at n
   >= 40 beats the DNN; Graham et al.'s "strengthened DNN" is the same
   relaxation (wil50 48,126).
2. **Step 8's projected RLT2 tested one of RLT2's three local projections.**
   All three were separated exactly here: facility triples (step 8's family),
   location triples, and a new coupled one (mixed triples). On nug12, on top
   of the DNN, they recover 29-62% of the gap; on wil50 node 13 the DNN already
   satisfies all three to within 0.5-4% of probability mass, and the facility
   agent found they add 0 to +0.5 in a Polyak pool. At m = 47 the useful part
   of 3-wise consistency is its regional aggregate, the threshold-triangle
   pool we already run.
3. **The Junger-Kaibel 1-box (count) inequalities are the strongest family
   found.** On nug12 they close the DNN gap on top of the DNN (SCS, uncertified)
   using only the 5 threshold regions, where the regional triangles add 0. On
   wil50, at a fixed multiplier mu = 0.3 on top of the triangle phase, they
   raise the Lagrangian primal by +17.0 on node 13 and +16.8 on d2; certified,
   +14.3 on node 13 and -4.0 on d2 (the dense external term slows the
   certificate). Size: about +15-20, not +40.
4. **Recommendation.** One more experiment is justified: large-Q boxes as an
   in-pool cut type (about a day of work, one protocol run), killed if the
   median gain is below +15. After that, stop looking for linear families on
   Y. The remaining gap is coupling between thresholds (the grid agent's 117
   units at node 13 after the triangles); in RLT terms it is RLT2's
   location-side equations applied to coarse super-cells.

## (a) The literature's strongest relaxations, against the DNN

Gap units = (bound - DNN) / (OPT - DNN). DNN from [OWX18] Table 1 column 4
and [GHILW22]; RLT2 the better of [HZGHS12] Table 2 and [DN17] Table 5; RLT3
from [HZGHS12] Table 2; "JK" = LP + 1-box cutting planes [JK01, report
97.285, Tables 1-3]. All RLT values are dual-ascent lower bounds on the LP.

| inst | OPT | DNN (gap %) | RLT2 | gap units | RLT3 | gap units | JK LP+box | gap units |
|---|---|---|---|---|---|---|---|---|
| nug12 | 578 | 568 (1.7) | 578 | 1.00 | 577.2 | 0.92 | 576.3 | 0.83 |
| nug15 | 1150 | 1141 (0.8) | 1150 | 1.00 | 1149.7 | 0.97 | 1129.4 | -1.29 |
| nug18 | 1930 | 1894 (1.9) | 1905 | 0.31 | 1930 | 1.00 | 1809.3 | -2.35 |
| nug20 | 2570 | 2507 (2.5) | 2520 | 0.21 | 2569.1 | 0.98 | 2313.5 | -3.07 |
| nug22 | 3596 | 3529 (1.9) | 3557 | 0.42 | 3594.0 | 0.97 | - | - |
| nug24 | 3488 | 3402 (2.5) | 3369 | -0.38 | 3486.1 | 0.98 | - | - |
| nug25 | 3744 | 3626 (3.2) | 3610 | -0.14 | 3723.9 | 0.83 | - | - |
| rou20 | 725520 | 695181 (4.2) | 699390 | 0.14 | 725314 | 0.99 | - | - |
| tai20a | 703482 | 671675 (4.5) | 675870 | 0.13 | 703482 | 1.00 | - | - |
| tai25a | 1167256 | 1096657 (6.0) | 1091480 | -0.07 | 1133716 | 0.52 | - | - |
| nug27 / nug30 | 5234 / 6124 | 5130 / 5950 | 5076 / 5846 | -0.52 / -0.60 | - | - | - | - |
| tai30a / 35a / 40a | | | | -0.17 / -0.19 / -0.22 | | | | |
| tho30 / tho40 / sko42 | | | | -0.43 / -0.94 / -1.25 | | | | |
| wil50 | 48816 | 48126 (1.4) | - | | - | | - | |

- **The DNN equivalences.** The DNN is ZKRW R3 with Y >= 0, Povh-Rendl, and
  Burer-Vandenbussche's lift-and-project [BV06] (BV: nug30 5934; tai30a
  1,706,875 against 1,706,872); it also equals SDP+RLT1 [dKNST14, Thm 4].
  [GHILW22]'s "strengthened" model only adds the redundant trace(R) = n+1
  inside the splitting subproblems (their section 2.2.2).
- **RLT3** closes about 0.97-1.00 of the nug gaps at n = 18-24 and 0.83 at
  nug25, at 5.6 M CPU-seconds on a 384 GB machine for nug25; 0.52 on tai25a;
  the authors call sizes above 25 "a challenge".
- **RLT2** declines by roughly 0.1 gap units per unit of n between n = 15 and
  30 (partly the dual ascent: the RLT2 LP has never been solved exactly at n
  >= 30). At m = 47 one dual-ascent iteration needs 2.3 M LAPs of size 45 and
  about 1.6e9 cubic multipliers; [DN17] used 95 GPUs at n = 42 (21.5 s per
  iteration, 2,000 iterations). On the sparse taiXXb instances RLT2 closes
  everything where Burer-Vandenbussche leaves 4-12% [DN17, Table 6]; wil50 is
  not of that kind.
- **Box inequalities at the LP level** beat the DNN only at n = 12 and on esc16;
  at n >= 15 they are 1.3-3.1 gap units below it. No prior use of boxes on
  top of an SDP for the QAP was found.
- **Not stronger than the DNN, or not new for us:** PB and QPB (dominated);
  Mittelmann-Peng matrix splitting [MP10] (wil50 46,526); Ding-Wolkowicz and
  de Klerk-Sotirov-Truetsch low-dimensional SDPs (weaker); vertex fixing
  [dKS12] (depth-1 branching); Kim-Kojima successive cuts (child bounds);
  Dong-Anstreicher 5x5 DNN-versus-CP cuts [DA13] (atomic level; atomic
  triangles gave +0.2 on nug30); de Meijer-Sotirov Chvatal-Gomory cuts for
  integer SDPs [dMS24] (no QAP computations); Qu-Sun-Xu DNN cutting planes
  (continuous QPs); 2023-2026 work (BBCPOP/Newton-bracket bounds, RiNNAL+,
  tai256c as a BQOP) solves the same relaxation faster.
- **Theory.** For graph isomorphism, Sherali-Adams level 3t is as strong as
  Lasserre level t, and 3t cannot be lowered [RS23]: positive
  semidefiniteness is worth a bounded number of RLT levels, which fits DNN
  about RLT2 in the table.
- **What makes RLT-d strong.** Its LP requires every subproblem conditioned on
  d assignments to be an assignment problem, consistently across the
  symmetric copies of each (d+1)-tuple; its dual is a nested Gilmore-Lawler
  computation with costs redistributed among the copies [HG98, AGHH07,
  HZGHS12]. RLT3 is Gilmore-Lawler at depth 3 with optimal cost splitting: at
  n <= 25 that covers much of the search tree, at n = 47 it is out of reach.

## (b) Why the projected-RLT2 attempt failed, and what would capture RLT strength

Y is RLT2-feasible if there is a cubic tensor z >= 0 (one entry per
assignment triple) with facility-side and location-side marginal equations
equal to y_ab. Three local relaxations can be separated triple by triple by
an exact LP over the m(m-1)(m-2) injective placements (97,290 at m = 47):
facility triples (step 8's family, containing every regional triangle on
every region), location triples, and mixed triples (new: two facilities and
one location or the reverse; every RLT2 equation is a marginal of some mixed
triple, so it is the local piece of RLT2's coupling).

**Table B1. nug12 root, grid on the location side as in wil50**, DNN plus
each family by exact cutting planes (SCS 1e-5; values still rising are lower
estimates; DNN gap 10.0; `small/logs/nug12_*.log`).

| relaxation | value | gap units |
|---|---|---|
| RLT1 LP (HiGHS; Resende et al. report 523) | 522.89 | -4.5 |
| DNN (dnntest certifies 567.98) | 567.99 | 0 |
| DNN + regional triangles (none violated at this root, step 7) | 567.99 | 0 |
| DNN + location triples, 3 rounds | 570.86 | >= 0.29 |
| DNN + mixed triples, 3 rounds | 572.00 | >= 0.40 |
| DNN + facility + location triples, 3 rounds | 573.73 | >= 0.57 |
| DNN + facility triples, 5 rounds | 574.18 | >= 0.62 |
| **DNN + 1-box on the 5 threshold regions, 2 rounds** | **578.0 (SCS, uncertified)** | **1.00** |
| RLT2 LP alone (dual ascent <= LP <= OPT) | 578 | 1.00 |

**Table B2. The DNN's local inconsistency against n** (value of the exact LP
separation, potentials in [-1, 1], max / median over 30 triples, at the root
iterate; `small/lmp_scale.py`, `small/logs/nug*_lmpscale.log`):

| | facility, heaviest-flow triples | facility, random triples | location, random triples |
|---|---|---|---|
| nug12 | 0.272 / 0.114 | 0.269 / 0.086 | 0.168 / 0.064 |
| nug20 | 0.398 / 0.134 | 0.193 / 0.020 | 0.171 / 0.026 |
| nug30 | 0.280 / 0.146 | 0.173 / 0.007 | 0.143 / 0.007 |
| wil50 node 13 (m = 47) | 0.027 (facility agent, uncut) | 0.002 (facility agent, uncut) | 0.016 / 0.011 (triangle iterate) |

On node 13's triangle iterate, central adjacent-cell location triples reach
0.044 (median 0.017) and mixed triples 0.0185.

1. **Step 8 tested the right object, but at m = 47 there is little left.**
   The DNN's pair blocks are consistent to within 0.5-4% of mass on every
   kind of triple, including wil50's heaviest (0.02-0.03 against 0.28-0.40 on
   nug's heaviest); random-triple inconsistency falls roughly as 1/n^2. In a
   Polyak pool on top of the triangle term, facility and location triples add
   0 and +0.3 to +0.5 (facility agent). Step 8's fixed-multiplier dense term
   cost more (-45 to -55) than any triple could give.
2. **The useful part of 3-wise consistency at this size is the regional
   triangle pool** (every regional triangle is a facility-triple inequality,
   +130 to +160 on the hardest nodes). The exact projection adds beyond the
   triangles only at small n (nug12: +6 of 10 where the triangles give 0).
3. **"Cannot be had from Y at all" is false in principle, true in practice at
   our size.** Beyond the local pieces, RLT2's projected inequalities are
   combinations over many triples whose cubic terms cancel by symmetry, with
   support on all ~2.3 M assignment pairs; separating them means solving the
   RLT2 LP, and every published RLT2 value at n >= 24 is below the DNN.
4. **What the literature measures above the DNN is level 3** (4-wise
   conditional consistency). Our regional proxy for it (local-7) was negative
   in the pool (step 11). The grid agent finds every single threshold at or
   above its own bisection optimum after the triangles; what remains is
   coupling between thresholds through the cells' capacities (mid5: 117 at
   node 13). In RLT terms that is the location-side equations applied to
   coarse super-cells: a coarse RLT2, not a projection of the fine one.
5. **What closed nug12's gap was integrality of counts, not RLT.** The nug12
   box cuts are large and balanced (|P| = 4-8 cells, |Q| = 5-9 facilities, b
   about |Q|/2); step 7's 5-point family is only the |Q| = 5 slice. All 200
   were checked on 100,000 random permutations (worst value minus rhs 0, each
   tight on some permutation).

## (c) Ranked candidate families (on top of the converged triangle phase)

| # | family | gain at depths 3-5 | separation cost | validity | effort | cheapest test -> kill |
|---|---|---|---|---|---|---|
| 1 | large-Q 1-box count inequalities on the 13 threshold regions [JK01] | [m] fixed mu 0.3, 200 cuts: primal +17.0 / +16.8, certified +14.3 (node 13) / -4.0 (d2); mu 1 overshoots (-58); nug12 0 -> 1.00 of the gap (SCS). [e] in the pool +10 to +25 certified; chance of a +40 median about 5-10% | local search over Q per region and b: 1-4 s per iterate in Python; term +3-8% per iteration; needs the 1,500-iteration quiet round | (q-b)(q-b-1) >= 0 for the integer count q; every cut checked and tight | ~1 day: generalize cut type 4 in dnncut.C to the add/drop Q search of box_gate.py | step-7 protocol on nodes 10, 13, d2 (1,000 triangles + 200 boxes per round, quiet 1,500); kill if median gain < +15 |
| 2 | coarse threshold coupling: RLT2 location-side equations on super-cells | [m, grid agent] at most 117 (mid5, node 13); [e] 0 to +30 | LP with 1e6-1e7 coarse cubic variables per separation | RLT2 equations of the coarse problem | 1-2 weeks | 4 super-cells at the node-13 triangle iterate; kill if separation value < 0.05 |
| 3 | facility, location and mixed triple polytopes | [m] violations <= 0.044; pool +0 to +0.5 | 4-40 s LP per triple | exact LP | done | closed |
| 4 | per-region Lasserre-2 (M2) block | [e] 0 to +10 (single thresholds already saturated) | +60-200% per iteration | moment PSD | 2-3 weeks | gate as in bound_strength item 3 |
| 5 | full RLT2 via dual ascent inside the certificate | [literature] RLT2 <= DNN at n >= 24 | 2.3 M LAPs per iteration | valid at every dual point | weeks | out |
| 6 | RLT3, SDP+RLT2, Lasserre-2 on Y | [literature] RLT3 0.52-1.0 at n <= 25 | O(n^8) memory / n^2 PSD blocks of order n^2 | - | - | out |
| 7 | PB, QPB, matrix splitting, low-dimensional SDPs, vertex fixing, Kim-Kojima, 5x5 CP cuts | 0 | - | - | - | out |

## (d) Tests on wil50 nodes (logs under this directory)

- **Box gate** (`wil50/box_gate.py`): node 13 max violation 0.150 uncut and
  0.152 at the grid agent's triangle iterate (raw units, about 2x the triangle
  normalisation); grid rectangles instead of threshold regions only 0.057; top
  cuts on the edge regions (col<1-3, row<1: 5-14 cells), median |Q| 16-17,
  where the triangles are weak (`wil50/n13_{uncut,tri5}.box.log`,
  `n13_tri5.boxrect.log`).
- **Location triples** (`wil50/lmp_loc_gate.py`): central max 0.044, median
  0.017; random max 0.016. **Mixed triples** (`wil50/mix_gate.py`): max 0.0185.
- **Plumbing and validity** (`wil50/plumb/`): nug12 root through `dnntest
  -Ladd`, cutoff 578: 567.98 uncut -> 568.81 at mu 0.5 and 570.29 at mu 2,
  never above 578; `make_box_ladd.py` checks every cut on 300-3,000 random
  permutations plus one per count value.
- **Box cuts on the uncut node-13 DNN** (`wil50/runs/n13_unc_*.log`, 1,500
  iterations): reference 48,627.29; 200 cuts at mu 1: 48,599.43 (-27.9),
  primal +10.2.
- **Box cuts on top of the triangle phase** (400-iteration give-up point, 5
  Polyak rounds of 300 iterations, 1,000 triangles per round, uncapped delta
  50, quiet round of 600 from the best round):

  | node | reference: certified / primal | box mu = 0.3: certified / primal | box mu = 1 |
  |---|---|---|---|
  | 13 | 48,713.85 / 48,735.51 (grid agent's n13_tri5.log) | 48,728.10 (+14.3) / 48,752.51 (+17.0) | 48,655.94 (-57.9) / 48,739.90 (+4.4) |
  | d2 | 48,749.62 / 48,778.16 (`runs/d2_tri.log`) | 48,745.63 (-4.0) / 48,794.96 (+16.8) | - |

  The term has about 0.83-0.86 M nonzeros and widens the certificate's lag
  behind the primal after the quiet round from 28.5 to 49.3 on d2. Confound on
  node 13: the reference's quiet round started from round 3's multipliers, the
  box arm's from round 5's. The Lagrangian moves +17 on both nodes;
  certification is what differs.
- **Box family at the roots, fixed mu, certified, as a fraction of the DNN
  gap** (`small/logs/nug{20,30}_sw_*.log`): nug12 0.23, nug20 0.07, nug30
  0.02, shrinking with n like the local projections.
- **CPU:** each wil50 protocol run 1,700-1,900 CPU-s; the box term +3-8% per
  iteration; separation seconds.

## (e) Recommendation

1. Do not pursue RLT2, RLT3, SDP+RLT2, Lasserre-2 or projected-RLT2 families
   for wil50: RLT2 is not stronger than the DNN at this size in any published
   value; RLT3 and SDP+RLT2 are 1e3-1e4 times beyond the per-node budget; all
   three local projections are already satisfied at our iterates.
2. One more bound experiment: large-Q 1-box count inequalities as an in-pool
   cut type (~1 day). The only literature family strictly stronger than our
   pool on a small grid instance (nug12: triangles 0%, boxes 100% of the gap,
   uncertified), outside the triangles' territory on wil50 (edge regions, |Q|
   about 16), and with a positive Lagrangian move on both test nodes. Expect
   +10 to +25 certified with optimised multipliers and a 1,500-iteration
   quiet round: a fifth of a level, roughly 1.3-2x on the tree by the
   47^(d/130) rule and less in practice. Kill if the median gain over the
   triangle reference on nodes 10, 13 and d2 is below +15.
3. Otherwise stop adding linear families on Y. The single-region families
   (triangles, 5-point, local-7, boxes, per-region M2) are bounded by
   threshold saturation, which the triangles already reach; the residual is
   threshold coupling, whose literature form is a coarse RLT2 on super-cells.
   Bound-side effort is better spent on the search, especially cut-pool
   inheritance (bound_strength.md item 1).

## References

- [AGHH07] W. Adams, M. Guignard, P. Hahn, W. Hightower. A level-2 RLT bound for the QAP. EJOR 180(3):983-996, 2007.
- [BV06] S. Burer, D. Vandenbussche. Solving lift-and-project relaxations of binary integer programs. SIAM J. Optim. 16(3):726-750, 2006.
- [DA13] H. Dong, K. Anstreicher. Separating doubly nonnegative and completely positive matrices. Math. Prog. 137:131-153, 2013.
- [DN17] K. Date, R. Nagi. RLT2-based parallel algorithms for solving large QAPs on GPU clusters. arXiv:1710.03732; INFORMS J. Comput. 2019. Tables 5 and 6.
- [dKNST14] E. de Klerk, M. Nagy, R. Sotirov, U. Truetsch. Symmetry in RLT-type relaxations for the QAP and standard quadratic optimization problems. EJOR, 2014.
- [dKS12] E. de Klerk, R. Sotirov. Improved semidefinite bounds for QAP with suitable symmetry. Math. Prog., 2012 (from memory).
- [dMS24] F. de Meijer, R. Sotirov. The Chvatal-Gomory procedure for integer SDPs. Math. Prog. 209:323-395, 2024.
- [GHILW22] N. Graham, H. Hu, H. Im, X. Li, H. Wolkowicz. A restricted dual Peaceman-Rachford splitting method for a strengthened DNN relaxation for QAP. INFORMS J. Comput. 34(4):2125-2143, 2022.
- [HG98] P. Hahn, T. Grant. Lower bounds for the QAP based upon a dual formulation. Oper. Res. 46, 1998 (from memory).
- [HZGHS12] P. Hahn, Y.-R. Zhu, M. Guignard, W. Hightower, M. Saltzman. A level-3 RLT-based bound for the QAP. INFORMS J. Comput. 24(2):202-209, 2012. Table 2.
- [JK01] M. Junger, V. Kaibel. Box-inequalities for quadratic assignment polytopes. Math. Prog. 91:175-197, 2001 (Koln report 97.285, Tables 1-3).
- [MP10] H. Mittelmann, J. Peng. Estimating bounds for QAPs with Hamming and Manhattan distance matrices. SIAM J. Optim. 20:3408-3426, 2010.
- [OWX18] D. Oliveira, H. Wolkowicz, Y. Xu. ADMM for the SDP relaxation of the QAP. Math. Prog. Comput. 10(4):631-658, 2018.
- [RRD95] M. Resende, K. Ramakrishnan, Z. Drezner. LP lower bounds for the QAP with an interior point method. Oper. Res. 43(5):781-791, 1995.
- [RS23] D. Roberson, T. Seppelt. Lasserre hierarchy for graph isomorphism and homomorphism indistinguishability. ICALP 2023.
- [ZKRW98] Q. Zhao, S. Karisch, F. Rendl, H. Wolkowicz. Semidefinite programming relaxations for the QAP. J. Comb. Optim. 2:71-109, 1998.

## Files

Under this directory: `NOTES.md` (full running log); `small/` (small-instance
work: `relax_small.py`, `cp_small.py`, `lmp_scale.py`, `check_box_nug12.py`,
the queue scripts, `logs/`); `wil50/` (node tests: `box_gate.py`,
`make_box_ladd.py`, `run_box*.sh`, `lmp_loc_gate.py`, `mix_gate.py`, `runs/`,
`plumb/`). Downloaded preprints and their extracted text (`papers/`) are kept
locally and not committed.
