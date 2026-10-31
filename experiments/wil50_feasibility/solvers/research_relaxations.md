# Cheaper relaxations of QAP: is anything within 0.1% of the DNN bound at m = 45-48?

Research note, Oct 31 2026. Question: is there a relaxation or reformulation whose
bound on wil50/sko49/tho40 nodes is within about 60 units (0.1%) of the facially
reduced DNN bound (ZKRW, lifted order (m-1)^2+1) but 5-20x cheaper per node?
Sources: literature from memory plus web-search snippets; every paper host
(arxiv, optimization-online, uwaterloo, opus4.kobv.de, QAPLIB, coral.lehigh,
plato.asu.edu, Tilburg) is blocked for direct fetch in this session, so values
marked (u) are from memory or inferred from snippets and should be checked.

Short answer: no. Every relaxation that is cheaper than the lifted DNN is weaker
by several percent of the objective, i.e. 20-100x the 60-unit budget. The DNN
relaxation (ZKRW R3 + Y >= 0, equivalently Povh-Rendl, Burer-Vandenbussche,
Kim-Kojima-Toh's Lagrangian DNN) is also the source of QAPLIB's current best
lower bounds on exactly these instances (sko49 22,650, wil50 48,121, tho40 226,490
in the QAPLIB mirror), so nothing published is stronger at n = 50 either.

## Reference values

- Oliveira, Wolkowicz, Xu (MPC 2018, ADMM on the same facially reduced DNN),
  results table (headers inferred from the snippet): sko49: previous best
  (Rendl-Sotirov bundle R3) 22,004 -> ADMM-DNN 22,653 (gap 3.13%); wil50:
  47,098 -> 48,125 (1.42%). These match our certified root values (wil50 1.4%,
  sko49 3.1%) and confirm that plain SDP R3 without Y >= 0, solved by the
  bundle method, is about 1,000 (2.1%) below DNN on wil50 and 650 (2.8%) on sko49.
- nug30 (our Table 4 and QAPLIB): GLB 4,539, PB 5,266, QPB 5,365, SDPB1 5,413,
  bundle R3 (SDP, no nonnegativity) 5,803 (u), DNN 5,944-5,950, OPT 6,124.
- Our wil50 nodes: root DNN about 48,130; depth-3 near misses converge to
  48,731-48,771 (43-84 short of 48,814); QPB at the same nodes 46,928-47,081
  (about 1,870 short, tuning/README.md section 1).

The scale to keep in mind: removing only the nonnegativity from the DNN
(SDP R3) already loses 2% on wil50. Every cheaper relaxation below is weaker
than SDP R3.

## The relaxations

1. Ding-Wolkowicz matrix-lifting MSDR1-3 (MOR 2009). Relaxes the orthogonal
   similarity set {X A X^T} by majorization; O(n^2) variables, PSD blocks of order
   about 2n, plus nonnegativity in MSDR3. Provably at least PB; the authors
   describe MSDR3 as comparable to SDR1 (= SDPB1) and QPB, sometimes below QPB.
   On nug30 that is roughly 5,365-5,450 (u) against DNN 5,950: about 9% of the
   objective, 500 units short on a 6,124 objective. Cost: seconds at n = 50.
   Certified bounds from the dual; fixed assignments and linear terms fit (it
   is a relaxation of the full QAP with a linear term). Verdict: 50x too weak.

2. Matrix splitting SDRMS (Peng, Mittelmann, Li, MPC 2010; Mittelmann-Peng,
   SIOPT 2010 for Hamming/Manhattan; Peng, Zhu, Luo, Toh, COAP 2015, SDRMS-SUM).
   Split B = B1 - B2 with B1, B2 PSD and relax X B1 X^T, X B2 X^T separately;
   O(n^2) variables, PSD matrices of order n to 2n. The Manhattan version uses
   the fact that a grid's Manhattan matrix is the sum of two path metrics
   (conditionally negative definite), so the splitting is natural for sko, wil,
   tho. Dominates PB, comparable to or a little better than MSDR3, dominated by
   SDP R3. Their contribution was size: first bounds for n up to 200 (tho150,
   sko100, wil100), and SDRMS held 18 of 33 best-known bounds on the open
   instances around 2010-2015 (Mittelmann's statement, snippet), all on
   instances too large for the lifted SDP. On sko49/wil50 the bundle R3
   values (22,004, 47,098) were already better (u). Gap on grid instances of
   this size roughly 4-8% (u). Cost: seconds to a minute at n = 50. Certified
   via the dual; works with fixed assignments. Verdict: 30-80x too weak.

3. de Klerk-Sotirov-Truetsch (IJOC 2015), Truetsch's thesis (Tilburg 2014).
   A relaxation with matrix variables of order n that dominates both MSDR3 and
   SDRMS and is dominated by ZKRW (u: I recall that ordering being proved);
   symmetry reduction applies; two new best-known QAPLIB bounds (on
   esc-type instances with large automorphism groups, u). Truetsch's SDP-based
   branch and bound reduced node counts against GLB/PB-based trees, not against
   DNN-based ones. I could not recover its nug30 or tho30 value; from the
   dominance chain it lies between QPB and bundle R3, so a gap of 5-10% at
   n = 30 (u). Verdict: the best of the O(n)-matrix family, still an order of
   magnitude too weak for this purpose.

4. Eigenvalue and grid-structure bounds. PB (Hadley-Rendl-Wolkowicz) and
   QPB (Anstreicher-Brixius) are already in the code and 1,870 short at wil50
   depth 3. Karisch-Rendl triangle decomposition bound (MP 1995) and the
   Clausen-Karisch-Perregaard-Rendl rectilinear B&B (1998): the strongest
   LP/decomposition bounds for grids in the 1990s, still weaker than
   SDPB1-class bounds at n = 30 (u). Gueye-Michelon (EJCO 2015) Manhattan LP
   with O(n^2) distance variables improved seven best-known bounds on large
   Manhattan instances in 2015, all superseded by Mittelmann's BBCPOP
   (Lagrangian-DNN) values in 2018. Nothing here is within reach of SDP R3.
   Verdict: out.

5. Lagrangian-DNN (Kim-Kojima-Toh, MP 2016; BBCPOP; Newton-bracketing,
   Kim-Kojima-Toh 2020; used by Fujii et al., ZIB 21-20). Same ZKRW-type lifted
   DNN cone (order n^2+1), with the assignment/gangster equalities moved into a
   Lagrangian penalty and the projection onto PSD intersect nonneg done by
   accelerated projected gradient; the bracketing yields a valid lower bound
   from an inexact solution (that is its selling point, like our certificate).
   Strength: the same relaxation, and their wil50 root gap and tree estimates
   (2.1e7-3.9e7 nodes) agree with ours. Per iteration: one eigendecomposition
   of the n^2 matrix, the same cost driver as ours. Verdict: equal strength,
   no 5-20x saving; relevant only as a solver alternative (other agent's topic).

6. 2015-2026 reports on wil50/sko49/tho40/tai. Oliveira-Wolkowicz-Xu 2018
   (ADMM, values above); Graham-Hu-Im-Li-Wolkowicz 2022 (restricted dual PRSM,
   same relaxation, we already use PRSM); Mittelmann 2018 BBCPOP bounds (same
   relaxation, large instances); Fujii et al. 2021 (same relaxation, solved
   tai30a and sko42); Kim-Kojima-Mittelmann-Shinano 2022/2024 (tai256c as a
   256-variable BQOP; a weaker, symmetry-specific relaxation); Tang-Toh 2025
   RiNNAL+ (low-rank Riemannian ALM for SDP-RLT/DNN relaxations: solver speed,
   same relaxation); the 2024 "error bounds for rank-one DNN reformulations"
   paper (penalty heuristics, not bounds). No paper computes a relaxation
   stronger than DNN at n >= 40, and no paper reports a cheaper one within 1%
   of it.

## Assessment

The gaps between distinct QAP relaxations are 2-10% of the objective (QPB ->
SDPB1 -> SDP R3 -> DNN: 12.4%, 11.6%, 5.2%, 2.8% on nug30). The budget here is
0.1%, and losing even the nonnegativity block costs 2% on wil50. So no known
relaxation is within 60 units of the DNN on wil50 nodes at any price below the
DNN's own; the low-dimensional family (MSDR3, SDRMS, de Klerk-Sotirov-Truetsch)
is 100-1,000x cheaper per iteration and 30-80x too weak. The only objects with
DNN strength are reformulations of the same cone (Lagrangian-DNN, Povh-Rendl,
Burer-Vandenbussche), whose cost is the same O(n^6) eigendecomposition.

Ranked candidates (none is expected to pass; ranked by chance of staying within 60):

1. Lagrangian-DNN / Newton-bracketing: passes the strength test by construction,
   fails the cost test. Worth nothing new as a relaxation; its only use is as
   an independent check of our certificates.
2. A block-sparse DNN: keep Y >= 0 everywhere but impose PSD only on a chordal
   or clique cover of principal submatrices of the lifted matrix (so the
   projection costs k small eigendecompositions). No literature for QAP; the
   extreme case (diagonal blocks only) degenerates to the Ding-Wolkowicz level,
   so the loss is probably large, but it is the only "cheaper DNN" whose loss is
   unmeasured and testable in our ADMM with a changed projection step.
3. SDRMS-SUM (Manhattan splitting) or the de Klerk-Sotirov-Truetsch relaxation:
   cheap, certified, handles fixed assignments; expected 1,000-2,000 short at
   the wil50 root and the depth-3 nodes.

## Cheap first experiments

A. Upper-bound every cheaper candidate at once with the existing code: run the
   near-miss depth-3 nodes (6, 9, 14, 15, 16 of depth3/nodes.txt) with the
   nonnegativity of Y switched off (SDP R3 only; add a flag to the projection if
   none exists). If SDP R3 is more than 60 below the DNN value on these nodes,
   which the nug30 and wil50 root figures predict by a factor of 10-20, then
   every relaxation dominated by SDP R3 (all of items 1-4) is excluded without
   implementing any of them. Cost: five runs, under an hour.

B. If A surprises, implement SDRMS-SUM or MSDR3 (PSD of order 2m = 94, O(m^2)
   variables; CVXPY with SCS or MOSEK, a few hundred lines) and compute the
   root (compare 48,125-48,133) and the five near misses (compare
   48,731-48,771). Pass criterion: within 60 of DNN on the near misses.

C. Optional, the one untested idea: in dnntest replace the full PSD projection
   by projections onto the (m-1) principal blocks indexed by facility pairs
   sharing a facility (a clique cover), run the same five nodes and record the
   loss. Any loss above 60 ends it.

My expectation: A excludes everything in one hour, and the conclusion of
NOTES.md stands: at depths 3-5 of wil50 the relaxation is the limit, and no
cheaper relaxation is close to it.
