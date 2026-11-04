# literature/ working notes (append-only)

Agent: literature-driven search for relaxations / inequality families beyond
the DNN for wil50. Core 1 only. Deliverable: REPORT.md.

## 2026-11-04, start

Read: latest/dnn.C top comment (certificate: any P containing the lifted
permutations works in the polyhedral step), latest/dnncut.h (Lagrangian cuts,
L_ext), stronger_bounds/README.md steps 1-13, triples/README.md,
roadmap/bound_strength.md, solvers/research_relaxations.md, NOTES.md s.5,
research_bound_quality.md s.3.4-3.6.

Already known in the repo (do not repeat): RLT2 nug30 5770 (HH dual ascent) /
5846 (Date-Nagi) vs DNN 5944; nug25 RLT2 3577 vs DNN 3662; nug20 RLT2 2508 vs
DNN 2503; RLT3 memory 10 TB at n = 50; Kim-Kojima cutting inequalities;
MADAM/BiqBin; Jünger-Kaibel box = count inequalities; Laurent degree-4 vs
5-point.

Observation on step 8 (triples): Astra's triple separator is the
FACILITY-side 3-local marginal polytope (LMP3_fac): for facilities i,k,p,
a distribution over injective (j,l,q) with the three pair blocks of Y as
marginals. RLT2 projected onto Y is contained in LMP3_fac AND LMP3_loc
(location triples) AND the coupling (one z per assignment triple shared by
its facility triple and its location triple). So the trial tested one half
of RLT2's constraints, on 60 random triples. Every regional triangle is a
valid inequality of LMP3_fac for its facility triple, so LMP3_fac on the
right triples dominates the regional triangle pool.

Plan: (1) literature, values; (2) small-instance experiment (nug12-type
Manhattan instances and wil50 subinstances): DNN vs DNN+regional triangles vs
DNN+LMP3_fac vs DNN+LMP3_loc vs DNN+RLT2 vs RLT2 alone vs OPT, exact conic
solves (cvxpy + Clarabel, installed in a scratchpad venv); (3) wil50 node test
of whatever survives.

Tools: venv at /tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv
(numpy scipy clarabel highspy cvxpy 1.9.3).

## 2026-11-04 ~06:20, literature so far (papers/ holds the PDFs, read as data)

- RLT3 (Hahn, Zhu, Guignard, Hightower, Saltzman, IJOC 24(2) 2012; preprint
  optimization-online 2008/09/2094, papers/rlt3_hahn2012) Table 2, dual-ascent
  values: nug12 577.2 (RLT2 578), nug15 1149.74 (1150), nug18 1930 (1905),
  nug20 2569.05 (2508), nug22 3594.04 (3511), nug24 3486.12 (3369), nug25
  3723.9 (3577), rou20 725314 (699390), tai20a 703482 (675870), tai25a
  1133716 (1091165), had16-20 / rou15 exact. 384 GB for n = 25; "sizes > 25
  still a challenge".
- DNN values: Oliveira-Wolkowicz-Xu (MPC 2018, arXiv 1512.05448) Table 1 col 4
  and Graham et al. (IJOC 2022, arXiv 2006.01529): nug12 568, nug15 1141/1142,
  nug18 1894, nug20 2507/2508, nug22 3529/3530, nug24 3402, nug25 3626, nug27
  5130, nug30 5950, rou20 695181, tai20a 671675, tai25a 1096657, tai30a
  1706872, tho30 143576 (OWX mislabels the kra30a/tho30 rows), tho40 226522,
  sko42 15336, sko49 22654, wil50 48126. Graham et al.'s "strengthened" = the
  redundant trace(R) = n+1 constraint used inside the splitting (better
  subproblems and dual estimates), the SAME relaxation (sec. 2.2.2).
- RLT2 at n >= 27 (Date & Nagi, arXiv 1710.03732, IJOC 2019, Table 5, 2000
  dual-ascent iterations): nug27 5076, nug30 5846, tai30a 1687500, tai35a
  2178360, tai40a 2777580, tho30 140827, tho40 213372, sko42 14741. All
  below the DNN; their Table 6 says Burer-Vandenbussche L&P (= DNN class)
  dominates RLT2-DA for n > 22.
- RLT3 closes 83-100% of the DNN gap on nug18-25, rou20, tai20a; 52% on
  tai25a. RLT2 closes -42%..+31% at n >= 18 (100% at nug12/15).
- de Meijer-Sotirov CG cuts (MP 2024): no QAP computations. Huber-Riedl
  (arXiv 1710.02472): LP linearizations only.

Small-instance experiment (small/relax_small.py, cvxpy + SCS/Clarabel/HiGHS):
- Clarabel is too slow with the cubic variables (KKT fill-in); SCS works
  (n = 7 DNN+RLT2 in 61 s).
- DNN is EXACT on every depth-2 and depth-3 node of nug12 tried (10 + 12
  random nodes, logs/scan_nug12_k{2,3}.log) and on wil50 sub-instances
  2x3/2x4 (gap <= 0.45). Gaps only at the nug12 root (568 vs 578).
- Launched run_nug12.sh (background, core 1): nug12 root, grid on side 0:
  dnn, rlt1lp, rlt2lp, tri, lmp1, lmp0, rlt2, lmp01; logs/nug12_<v>.log.

## ~07:00 more literature, first gates

- Jünger-Kaibel box inequalities (MP 91, 2001; report zpr97-285,
  papers/kaibel_box_zpr97-285, Tables 1-3): LP (symmetric equation bound) +
  1-box cutting planes: nug12 576.3 (DNN 568, opt 578), nug15 1129.4 (DNN
  1141), nug16a 1567 (1600), nug16b 1208.2 (1219), nug17 1643.5 (1708),
  nug18 1809.3 (1894), nug20 2313.5 (2507); all esc16 except esc16a closed to
  the optimum (esc16a 66 of 68; DNN 64). 1-box beat 2-box and 4-box in their
  separation; b in {2, 3}. So at n = 12 the LP + boxes beats the DNN; at
  n = 20 it is 200 below. Nobody has run boxes on top of an SDP/DNN for QAP.
- de Klerk, Nagy, Sotirov, Truetsch, EJOR 2014 (optimization-online 2012/01/
  3303, papers/deklerk_symrlt_3303): QAP_SDP+RLT-2 (every slice Z[ij] DNN,
  RLT2 marginals) >= Povh-Rendl; computable only with large symmetry
  (Higman-Sims k-section: min 20-section 950 -> 975 (ub 980), min 25-section
  960 -> 1000 (opt), max 4-section 1100 -> 1048). n^2 PSD blocks of order n^2.
- RLT1 LP on nug12 computed here: 522.894 (HiGHS) = Resende et al.'s 523:
  the formulation checks out.
- Box gate on wil50 node 13 (wil50/box_gate.py, local search per b over Q,
  13 threshold regions), raw violation = -E[(q-b)(q-b-1)] (= 2x the
  dnncut.h normalisation): triangle-capped iterate (grid agent's
  runs/n13_tri5.Y.bin, 5 triangle rounds + 600 quiet) max 0.152, top-20
  0.12-0.15, from |Q| = 4..20, b = 1..7; uncut iterate (runs/n13_uncut.Y.bin)
  max 0.150. So ~0.07 in triangle units, the same scale as the 5-point
  family's violations; the triangles barely reduce them. Same shape as the
  5-point and local-7 gates that passed but washed out as pools.

## ~07:30 note on the deliverable

The harness refuses to let this subagent Write a file named REPORT.md
("subagents should return findings as text"). The report will be delivered
as the agent's final message to the coordinator; the coordinator can save
it as REPORT.md. NOTES.md holds everything needed to reconstruct it.

Gap-unit reading of the literature tables ((bound - DNN)/(OPT - DNN)):
  RLT2 (best of Hahn / Date-Nagi): nug12 +1.00, nug15 +1.00, nug18 +0.31,
  nug20 +0.21, nug22 +0.42, nug24 -0.38, nug25 -0.14, nug27 -0.52, nug30
  -0.60, rou20 +0.14, tai20a +0.13, tai25a -0.08, tai30a -0.17, tai35a -0.19,
  tai40a -0.22, tho30 -0.43, tho40 -0.94, sko42 -1.25.
  RLT3 (Hahn et al. 2012): nug12 0.92 (verifies opt), nug15 0.97, nug18 1.00,
  nug20 0.98, nug22 0.97, nug24 0.98, nug25 0.83, rou20 0.99, tai20a 1.00,
  tai25a 0.52.
  JK LP + 1-box: nug12 +0.83, nug15 -1.29, nug16a -3.3, nug16b -1.04,
  nug17 -2.8, nug18 -2.4, nug20 -3.1; esc16 closed where the DNN is not.
Our regional triangles on the wil50 test nodes, in units of the gap to the
cutoff (a lower bound on the node's true DNN gap): node 10 +88/158 = 0.56,
node 13 +123/182 = 0.68 (slow30: +0.72), d2 +85/109 = 0.78 (slow30 0.83).
So DNN + regional triangles at m = 46-47 already closes a fraction of the
gap comparable to RLT3 at n = 25 (0.52-0.83), and far more than RLT2 does at
any n >= 20.

## ~06:45 wil50 gates (node 13, grid agent's triangle iterate runs/n13_tri5.Y.bin)

- Location-triple LMP3 (the RLT2 half step 8 did not test; wil50/lmp_loc_gate.py,
  exact LP over the 97,290 injective facility triples per location triple,
  potentials in [-1,1]): central adjacent-cell triples (L-shapes, segments)
  max 0.044 (cells 14,24,34), median 0.017; random triples max 0.016, median
  0.011 (wil50/n13_tri5.lmploc.log). The facility agent measured facility
  triples at 0.02-0.027 (heaviest) / 0.002 (random). So BOTH halves of
  projected RLT2 are consistent to 1-4% of probability mass at the triangle
  iterate. Each LP 10-20 s.
- Box (JK 1-box) on grid rectangles instead of threshold regions: max 0.057
  raw (wil50/n13_tri5.boxrect.log) vs 0.15 on the threshold regions: the
  DNN's count errors sit on the threshold regions only.
- Grid agent (their NOTES): at the triangle iterate every single threshold is
  at or above its own bisection optimum in the objective direction; what is
  left is row-column coupling (mid5 super-cell coupled optimum 117 above the
  iterate). => every SINGLE-REGION family (5-point, local-7, boxes with P a
  threshold region, the per-region M2/Lasserre block) is bounded by a
  quantity that is already met; the literature mechanism for the coupling is
  RLT2's location-side (conditional capacity) equations, at a coarse level.
- Queue d (small/run_queue_d.sh) runs the nug12 cutting planes, then
  wil50/run_box.sh d2 "1 4" 50 (box cuts at fixed mu on top of the triangle
  phase, 6 rounds + quiet 900), then the LP decompositions.

## ~06:50 box plumbing/validity check (wil50/plumb/)

nug12 root (grid on side 0, small/inst/nug12_sw.dat), dnntest 600 iterations,
cutoff 578: uncut 567.98; box_gate over the 5 threshold regions finds 1-box
violations up to 0.193 raw (no triangle is violated at this root); the top 10
cuts at fixed mu through -Ladd: mu 0.5 568.81 (+0.83), mu 2 570.29 (+2.31),
mu 8 490.4 (overshoot). Never above 578. make_box_ladd.py checks every cut on
3,000 random permutations plus one permutation per count value (all tight).

## ~07:05 mixed triples (coupled projection of RLT2)

Observation: every RLT2 equation is a marginal equation of some MIXED triple
(two facilities + one location, or two locations + one facility): for z with
facilities {i,k,p} at {j,l,q}, the mixed triple (i,k; q) has the marginals
sum_p z = y_(ij),(kl) (location-side RLT2 equation) and sum_l z = y_(ij),(pq),
sum_j z = y_(kl),(pq) (facility-side equations), all on the same z. So the
mixed 3-local polytope is a coupled projection of RLT2 that neither the step-8
triples (facility triples) nor location triples see. Implemented:
small/cp_small.py mix1/mix0 (validated on w8_1: reaches the optimum 332, never
above), wil50/mix_gate.py.
wil50 node 13 triangle iterate (wil50/n13_tri5.mix.log): type F (two
facilities + location) heaviest-flow pairs at central locations max 0.0185,
median 0.006; random 0.009/0.005; type L (two central locations + facility)
max 0.008. Smaller than either pure half. 10-40 s per LP.
=> On wil50 at the triangle iterate, every local projection of RLT2 we can
write (facility triples, location triples, mixed triples) is satisfied to
within 0.5-4% of probability mass.
nug12 root: facility-triple (lmp1) max violation 0.278 at the DNN optimum
(logs/nug12_cp_lmp1.log round 0).

## ~07:10 facts for the write-up

- Date & Nagi Table 5: F2-RLT2-DA, 2000 iterations; sko42 needs 1,482,642 LAPs
  per iteration on K = 95 GPUs at 21.5 s/iteration (12 h); tho40, tai40a K = 71.
  At m = 47: 2,337,122 LAPs of size 45 per iteration and ~1.6e9 cubic
  multipliers. Date & Nagi Table 6: Burer-Vandenbussche L&P (DNN class)
  dominates RLT2-DA for n > 22 except taiXXb (RLT2 0.00% gap on tai20b/25b,
  L&P 4.1%/11.7%: sparse, weakly planar, widely separated asymmetric data).
- RLT-d dual (Adams-Johnson / Hahn-Grant / Adams et al.): nested LAPs, i.e. a
  Gilmore-Lawler bound at depth d with costs split across the symmetric
  copies of each (d+1)-tuple; the LP says every depth-d conditional
  subproblem is an assignment problem, consistently.
- Roberson & Seppelt (ICALP 2023, arXiv 2302.10538): for graph isomorphism
  SA level 3t is as strong as Lasserre level t, tight; i.e. the PSD condition
  is worth a bounded number of SA levels, consistent with DNN ~ RLT2.
- nug12 lmp1 (facility-triple LMP3) cutting planes: round 0 567.99 (max viol
  0.278), round 1 571.21 with 200 cuts.

## ~07:25 nug12 DNN + facility-triple LMP3 (Astra's family), exact cutting planes

logs/nug12_cp_lmp1.log (SCS 1e-5, 200 cuts a round over all 220 facility
triples): 567.99 -> 571.21 -> 572.63 -> 573.33 -> 573.82 (rounds 0-4), max
violation 0.278 -> 0.185, still rising ~0.5 a round. So at n = 12 the
step-8 family closes >= 58% of the DNN gap ON TOP of the DNN (RLT2 alone
closes 100% here). At wil50 m = 47 the same family is violated by 0.02
(facility agent, heaviest triples, triangle iterate) against 0.28 here:
the per-triple inconsistency of the DNN shrinks with n, the same trend as
the literature's RLT2-vs-DNN gap units (+1.0 at n = 12 to -0.6 at n = 30,
-1.25 at n = 42). At m = 47 the coherent part of 3-facility consistency is
the regional aggregate, i.e. the threshold triangles we already have.
- nug12 DNN + facility-triple LMP3, round 5 (1,000 cuts): 574.18 (max viol 0.16, still rising).
- nug12 DNN + LOCATION-triple LMP3 (logs/nug12_cp_lmp0.log, 3 rounds): 567.99 ->
  569.71 -> 570.41 -> 570.86 (max viol 0.208 -> 0.133). After 3 rounds the
  facility side had 573.33: the location half is weaker at nug12.
- nug12 DNN + MIXED triples (both types, coupled projection; logs/nug12_cp_mix.log,
  300 cuts a round): 567.99 -> 570.25 -> 571.28 -> 572.00 (3 rounds, max viol
  0.238 -> 0.127), between the facility half (573.33 at 3 rounds) and the
  location half (570.86).

## ~08:55 wil50 node 13, large-Q box cuts at a fixed multiplier on the UNCUT DNN

wil50/run_box_uncut.sh: dnntest -x 21:28,8:41,28:19 -u 48816 -d 3 -L 2 -r 40
-i 1500 -H 0 (the grid agent's flags). Reference (runs/n13_unc_ref.log):
48,627.29, primal 48,628.42, 934 s. 200 box cuts separated at the uncut
iterate (13 threshold regions, |Q| 4-40, b 1-7, raw violations
0.06-0.15; all checked on 300 random permutations + one per count value),
mu = 1 (runs/n13_unc_box_mu1_k200.log): certified 48,599.43 (-27.9),
Lagrangian primal 48,638.62 (+10.2 over the reference primal). The term
(859k nonzeros, |L| 4.7e6) lags the certificate by 39 at 1,500 iterations,
exactly the triples' failure mode; the Lagrangian moves +10 at most, the
same order as local-7's +12 at its gate. mu = 3 killed (would overshoot);
mu = 0.3 queued later.

## ~09:10 nug12: DNN + Jünger-Kaibel 1-box cuts on the 5 threshold regions CLOSES THE GAP

logs/nug12_cp_box1.log (cp_small.py box1: exhaustive over facility sets Q,
|Q| >= 3, b = floor(E q), P = the 5 threshold regions of the 3x4 grid, 200
cuts a round, SCS 1e-5): 567.99 -> 576.83 -> 578.011 -> 578.011 (round 2-3,
max viol 2e-4; 578.011 vs OPT 578 is SCS tolerance). Triangles: none violated
at this root (stronger_bounds step 7). So at n = 12 the 1-box family alone
does everything RLT2 does, and more than any local RLT2 projection
(facility 0.62, mixed 0.40, location 0.29 of the gap after 3-5 rounds).
Validity: algebraic (integer count), make_box_ladd.py permutation checks,
dnntest -Ladd certified values never above 578; check_box_nug12.py re-checks
cp_small's rows on 100k permutations (running).
Consequence: the on-top-of-triangles test on wil50 is the one that matters;
queue h runs wil50/run_box_ontop.sh "0.3 1" 200 on node 13 with the grid
agent's n13_tri5 protocol (reference 48,713.85).
- check_box_nug12.py: the 200 box rows cp_small separates at the nug12 root
  are valid on 100,000 random permutations (worst value - rhs = 0) and each is
  tight on some permutation. Shapes: |P| = 6 (col<2) or 4/8, |Q| = 5-9,
  b = 2-4, i.e. b ~ E q ~ |Q|/2: large balanced boxes on the half-grid
  region, the cut-polytope clique/hypermetric inequalities of one threshold's
  bisection, not the 5-point family's |Q| = 5 alone.
- nug12 lmp1+lmp0 (decoupled both halves) round 1: 571.48 (running).
- nug12 DNN + both local halves, decoupled (logs/nug12_cp_lmp01.log, 300 cuts a
  round): 567.99 -> 571.48 -> 572.95 -> 573.73 (3 rounds; max viol still
  0.20): barely above the facility half alone at 3 rounds (573.33). All the
  local-projection runs are lower estimates (cutting planes not converged).

## ~09:50 scale of the DNN's 3-local inconsistency (small/lmp_scale.py, small/logs/nug*_lmpscale.log)

Root DNN iterates (dnntest -i 3000 -H 0): nug12 567.99, nug20 2506.34, nug30
5949.40 (QAPLIB orientation: side 0 = facilities (flow), side 1 = grid).
Exact LMP3 separation value (max / median over 30 triples):
  n   facility heaviest   facility random   location heaviest-distance   location random
  12  0.272 / 0.114       0.269 / 0.086     0.121 / 0.072                0.168 / 0.064
  20  0.398 / 0.134       0.193 / 0.020     0.055 / 0.015                0.171 / 0.026
  30  0.280 / 0.146       0.173 / 0.007     0.017 / 0.006                0.143 / 0.007
wil50 node 13 (m = 47): facility heaviest 0.027 / random 0.002 (facility
agent, uncut); location central 0.044 / 0.017, random 0.016 / 0.011, mixed
0.0185 max (mine, triangle iterate).
Random-triple medians fall ~1/n^2 (0.086 -> 0.020 -> 0.007). On nug the
heaviest-flow facility triples stay at 0.28-0.40 up to n = 30; on wil50 they
are 0.02-0.03: wil50's DNN is locally consistent even on its heaviest triples.
