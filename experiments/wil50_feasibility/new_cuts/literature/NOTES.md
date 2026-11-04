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
