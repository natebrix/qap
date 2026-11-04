# Facility side and disjunctive strengthening: working notes (agent 3, core 3)

Started Nov 4 2026, 05:27. Append-only log of what was read, found, run.

## 05:30 Read
- latest/dnn.C header, latest/dnncut.h, experiments/dnn_cuts/dnncuts.C header
- stronger_bounds/README.md steps 1-13, research_bound_quality.md (1.3, 1.4, 3.4-3.6, 5),
  roadmap/bound_strength.md, NOTES.md sections 1-5, triples/ and local7/ READMEs
- Key facts: regional triangles on the 13 grid thresholds (side 0 = first matrix =
  Manhattan grid) give +158.6 on node 13, +131 on node 10, +133 on d2 after the quiet
  round (laptop dual/tri_*.log); 5-point, triples, local-7 all washes or negative.
- research_bound_quality.md section 5 lists "facility-side regions on wil50" as
  do-not-bother (flow matrix random 0-9, no decomposition identity). To be tested.
- Machine: Xeon @2.8 GHz Cascade Lake (avx512, NO AMX), 4 cores, 15 GB.

## 05:35-05:45 Setup
- Conventions checked: first index = rows of A = wil50's 5x10 Manhattan grid (index r*10+c),
  second = rows of B = flows. Here "location" = first index, "facility" = second index;
  a fix a:b = facility b on location a; lifted index of x_(a,j) = 1 + a*m + j (local).
  scripts/ynode.py reproduces the logged primal from a dumped Y to 1e-6 (plumbing node).
- Speed here (no AMX on this Cascade Lake): dnntest node 13, -d 4: 477 ms/it with DNN_AMX=1;
  -d 3 -L 2 -r 40: 389 ms/it with DNN_AMX=1, 346 ms/it with DNN_AMX=0. With the grid agent's
  run on another core, 0.56 s/it measured in the phase. Using DNN_AMX=0 -d 3 -L 2 -r 40.
- Built src/fdrv (src/build.sh): links latest/*.o read-only + src/dnncut_f.C, a renamed copy of
  latest/dnncut.C made by src/make_dnncut_f.py (adds an export hook only). Modes: -phase
  (triangle protocol, stops before the quiet round, saves the pool's Lagrangian term + ADMM
  state), -cont (continue from a saved state with extra fixed terms), -pool (Polyak on an
  external sparse cut list on top of the frozen triangle term). Plumbing (m=10 node):
  fdrv -phase + -cont reproduces dnncuts' quiet round to all printed digits (44974.9659).
- Flow matrix (logs/flow_analysis.log): statistically a uniform random symmetric 0-9 matrix.
  Spectrum top 223.8 then 34.4 (random matrix: 227.4, 31.7); centered top 39.3 vs 36.9 random
  (semicircle edge 41.3); triangles with all flows >= 9: 32 (random expectation 29), >= 8: 168
  (165); densest 4/6/10/15-subsets as dense as in a random 0-9 matrix (mean flow 9.0/8.13/7.0/6.3
  vs 8.83/8.0/6.93/6.16). No clusters, no block structure: facility "clusters" are noise.
- Literature fetches: arxiv.org and kups.ub.uni-koeln.de do not resolve from this session
  (WebFetch blocked); WebSearch snippets work. Literature is cited from search results + memory.
- 05:44 launched runs/n13 (run_phase.sh, 8 rounds; 12-round run killed: at 0.56 s/it it would
  have taken ~37 min, above the 30-min limit).

## 05:50-06:00 Diagnostics at node 13's uncut converged iterate
Iterate: the grid agent's dump new_cuts/grid/runs/n13_uncut.Y.bin (dnntest 1,500 it, -d 3 -L 2
-r 40, lb 48,627.29, primal 48,628.38), copied to runs/n13_uncut1500.Y.bin (read-only use).
scripts/ydiag.py (logs/ydiag_n13_uncut1500.log), box_sep.py (logs/box_n13_uncut1500.log),
sep_dist.py --fam hull3 (logs/hull3_n13_uncut1500.log), gapdecomp.py (logs/gapdecomp_n13_uncut1500.log):
- Facility triples: min over all 16,215 triples of Delta_jk + Delta_jl + Delta_kl = 4.79 (valid >= 4,
  grid bipartite + injective): parity cut never violated. Checkerboard region: no triangle violated
  (max -0.0067). Facility regions (random facility sets, s over locations): no triangle and no box
  inequality violated (max -0.03..-0.10).
- Threshold regions: max triangle violation 0.034-0.080, 41-2,530 violated per region (as in the record).
  General box inequalities (any support, local search) on the thresholds: best normalized violation
  (half the box measure = triangle units) 0.03-0.064 with supports 3-19, i.e. no larger than the
  triangles. Parity regions: checker 0, rowpar 0, colpar 0.023 (support 3).
- Triple distance hull (LP over the 520 achievable distance triples): violated only by metric
  triangle facets H = (1,-1,1) (max 0.49, 20 of the 2,000 smallest-sum triples), which are sums of
  the 13 thresholds' apex triangles on the same facility triple: implied by the regional family.
- Gap vs a 49,036 local-search solution (scripts/node_ls.py): the gap (407.6) is a small difference
  of +6,424 / -5,807 over pairs; heavy pairs (B 8-9) sit at mean expected distance 3.46 in the DNN vs
  3.62 in the solution, light pairs (0-2) at 6.16 vs 6.01. Per-facility shares are uncorrelated with
  flow degree (corr 0.00): no facility group carries the gap; it is the diffuse bisection smoothing
  that the threshold cuts address.
- Lines study (tuning/lines/work, 300-it cold children): node 10's best lines row 49 and col 45 share
  their lowest child 49:45 (48,625.0 at 300 it); min child per line is parent -37..+57 at 300 it.

## 06:00 More diagnostics at node 13 uncut converged (same iterate)
- Facility-region triangles with T OPTIMIZED (scripts/facreg_search.py: for each location triple
  the most violated T is a 47-variable QUBO, multi-start local search; 60 location triples x 4
  types, half of them compact): max violation 0.0062 (logs/facreg_n13_uncut1500.log). So no
  facility region at all, including any region of a cut decomposition of B + c(J - I) (which
  would make the objective linear in facility-regional moments, as the 13 thresholds do on the
  location side), has a violated triangle.
- Dispersion (scripts/disp_check.py, logs/disp_n13_uncut1500.log): smallest DNN sum of expected
  pairwise distances over k facilities vs the grid minimum sigma_k: k=4 10.74 vs 8, k=5 20.62 vs 16,
  k=6 33.33 vs 25, k=7 50.05 vs 38. Slack 2.7-12: the DNN never places a group near compactness.
  "Facility groups placed too compactly" is refuted at this iterate.
- Lines study x-weighted (logs/lines_xweighted.log): on every line of the three parents, the
  x-weighted mean of the children's 300-iteration bounds is above the parent's converged bound on
  most lines (parent 1: +66..+124; parent 3: -21..+56; parent 10: -12..+34), so the parent's
  optimum is outside the hull of the line's children (if Y* = sum x_c Y_c then V >= sum x_c V_c):
  one level of lift-and-project on a line does cut it off, but the closure's value is the MIN
  child, which is the child the DNN favours most (x = 0.30-0.47 on the best lines; node 10's best
  lines row 49 and col 45 share the bottleneck child 49:45, x = 0.47).
- Box inequalities with BOTH sides free (scripts/box2_search.py, logs/box2_n13_uncut1500.log):
  count c = |pi(Q) cap P| for any location set P and facility set Q (Juenger-Kaibel boxes,
  unsigned), alternating QUBO search from 13 threshold starts + 17 random starts: violated only
  when P is a threshold (max 0.031 triangle units, |Q| = 15 and 12); every random-start search
  ends at zero violation. Together with facreg_search: at this iterate the DNN is consistent on
  every count except those of the 13 threshold regions, i.e. it "cheats" only in the variables
  the objective is written in.
- Cross-check: the grid agent (new_cuts/grid/NOTES.md, 05:30-05:55) finds on the same iterate the
  same triple minimum 4.79, zero violated parity-region triangles, col-parity max 0.023, and
  flow-weighted sub-QAP slack 140-1,700: independent agreement.
- Facility 45, the xbar line at all three lines-study parents (tuning/README.md calls it
  "location 45"; in the lines code "col" is the second index, i.e. a FACILITY), is the LIGHTEST
  facility of wil50 (flow row sum 178, rank 50 of 50; mean 222). The DNN sends it to the
  periphery: at node 10 its largest marginal is location 49 (the corner (4,9)), x = 0.47, which is
  also the bottleneck child of the best lines. Facility-side reading of the branching data: the
  lightest facility's line is the most decided one; rule 8/9 (xbar) already finds it.
- Node 13 uncut: corr(flow degree, DNN expected eccentricity of the facility's location) = -0.90,
  corr(flow degree, sum_a x_aj^2) = -0.48: heavy facilities central, light ones peripheral and more
  decided. This first-order structure is the GL/star part, already exact in the DNN.

## 06:07 Note on the deliverable
- Writing REPORT.md with the Write tool is refused by the harness for subagents ("return findings
  as text"). The report will be delivered as the text of the final hand-back to the coordinator,
  who can commit it as REPORT.md. NOTES.md (this crash-recovery log) is kept as asked.

## 06:12 Queue (runq.sh, jobs.txt) and plan
- runs/n13: triangle phase (8 rounds, -d 3 -L 2 -r 40; rounds track dual/tri_13.log to +-1:
  r1 +39.4, r2 +68.9, r3 +89.2, r6 +96.9, r7 +103.6), then quiet 1,500 with -C (child bounds
  from the cut certificate), dumping quiet.Y.bin.
- runs/valid: validate_fdrv.py (extra-region path, n=8 enumeration).
- runs/n10c: the cell disjunction at node 10's bottleneck cell 49:45 (x = 0.47): child 49:45,
  parent with 49:45 eliminated (-E), parent itself (-C 2 -X), 1,500 iterations each. One level of
  Sherali-Adams/N+ on the variable x_(49,45) has value min(child, eliminated); every disjunction on
  row 49 or column 45 (incl. the full line closures) is bounded by the child.
- next: count-split disjunction on node 13 on top of the triangle state (scripts/pick_split.py,
  count_term.py: both sides by Lagrangian terms at a few multipliers from the quiet state);
  runs/d2 phase + quiet for a second node's triangle-phase iterate.
- Star inequalities (scripts/star_check.py, logs/star_n13_uncut1500.log): hub j with weights
  B_j.: DNN value minus the conditional GL bound sum_a x_aj LAP_a is 110-264 per hub (minus
  min_a LAP_a: 144+). Implied (proof in the report) and far from binding.
